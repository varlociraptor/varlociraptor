use std::collections::BTreeMap;
use std::thread;

use anyhow::Result;
use crossbeam::channel::{bounded, Receiver, Sender};
use crossbeam::thread::{scope, Scope, ScopedJoinHandle};

pub(crate) enum WorkItem<T> {
    Item(T),
    Stop,
}

/// Create and execute a worker pool.
/// # Arguments
/// * `workers` - Closures that execute the work.
pub(crate) fn worker_pool<Post, Pre, Workers, W, U, T>(
    preprocessor: Pre,
    workers: Workers,
    postprocessor: Post,
    in_capacity: usize,
    out_capacity: usize,
) -> Result<()>
where
    Post: FnOnce(Box<T>) -> Result<()>,
    Post: Send,
    Pre: FnOnce(Sender<U>) -> Result<()>,
    Pre: Send,
    Workers: Iterator<Item = W>,
    W: FnOnce(Receiver<U>, Sender<Box<T>>) -> Result<()>,
    W: Send,
    T: Send + Orderable,
    U: Send,
{
    scope(|scope| -> Result<()> {
        let (in_sender, in_receiver) = bounded(in_capacity);
        let (out_sender, out_receiver) = bounded(out_capacity);

        let preprocessor = scope.spawn(move |_| {
            let ret = preprocessor(in_sender);
            // tell consuming threads that we are done
            drop(in_sender);
            ret
        });

        let workers: Vec<_> = workers
            .map(|worker: W| {
                scope.spawn(move |_| {
                    let ret = worker(in_receiver, out_sender);
                    // tell consuming threads that we are done
                    drop(out_sender);
                    ret
                })
            })
            .collect();
        let postprocessor = scope.spawn(move |_| -> Result<()> {
            let mut items = OrderedContainer::new();
            let last_index = None;

            for item in out_receiver {
                items.insert(item.index(), item);

                // Find continuous prefix, postprocess in order.
                for item in items.remove_continuous_prefix(&mut last_index) {
                    postprocessor(item)?;
                }
            }

            Ok(())
        });

        let mut errors = Vec::new();

        let ret = postprocessor.join().unwrap();
        if ret.is_err() {
            errors.push(ret);
        }

        for worker in workers {
            let ret = worker.join().unwrap();
            if ret.is_err() {
                errors.push(ret);
            }
        }

        let ret = preprocessor.join().unwrap();
        if ret.is_err() {
            errors.push(ret);
        }

        if !errors.is_empty() {
            errors[0]
        } else {
            Ok(())
        }
    })
    .unwrap()?;

    Ok(())
}

// #[derive(Getters)]
// pub(crate) struct WorkerPool<'a, U> {
//     in_sender: Sender<U>,
//     workers: Vec<ScopedJoinHandle<'a, Result<()>>>,
//     postprocessor: ScopedJoinHandle<'a, Result<()>>,
// }

// impl<'a, U> WorkerPool<'a, U>
// where
//     U: Send + 'a,
// {
//     pub(crate) fn new<P, W, WS, T>(threads: usize, in_capacity: usize, out_capacity: usize, workers: WS, postprocessor: P) -> Self
//     where
//         P: FnOnce(Box<T>) -> Result<()>,
//         P: Send,
//         WS: Iterator<Item=W>,
//         W: FnOnce(Receiver<U>, Sender<Box<T>>) -> Result<()>,
//         W: Send,
//         T: Send + 'a + Orderable,
//     {
//         let (in_sender, in_receiver) = bounded(in_capacity);
//         let (out_sender, out_receiver) = bounded(out_capacity);

//         let workers = workers.map(|worker: W| scope.spawn(move |_| worker(in_receiver, out_sender))).collect();
//         let postprocessor = scope.spawn(move |_| -> Result<()> {
//             let mut items = OrderedContainer::new();
//             let last_index = None;

//             for item in out_receiver {
//                 items.insert(item.index(), item);

//                 // Find continuous prefix, postprocess in order.
//                 for item in items.remove_continuous_prefix(&mut last_index) {
//                     postprocessor(item)?;
//                 }
//             }

//             Ok(())
//         });

//         WorkerPool {
//             in_sender,
//             workers,
//             postprocessor
//         }
//     }

//     pub(crate) fn push(&self, item: U) -> Result<()> {
//         self.in_sender.send(item).unwrap();

//         Ok(())
//     }
// }

pub(crate) trait Orderable {
    fn index(&self) -> usize;
}

struct OrderedContainer<T> {
    inner: BTreeMap<usize, Box<T>>,
}

impl<T> OrderedContainer<T>
where
    T: Orderable,
{
    fn new() -> Self {
        OrderedContainer {
            inner: BTreeMap::new(),
        }
    }

    fn insert(&self, key: usize, value: Box<T>) {
        self.inner.insert(key, value);
    }

    fn remove_continuous_prefix(&mut self, last_idx: &mut Option<usize>) -> Vec<Box<T>> {
        // TODO replace with drain_filter once stable.
        let items = Vec::new();

        for i in self.inner.keys().cloned().collect::<Vec<_>>() {
            if let Some(last_idx) = last_idx {
                if i == *last_idx + 1 {
                    items.push(self.inner.remove(&i).unwrap());
                }
            } else {
                items.push(self.inner.remove(&i).unwrap());
            }
            last_idx.replace(i);
        }

        items
    }
}
