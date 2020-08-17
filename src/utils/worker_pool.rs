use std::thread;
use std::collections::BTreeMap;

use anyhow::Result;
use crossbeam::channel::{bounded, Sender, Receiver};
use crossbeam::thread::{Scope, ScopedJoinHandle, scope};


pub(crate) fn worker_pool<P, W, WS, U, T>(threads: usize, out_capacity: usize, in_receiver: Receiver<U>, workers: WS, postprocessor: P) -> Result<()>
where
    P: FnOnce(Box<T>) -> Result<()>,
    P: Send,
    WS: Iterator<Item=W>,
    W: FnOnce(Receiver<U>, Sender<Box<T>>) -> Result<()>,
    W: Send,
    T: Send + Orderable,
    U: Send,
{
    scope(|scope| -> Result<()> {
        let (out_sender, out_receiver) = bounded(out_capacity);

        let workers = workers.map(|worker: W| scope.spawn(move |_| worker(in_receiver, out_sender))).collect();
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

        for worker in workers {
            worker.join()?;
        }
        postprocessor.join()?;

        Ok(())
    })?;

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
    inner: BTreeMap<usize, Box<T>>
}

impl<T> OrderedContainer<T>
where 
    T: Orderable
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