use std::collections::BTreeMap;

use anyhow::Result;
use crossbeam::channel::{bounded, Receiver, Sender};
use crossbeam::thread::scope;

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
    Post: FnOnce(Receiver<Box<T>>) -> Result<()>,
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
        let (buffer_sender, buffer_receiver) = bounded(out_capacity);
        let (out_sender, out_receiver) = bounded(out_capacity);

        let preprocessor = scope.spawn(move |_| preprocessor(in_sender));

        let workers: Vec<_> = workers
            .map(|worker: W| {
                let in_receiver = in_receiver.clone();
                let buffer_sender = buffer_sender.clone();
                scope.spawn(move |_| worker(in_receiver, buffer_sender))
            })
            .collect();

        // drop channel clones that are not needed anymore
        drop(in_receiver);
        drop(buffer_sender);

        let buffer_processor = scope.spawn(move |_| {
            let mut items = OrderedContainer::new();
            let mut last_index = None;

            for item in buffer_receiver {
                items.insert(item.index(), item);

                // Find continuous prefix, postprocess in order.
                for item in items.remove_continuous_prefix(&mut last_index) {
                    out_sender.send(item).unwrap();
                }
            }
        });

        let postprocessor = scope.spawn(move |_| -> Result<()> { postprocessor(out_receiver) });

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

        buffer_processor.join().unwrap();

        let ret = preprocessor.join().unwrap();
        if ret.is_err() {
            errors.push(ret);
        }

        if !errors.is_empty() {
            errors.remove(0)
        } else {
            Ok(())
        }
    })
    .unwrap()?;

    Ok(())
}

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

    fn insert(&mut self, key: usize, value: Box<T>) {
        self.inner.insert(key, value);
    }

    fn remove_continuous_prefix(&mut self, last_idx: &mut Option<usize>) -> Vec<Box<T>> {
        // TODO replace with drain_filter once stable.
        let mut items = Vec::new();

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
