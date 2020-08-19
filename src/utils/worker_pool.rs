use std::cmp;
use std::collections::BTreeMap;
use std::sync::{Arc, Condvar, Mutex};

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
    buffer_capacity: usize,
) -> Result<()>
where
    Post: FnOnce(Receiver<Box<T>>) -> Result<()>,
    Post: Send,
    Pre: FnOnce(Sender<U>, Arc<BufferGuard>) -> Result<()>,
    Pre: Send,
    Workers: Iterator<Item = W> + ExactSizeIterator,
    W: FnOnce(Receiver<U>, Sender<Box<T>>) -> Result<()>,
    W: Send,
    T: Send + BufferItem,
    U: Send,
{
    scope(|scope| -> Result<()> {
        let (in_sender, in_receiver) = bounded(1);
        let (buffer_sender, buffer_receiver) = bounded(1);
        let (out_sender, out_receiver) = bounded(1);
        let buffer_guard = Arc::new(BufferGuard::new(buffer_capacity));
        let buffer_guard_preprocessor = Arc::clone(&buffer_guard);

        let preprocessor = scope.spawn(move |_| preprocessor(in_sender, buffer_guard_preprocessor));

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
            let mut items = Buffer::new();
            let mut last_index = None;
            buffer_guard.set_size(items.len());
            let size = 0;

            for item in buffer_receiver {
                items.insert(item.index(), item);
                if !item.skip_capacity() {
                    size += 1;
                }

                // Find continuous prefix, postprocess in order.
                for item in items.remove_continuous_prefix(&mut last_index) {
                    out_sender.send(item).unwrap();
                    size -= 1;
                }
                buffer_guard.set_size(size);
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

pub(crate) trait BufferItem {
    /// Index for buffer.
    fn index(&self) -> usize;

    /// Skip item when calculating buffer capacity.
    fn skip_capacity(&self) -> bool;
}

struct Buffer<T> {
    inner: BTreeMap<usize, Box<T>>,
}

impl<T> Buffer<T>
where
    T: BufferItem,
{
    fn new() -> Self {
        Buffer {
            inner: BTreeMap::new(),
        }
    }

    fn insert(&mut self, key: usize, value: Box<T>) {
        self.inner.insert(key, value);
    }

    fn len(&self) -> usize {
        self.inner.len()
    }

    fn remove_continuous_prefix(&mut self, last_idx: &mut Option<usize>) -> Vec<Box<T>> {
        // TODO replace with drain_filter once stable.
        let mut items = Vec::new();

        for i in self.inner.keys().cloned().collect::<Vec<_>>() {
            if let Some(last_idx) = last_idx {
                if i == *last_idx + 1 {
                    items.push(self.inner.remove(&i).unwrap());
                } else {
                    break;
                }
            } else if i == 0 {
                items.push(self.inner.remove(&i).unwrap());
            } else {
                break;
            }
            last_idx.replace(i);
        }
        items
    }
}

#[derive(new)]
#[allow(clippy::mutex_atomic)]
pub(crate) struct BufferGuard {
    #[new(default)]
    size: Mutex<usize>,
    #[new(default)]
    condvar: Condvar,
    max_capacity: usize,
}

impl BufferGuard {
    #[allow(clippy::mutex_atomic)]
    pub(crate) fn wait_for_free(&self) {
        let mut size = self.size.lock().unwrap();
        while *size > self.max_capacity {
            size = self.condvar.wait(size).unwrap();
        }
    }

    #[allow(clippy::mutex_atomic)]
    pub(crate) fn set_size(&self, size: usize) {
        {
            let mut s = self.size.lock().unwrap();
            *s = size;
        }
        if size <= self.max_capacity {
            self.condvar.notify_one();
        }
    }
}
