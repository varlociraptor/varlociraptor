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
) -> Result<()>
where
    Post: FnOnce(Receiver<T>) -> Result<()>,
    Post: Send,
    Pre: FnOnce(Sender<U>) -> Result<()>,
    Pre: Send,
    Workers: Iterator<Item = W> + ExactSizeIterator,
    W: FnOnce(Receiver<U>, Sender<T>) -> Result<()>,
    W: Send,
    T: Send,
    U: Send,
{
    scope(|scope| -> Result<()> {
        let (in_sender, in_receiver) = bounded(1);
        let (out_sender, out_receiver) = bounded(1);

        let preprocessor = scope.spawn(move |_| preprocessor(in_sender));

        let workers: Vec<_> = workers
            .map(|worker: W| {
                let in_receiver = in_receiver.clone();
                let out_sender = out_sender.clone();
                scope.spawn(move |_| worker(in_receiver, out_sender))
            })
            .collect();

        // drop channel clones that are not needed anymore
        drop(in_receiver);
        drop(out_sender);

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
