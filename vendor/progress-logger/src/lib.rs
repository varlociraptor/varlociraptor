#[macro_use]
extern crate log;

use std::time::{Duration, Instant};
use sysinfo::{System, SystemExt};

/// A tool to report the progress of computations. It can be built and configured
/// using the `builder` function. If given the expected number of updates,
/// reports the expected time to completion, based on the current throughtput.
///
/// Progress is reported every 10 seconds by default. See the examples about how
/// to change it.
///
/// There are three methods to update the internal counter:
///
///  - `update`, for events that don't happen frequently
///  - `update_light`, which tries to report (by checking the configured frequency
///     of updates) only once every million updates. To be used in situations where
///     updates are frequent: it's an order of magnitude faster than `update`.
///
/// Reports are issued on the console using the `info!()` macro from the `log` crate.
/// Therefore, the reports depend on your logging configuration.
///
/// Inspired by `ProgressLogger` in the [`dsiutil`](http://dsiutils.di.unimi.it/docs/it/unimi/dsi/logging/ProgressLogger.html) Java library.
///
/// # Examples
///
/// ## Basic usage
/// ```
/// use progress_logger::ProgressLogger;
///
/// let mut pl = ProgressLogger::builder().start();
/// let mut cnt = 0;
/// for i in 0..10000 {
///     cnt += 1;
///     pl.update(1u32);
/// }
/// pl.stop();
/// ```
///
/// ## Reporting every 5 seconds
/// ```
/// use progress_logger::ProgressLogger;
/// use std::time::Duration;
///
/// let mut pl = ProgressLogger::builder()
///     .with_frequency(Duration::from_secs(5))
///     .start();
/// let mut cnt = 0;
/// for i in 0..10000 {
///     cnt += 1;
///     pl.update(1u32);
/// }
/// pl.stop();
///
/// ```
///
/// ## Changing the names of updates
/// ```
/// use progress_logger::ProgressLogger;
///
/// let mut pl = ProgressLogger::builder()
///     .with_items_name("points")
///     .start();
/// let mut cnt = 0;
/// for i in 0..10000 {
///     cnt += 1;
///     pl.update(1u32);
/// }
/// pl.stop();
/// ```
pub struct ProgressLogger {
    start: Instant,
    count: u64,
    expected_updates: Option<u64>,
    items: String,
    last_logged: Instant,
    /// the estimated time to completion, in seconds
    ettc: Option<f64>,
    throughput: Option<f64>,
    frequency: Duration,
    system: System,
}

impl ProgressLogger {
    /// Creates a builder to configure a new progress logger
    pub fn builder() -> ProgressLoggerBuilder {
        ProgressLoggerBuilder {
            expected_updates: None,
            items: None,
            frequency: None,
        }
    }

    fn log(&mut self) {
        let elapsed = Instant::now() - self.start;
        let throughput = self.count as f64 / elapsed.as_secs_f64();
        self.throughput.replace(throughput);
        self.system.refresh_memory();
        let used_kb = PrettyNumber::from(self.system.get_used_memory());
        let used_swap_kb = PrettyNumber::from(self.system.get_used_swap());
        if let Some(expected_updates) = self.expected_updates {
            let prediction = (expected_updates - self.count) as f64 / throughput;
            self.ettc.replace(prediction);
            info!(
                "[mem: {} kB, swap: {} kB] {:.2?} {} {}, {:.2} s left ({:.2} {}/s)",
                used_kb,
                used_swap_kb,
                elapsed,
                PrettyNumber::from(self.count),
                self.items,
                prediction,
                PrettyNumber::from(throughput),
                self.items
            );
        } else {
            info!(
                "[mem: {} kB, swap: {} kB] {:.2?} {} {} ({:.2} {}/s)",
                used_kb,
                used_swap_kb,
                elapsed,
                PrettyNumber::from(self.count),
                self.items,
                PrettyNumber::from(throughput),
                self.items
            );
        }
    }

    /// Get the estimated time to completion, if such prediction is available
    pub fn time_to_completion(&self) -> Option<Duration> {
        self.ettc.map(Duration::from_secs_f64)
    }

    pub fn throughput(&self) -> Option<f64> {
        self.throughput
    }

    /// Try to report progress only once every million updates
    #[inline]
    pub fn update_light<N: Into<u64>>(&mut self, cnt: N) {
        self.count += cnt.into();
        if self.count % 1_000_000 == 0 {
            let now = Instant::now();
            if (now - self.last_logged) > self.frequency {
                self.log();
                self.last_logged = now;
            }
        }
    }

    /// Update the internal counter and report progress if the time
    /// since the last report is greater than the configured duration
    #[inline]
    pub fn update<N: Into<u64>>(&mut self, cnt: N) {
        let cnt: u64 = cnt.into();
        self.count += cnt;
        let now = Instant::now();
        if (now - self.last_logged) > self.frequency {
            self.log();
            self.last_logged = now;
        }
    }

    /// Stops and drops the progress logger, logging the completion statement
    pub fn stop(self) {
        let elapsed = Instant::now() - self.start;
        let throughput = self.count as f64 / elapsed.as_secs_f64();
        info!(
            "Done in {:.2?}. {} {} ({:.2} {}/s)",
            elapsed,
            PrettyNumber::from(self.count),
            self.items,
            PrettyNumber::from(throughput),
            self.items
        );
    }
}

/// Builds a new progress logger. All the configurations are optional,
/// To obtain a builder, use `ProgressLogger::builder()`.
pub struct ProgressLoggerBuilder {
    expected_updates: Option<u64>,
    items: Option<String>,
    frequency: Option<Duration>,
}

impl ProgressLoggerBuilder {
    /// Configure the expected number of updates.
    pub fn with_expected_updates<N: Into<u64>>(mut self, updates: N) -> Self {
        self.expected_updates = Some(updates.into());
        self
    }
    /// Set the name of the items being counted.
    pub fn with_items_name<S: Into<String>>(mut self, name: S) -> Self {
        self.items = Some(name.into());
        self
    }
    /// Set the frequency of reports on the console.
    pub fn with_frequency(mut self, freq: Duration) -> Self {
        self.frequency = Some(freq);
        self
    }
    /// Builds the `ProgressLogger`, starting the internal timer.
    pub fn start(self) -> ProgressLogger {
        let now = Instant::now();
        ProgressLogger {
            start: now,
            count: 0,
            expected_updates: self.expected_updates,
            items: self.items.unwrap_or_else(|| "updates".to_owned()),
            last_logged: now,
            ettc: None,
            throughput: None,
            frequency: self.frequency.unwrap_or_else(|| Duration::from_secs(10)),
            system: System::default(),
        }
    }
}

struct PrettyNumber {
    rendered: String,
}

impl std::fmt::Display for PrettyNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.rendered)
    }
}

impl std::fmt::Debug for PrettyNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.rendered)
    }
}

impl From<u64> for PrettyNumber {
    fn from(n: u64) -> PrettyNumber {
        let s = format!("{}", n);
        let tmp: Vec<char> = s.chars().rev().collect();
        let mut chunks: Vec<&[char]> = tmp.chunks(3).collect();

        let mut rendered = String::new();
        let mut ul = chunks.len() % 2 == 1;
        while let Some(chunk) = chunks.pop() {
            let mut chunk = Vec::from(chunk);
            if ul {
                rendered.push_str("\x1B[0m");
            } else {
                rendered.push_str("\x1B[4m");
            }
            ul = !ul;
            while let Some(c) = chunk.pop() {
                rendered.push(c);
            }
        }
        if ul {
            rendered.push_str("\x1B[0m");
        }

        PrettyNumber { rendered }
    }
}

impl From<f64> for PrettyNumber {
    fn from(x: f64) -> PrettyNumber {
        assert!(x >= 0.0, "only positive number are supported for now");
        let s = format!("{:.2}", x);
        let mut parts = s.split(".");
        let s = parts.next().expect("missing integer part");
        let decimal = parts.next();
        let tmp: Vec<char> = s.chars().rev().collect();
        let mut chunks: Vec<&[char]> = tmp.chunks(3).collect();

        let mut rendered = String::new();
        let mut ul = chunks.len() % 2 == 1;
        while let Some(chunk) = chunks.pop() {
            let mut chunk = Vec::from(chunk);
            if ul {
                rendered.push_str("\x1B[0m");
            } else {
                rendered.push_str("\x1B[4m");
            }
            ul = !ul;
            while let Some(c) = chunk.pop() {
                rendered.push(c);
            }
        }
        if ul {
            rendered.push_str("\x1B[0m");
        }
        if let Some(decimal) = decimal {
            rendered.push('.');
            rendered.push_str(decimal);
        }

        PrettyNumber { rendered }
    }
}
