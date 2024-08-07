pub mod common;

#[macro_export]
macro_rules! testcase {
    ($name:ident, $($pairhmm_mode:ident),+) => {
        paste! {
            lazy_static! {
                static ref [<$name:upper _MUTEX>]: Mutex<()> = Mutex::new(());
            }

            $(
                #[test]
                fn [<$name _ $pairhmm_mode _mode>]() {
                    use crate::testcase::runner::common::load_testcase;
                    // Poison error can be ignored here, because it just means that the other test failed
                    // and we are safe to go on.
                    let _guard = [<$name:upper _MUTEX>].lock();
                    let name = stringify!($name);
                    let testcase = load_testcase(
                        &Path::new(file!())
                            .parent()
                            .unwrap()
                            .join("resources/testcases")
                            .join(name),
                    )
                    .expect("Failed to load testcase");
                    let mode = stringify!($pairhmm_mode);

                    // setup logger
                    // fern::Dispatch::new()
                    // .level(log::LevelFilter::Info)
                    // .chain(std::io::stderr())
                    // .apply()
                    // .unwrap();

                    testcase.run(mode).expect("Failed to run testcase");
                    testcase.check();
                }
            )*
        }
    };
}

#[macro_export]
macro_rules! testcase_should_panic {
    ($name:ident, $($pairhmm_mode:ident),+) => {
        paste! {
            lazy_static! {
                static ref [<$name:upper _MUTEX>]: Mutex<()> = Mutex::new(());
            }

            $(
                #[should_panic]
                #[test]
                fn [<$name _ $pairhmm_mode _mode>]() {
                    use crate::testcase::runner::common::load_testcase;
                    // Poison error can be ignored here, because it just means that the other test failed
                    // and we are safe to go on.
                    let _guard = [<$name:upper _MUTEX>].lock();
                    let name = stringify!($name);
                    let testcase = load_testcase(
                        &Path::new(file!())
                            .parent()
                            .unwrap()
                            .join("resources/testcases")
                            .join(name),
                    )
                    .unwrap();
                    let mode = stringify!($pairhmm_mode);
                    testcase.run(mode).unwrap();
                    testcase.check();
                }
            )*
        }
    };
}

pub use testcase;
pub use testcase_should_panic;
