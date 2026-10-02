# Contributing
Contributions are welcome!

Here are some basic guidelines:

* Open an issue for questions or bug reports. Be detailed about the platform and CPU features so the
bug can be reproduced.
* Make sure all tests pass for pull requests. Note that both the SIMD and scalar variants of the code
should return the same result. Remember to add tests for new features!
* Performance regressions due to code changes should be reported in the pull requests. New benchmarks
should be added, if necessary. Performance is very important in this library!
* Use a similar code style as the current code for pull requests.
* It may be helpful to inspect the LLVM-IR or assembly output by using `build_ir_asm.sh`.
* It may be helpful to use the debug feature flag through `--features "debug"` to get debug output.
