# triple_accel
![Test](https://github.com/Daniel-Liu-c0deb0t/triple_accel/workflows/Test/badge.svg)
![GitHub](https://img.shields.io/github/license/Daniel-Liu-c0deb0t/triple_accel)
![Crates.io](https://img.shields.io/crates/v/triple_accel)
![Docs.rs](https://docs.rs/triple_accel/badge.svg)

Rust edit distance routines accelerated using SIMD. Supports fast Hamming, Levenshtein,
restricted Damerau-Levenshtein, etc. distance calculations and string search.

Although vectorized SIMD code allows for up to 20-30x speedups over their scalar counterparts,
the difficulty of handling platform-dependent SIMD code makes SIMD routines less attractive.
The goal of this library is to provide an easy-to-use abstraction over SIMD edit distance routines
that fall back to scalar routines if the target CPU architecture is not supported.
Additionally, all limitations and tradeoffs of the edit distance routines should be provided upfront
so the user knows exactly what to expect.
Finally, this library should lead to performance boosts on both short and longer strings, so it
can be used for a variety of tasks, from bioinformatics to natural language processing.
`triple_accel` is very lightweight: it only has dependencies on other crates for benchmarking.
It can be built on machines without CPUs that have AVX2 or SSE4.1 support. It can also run on
machines without SIMD support by automatically using scalar alternatives.

## Install
Add
```
triple_accel = "*"
```
to the `[dependencies]` section of your `Cargo.toml`. This library is available
[here](https://crates.io/crates/triple_accel) on crates.io.

Alternatively, you can clone this repository and run
```
cargo build --release
```
In general, for maximum efficiency, use `RUSTFLAGS="-C target-cpu=native"` if portability is not an issue.

## Tests
You can run tests with
```
cargo test
```
after cloning the repository.

Continuous integration is used to ensure that the code passes all tests on the latest Linux, Windows,
and Mac platforms. Additionally, crate feature flags like `jewel-sse`, `jewel-avx`, `jewel-8bit`,
`jewel-16bit`, and `jewel-32bit` are used to override the default automatic detection of CPU features,
so all features can be thoroughly tested in continuous integration. The `debug` feature flag is specified,
so the exact underlying vector type that is used is printed.

## Benchmarks
Benchmarks can be ran with
```
cargo bench
```

## Docs
The docs are available [here](https://docs.rs/triple_accel). To build them on
your machine, run
```
cargo doc
```

## Features
This library provides routines for both searching for some needle string in a haystack string
and calculating the edit distance between two strings. Hamming distance (mismatches only),
Levenshtein distance (mismatches + gaps), and restricted Damerau-Levenshtein distance
(transpositions + mismatches + gaps) are supported, along with arbitrary edit costs. This
library provides a simple interface, in addition to powerful lower-level control over the edit
distance calculations.

At runtime, the implementation for a certain algorithm is selected based on CPU support, going
down the list:

1. Vectorized implementation with 256-bit AVX vectors, if AVX2 is supported.
2. Vectorized implementation with 128-bit SSE vectors, if SSE4.1 is supported.
3. Scalar implementation.

Currently, vectorized SIMD implementations are only available for x86 or x86-64 CPUs. However,
after compiling this library on a machine that supports those SIMD intrinsics, the library can
be used on other machines.
Additionally, the internal data structure for storing vectors and the bit width of the values
in the vectors are selected at runtime for maximum efficiency and accuracy, given the lengths
of the input strings.

## Limitations
Due to the use of SIMD intrinsics, only binary strings that are represented with `u8` bytes
are supported. Unicode strings are not currently supported.

## Examples
`triple_accel` provides a very simple and easy to use framework for common edit distance operations.
Calculating the Hamming distance (number of mismatches) between two strings is extremely simple:
```Rust
use triple_accel::*;

let a = b"abcd";
let b = b"abcc";

let dist = hamming(a, b);
assert!(dist == 1);
```
By default, SIMD will be used if possible.
Similarly, we can easily calculate the Levenshtein distance (character mismatches and gaps all have
a cost of 1) between two strings with the following code:
```Rust
use triple_accel::*;

let a = b"abc";
let b = b"abcd";

let dist = levenshtein_exp(a, b);
assert!(dist == 1);
```
This uses exponential search to estimate the number of edits between `a` and `b`, which makes it
more efficient than the alternative `levenshtein` function when the number of edits between `a`
and `b` is low.

In addition to edit distance routines, `triple_accel` also provides search routines. These routines
return an iterator over matches that indicate where the `needle` string matches the `haystack` string.
`triple_accel` will attempt to maximize the length of matches that end at the same position and remove
shorter matches when some matches fully overlap.
```Rust
use triple_accel::*;

let needle = b"helllo";
let haystack = b"hello world";

let matches: Vec<Match> = levenshtein_search(needle, haystack).collect();
// note: start index is inclusive, end index is exclusive!
assert!(matches == vec![Match{start: 0, end: 5, k: 1}]);
```
Sometimes, it is necessary to use the slightly lower level, but also more powerful routines that
`triple_accel` provides. For example, it is possible to allow transpositions (character swaps) that
have a cost of 1, in addition to mismatches and gaps:
```Rust
use triple_accel::levenshtein::*;

let a = b"abcd";
let b = b"abdc";
let k = 2; // upper bound on allowed cost
let trace_on = false; // return edit traceback?

let dist = levenshtein_simd_k_with_opts(a, b, k, trace_on, RDAMERAU_COSTS);
// note: dist may be None if a and b do not match within a cost of k
assert!(dist.unwrap().0 == 1);
```
Don't let the name of the function fool you! `levenshtein_simd_k_with_opts` will still fall back to
the scalar implementation if AVX2 or SSE4.1 support is not available. It just prefers to use SIMD
where possible.

For most common cases, the re-exported functions are enough, and the low level functions do not
have to be used directly.

## License
[MIT](LICENSE)

## Contributing
Read the contributing guidelines [here](CONTRIBUTING.md).

## Code of Conduct
Read the code of conduct [here](CODE_OF_CONDUCT.md).

## Why the name "triple_accel"?
Because "Time Altar - Triple Accel" is a magical ability used by Kiritsugu Emiya to boost his speed
and reaction time in Fate/Zero. There are also some other references to the Fate series...
