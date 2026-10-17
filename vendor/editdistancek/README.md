# editdistancek

This library calculates the Levenshtein edit distance between two strings. The edit distance represents the minimum
number of edits required to transform one string into the other.

Our method is inspired by the algorithms of Ukkonen and Landau, Myers, and Schmidt. These strategies offer efficient
ways of calculating the edit distance.

The algorithm operates with a running time of O(d*min(n,m)) and O(d^2 + n) in expectation. Here, 'd' stands for the edit distance between the two
strings, while 'n' and 'm' denote the lengths of these strings, respectively. This efficient running time ensures the
algorithm remains performant even for larger strings.

## Motivation

The primary objective of this library is to devise a quick solution for the bounded version of the edit distance
problem. Specifically, it aims to handle scenarios where, given a threshold 'k', it returns the edit distance 'd' if 'd'
is less than or equal to 'k'. Otherwise, it indicates that the distance exceeds 'k'. This approach addresses the typical
shortcoming of the classic edit distance algorithm, which tends to be inefficient in these cases. 

The 'bounded' aspect of this library makes it particularly useful in scenarios where there's a predefined threshold of acceptable dissimilarity. This might be beneficial in time-sensitive applications or when working with very large strings, where calculating the full Levenshtein distance would be computationally prohibitive. 

P.S. The 'k' in editdistancek refers to the upper bound for the algorithm. In practical terms, this means that the algorithm only computes the Levenshtein edit distance up to this 'k' value.

## Installation

Add to Cargo.toml

```toml
[dependecies]
editdistancek = "1.*"
```

## Usage

```rust
edit_distance("kitten".as_bytes(), "sitting".as_bytes()); // => 3
edit_distance_bounded("kitten".as_bytes(), "sitting".as_bytes(), 3); // => Some(3)
edit_distance_bounded("kitten".as_bytes(), "sitting".as_bytes(), 2); // => None
```

## License

MIT License
