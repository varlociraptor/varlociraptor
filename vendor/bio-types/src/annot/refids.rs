//! Intern reference sequence (e.g., chromosome) names
use std::collections::HashMap;
use std::ops::Deref;

/// Data structure for interning sequence names efficiently.
///
/// The structure is parameterized over the reference type `R` used to
/// intern strings. Typically, this would be `Rc` for single-threaded
/// access or `Arc` for multi-threaded access. These reference types
/// provide fast, reference-counted cloning with no new allocation,
/// which can make sequence location calculations faster as well as
/// reducing the memory footprint required.
///
/// ```
/// use std::rc::Rc;
/// use bio_types::strand::ReqStrand;
/// use bio_types::annot::contig::Contig;
/// use bio_types::annot::loc::Loc;
/// use bio_types::annot::refids::RefIDSet;
/// let mut refids: RefIDSet<Rc<String>> = RefIDSet::new();
/// let pau8 = Contig::new(refids.intern("chrI"), 1807, 2170 - 1807, ReqStrand::Reverse);
/// {
///   let chr_i = refids.intern("chrI");
///   // One reference for the RefIDSet itself, one for the pau8 Contig, one for chr_i
///   assert_eq!(Rc::strong_count(&chr_i), 3);
/// }
/// let seo1 = Contig::new(refids.intern("chrI"), 7235, 9017 - 7235, ReqStrand::Reverse);
/// let tda8 = Contig::new(refids.intern("chrI"), 13363, 13744 - 13363, ReqStrand::Reverse);
/// {
///   let chr_i = refids.intern("chrI");
///   assert_eq!(Rc::strong_count(&chr_i), 5);
/// }
/// let seo1_beginning = seo1.first_pos();
/// let seo1_ending = seo1.last_pos();
/// {
///   let chr_i = refids.intern("chrI");
///   assert_eq!(Rc::strong_count(&chr_i), 7);
/// }
/// ```
pub struct RefIDSet<R> {
    refids: HashMap<String, R>,
}

impl<R> Default for RefIDSet<R> {
    fn default() -> Self {
        Self::new()
    }
}

impl<R> RefIDSet<R> {
    /// Create a new, empty table of interned reference names
    pub fn new() -> Self {
        RefIDSet {
            refids: HashMap::new(),
        }
    }

    /// Intern a reference name.
    ///
    /// This returns a shared reference of type `R` for the name. This
    /// reference will be shared with any other intern calls for the
    /// same name. The name is given originally as a reference, and it
    /// will be cloned into an owned `String` only when the name is
    /// new for the data type.
    pub fn intern(&mut self, id: &str) -> R
    where
        R: Deref<Target = String> + From<String> + Clone,
    {
        if self.refids.contains_key(id) {
            if let Some(ref r) = self.refids.get(id) {
                (*r).clone()
            } else {
                panic!("RefIDSet::ensure failed to get() after contains()");
            }
        } else {
            let r = R::from(id.to_owned());
            self.refids.insert(id.to_owned(), r.clone());
            r
        }
    }
}
