// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

#[derive(Debug, Clone)]
pub(crate) struct TumorNormalPair<T> {
    pub(crate) tumor: T,
    pub(crate) normal: T,
}

impl<T> From<TumorNormalPair<T>> for Vec<T> {
    fn from(tnp: TumorNormalPair<T>) -> Self {
        vec![tnp.tumor, tnp.normal]
    }
}
