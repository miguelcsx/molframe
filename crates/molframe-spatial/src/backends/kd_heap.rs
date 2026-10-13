//! Bounded max-heap of nearest candidates.

/// Fixed-capacity max-heap ordered by `(squared distance, atom)`.
#[derive(Debug)]
pub(super) struct KnnHeap {
    capacity: usize,
    limit: f32,
    unique: bool,
    items: Vec<(f32, u32)>,
}

#[inline]
fn less(a: (f32, u32), b: (f32, u32)) -> bool {
    a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).is_lt()
}

impl KnnHeap {
    pub(super) fn new(capacity: usize) -> Self {
        Self {
            capacity,
            limit: f32::INFINITY,
            unique: false,
            items: Vec::with_capacity(capacity.min(1 << 16) + 1),
        }
    }

    pub(super) fn reset(&mut self, limit: f32, unique: bool) {
        self.items.clear();
        self.limit = limit;
        self.unique = unique;
    }

    #[inline]
    pub(super) fn current_bound(&self) -> f32 {
        if self.items.len() >= self.capacity {
            self.items.first().map_or(self.limit, |item| item.0)
        } else {
            self.limit
        }
    }

    fn sift_up(&mut self, mut child: usize) {
        while child > 0 {
            let parent = (child - 1) / 2;
            if less(self.items[parent], self.items[child]) {
                self.items.swap(parent, child);
                child = parent;
            } else {
                break;
            }
        }
    }

    fn sift_down(&mut self, mut parent: usize) {
        let len = self.items.len();
        loop {
            let mut largest = parent;
            for child in [2 * parent + 1, 2 * parent + 2] {
                if child < len && less(self.items[largest], self.items[child]) {
                    largest = child;
                }
            }
            if largest == parent {
                break;
            }
            self.items.swap(parent, largest);
            parent = largest;
        }
    }

    #[inline]
    pub(super) fn push(&mut self, squared: f32, atom: u32) {
        if squared > self.limit {
            return;
        }
        if self.unique
            && let Some(position) = self.items.iter().position(|item| item.1 == atom)
        {
            if squared < self.items[position].0 {
                self.items[position].0 = squared;
                self.sift_down(position);
            }
            return;
        }
        let candidate = (squared, atom);
        if self.items.len() < self.capacity {
            self.items.push(candidate);
            let last = self.items.len() - 1;
            self.sift_up(last);
        } else if self.items.first().is_some_and(|&top| less(candidate, top)) {
            self.items[0] = candidate;
            self.sift_down(0);
        }
    }

    pub(super) fn items(&self) -> &[(f32, u32)] {
        &self.items
    }

    /// Items in increasing `(squared distance, atom)` order.
    pub(super) fn sorted(mut self) -> Vec<(f32, u32)> {
        self.items
            .sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        self.items
    }

    /// Sorts in place and returns the ordered items.
    pub(super) fn sort_in_place(&mut self) -> &[(f32, u32)] {
        self.items
            .sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        &self.items
    }
}
