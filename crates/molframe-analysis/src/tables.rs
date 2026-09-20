//! Shared implementation for compact structure-of-arrays result tables.

macro_rules! define_soa_table {
    ($(#[$meta:meta])* $visibility:vis struct $table:ident for $row:ident {
        $($(#[$field_meta:meta])* $field:ident: $field_type:ty),+ $(,)?
    }) => {
        $(#[$meta])*
        #[derive(Clone, Debug, Default, PartialEq)]
        $visibility struct $table {
            $($field: Vec<$field_type>),+
        }

        impl $table {
            /// Starts an empty table with capacity for `capacity` rows.
            #[must_use]
            pub fn with_capacity(capacity: usize) -> Self {
                Self { $($field: Vec::with_capacity(capacity)),+ }
            }

            /// Number of aligned rows.
            #[must_use]
            pub fn len(&self) -> usize {
                let mut lengths = [$(self.$field.len()),+].into_iter();
                let first = match lengths.next() {
                    Some(length) => length,
                    None => 0,
                };
                debug_assert!(lengths.all(|length| length == first));
                first
            }

            /// Whether the table has no rows.
            #[must_use]
            pub fn is_empty(&self) -> bool {
                self.len() == 0
            }

            $(
                $(#[$field_meta])*
                #[must_use]
                pub fn $field(&self) -> &[$field_type] {
                    &self.$field
                }
            )+

            /// Reconstructs one row without allocating.
            #[must_use]
            pub fn row(&self, index: usize) -> Option<$row> {
                Some($row { $($field: *self.$field.get(index)?),+ })
            }

            /// Iterates rows in deterministic table order.
            #[must_use]
            pub fn iter(&self) -> impl ExactSizeIterator<Item = $row> + '_ {
                (0..self.len()).map(|index| $row {
                    $($field: self.$field[index]),+
                })
            }

            /// Appends one row while preserving column alignment.
            pub fn push(&mut self, row: $row) {
                $(self.$field.push(row.$field));+
            }

            /// Moves a complete reducer block into this table.
            pub fn append(&mut self, other: &mut Self) {
                $(self.$field.append(&mut other.$field));+
            }
        }

        impl FromIterator<$row> for $table {
            fn from_iter<I: IntoIterator<Item = $row>>(rows: I) -> Self {
                let iterator = rows.into_iter();
                let (lower, _) = iterator.size_hint();
                let mut table = Self::with_capacity(lower);
                for row in iterator {
                    table.push(row);
                }
                table
            }
        }
    };
}
