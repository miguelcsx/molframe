//! A rooted tree and its Newick serialisation.
//!
//! Clustering methods build a tree of relationships; Newick is the text format
//! that records one. The tree here is strictly binary — every internal node
//! joins exactly two children with their branch lengths — which is what the
//! distance-based methods in this crate produce.

/// A rooted binary tree with named leaves and branch lengths.
#[derive(Clone, Debug, PartialEq)]
pub enum Tree {
    /// A named tip.
    Leaf {
        /// The taxon name.
        name: String,
    },
    /// An internal node joining two subtrees, each by a branch length.
    Clade {
        /// The left subtree.
        left: Box<Tree>,
        /// Length of the branch to the left subtree.
        left_length: f64,
        /// The right subtree.
        right: Box<Tree>,
        /// Length of the branch to the right subtree.
        right_length: f64,
    },
}

/// Child-order convention for deterministic ladderisation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LadderDirection {
    /// Smaller clades precede larger clades.
    Ascending,
    /// Larger clades precede smaller clades.
    Descending,
}

/// Depth-first node visitation order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraversalOrder {
    /// Visit a clade before its descendants.
    Preorder,
    /// Visit descendants before their clade.
    Postorder,
}

/// Why a tree cannot be rerooted at the requested leaf edge.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum RerootError {
    /// No leaf has the requested exact name.
    LeafNotFound,
    /// More than one leaf has the requested name.
    AmbiguousLeaf,
    /// Root position is non-finite or outside `[0, 1]` from leaf to neighbour.
    InvalidFraction,
    /// An existing branch length is negative or non-finite.
    InvalidBranchLength,
    /// Input topology cannot form a rooted binary tree.
    InvalidTopology,
}

impl std::fmt::Display for RerootError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LeafNotFound => formatter.write_str("reroot leaf not found"),
            Self::AmbiguousLeaf => formatter.write_str("reroot leaf name is ambiguous"),
            Self::InvalidFraction => {
                formatter.write_str("reroot fraction must be finite and in [0, 1]")
            }
            Self::InvalidBranchLength => {
                formatter.write_str("tree has a negative or non-finite branch")
            }
            Self::InvalidTopology => {
                formatter.write_str("tree topology cannot be rerooted as binary")
            }
        }
    }
}

impl std::error::Error for RerootError {}

impl Tree {
    /// Renders the tree in Newick format, terminated by a semicolon.
    ///
    /// Branch lengths are written with the default floating-point formatting, so
    /// a whole number appears without a decimal point.
    #[must_use]
    pub fn to_newick(&self) -> String {
        let mut out = String::new();
        self.write(&mut out);
        out.push(';');
        out
    }

    /// Number of leaves below this node.
    #[must_use]
    pub fn leaf_count(&self) -> usize {
        match self {
            Self::Leaf { .. } => 1,
            Self::Clade { left, right, .. } => left.leaf_count() + right.leaf_count(),
        }
    }

    /// Visits every node in an explicit depth-first order.
    #[must_use]
    pub fn traverse(&self, order: TraversalOrder) -> Vec<&Tree> {
        let mut nodes = Vec::new();
        collect_nodes(self, order, &mut nodes);
        nodes
    }

    /// Sorts sibling clades by size, breaking equal-size ties by leaf name.
    pub fn ladderize(&mut self, direction: LadderDirection) {
        let Self::Clade { left, right, .. } = self else {
            return;
        };
        left.ladderize(direction);
        right.ladderize(direction);
        let left_key = (left.leaf_count(), minimum_leaf_name(left));
        let right_key = (right.leaf_count(), minimum_leaf_name(right));
        let swap = match direction {
            LadderDirection::Ascending => {
                left_key.0 > right_key.0 || (left_key.0 == right_key.0 && left_key.1 > right_key.1)
            }
            LadderDirection::Descending => {
                left_key.0 < right_key.0 || (left_key.0 == right_key.0 && left_key.1 > right_key.1)
            }
        };
        if swap {
            std::mem::swap(left, right);
        }
    }

    /// Places a new root along the pendant edge of an exactly named leaf.
    ///
    /// `fraction_from_leaf` is zero at the leaf and one at its current
    /// neighbour. Pairwise path lengths are preserved exactly.
    ///
    /// # Errors
    ///
    /// Returns [`RerootError`] for absent/duplicate names, invalid fractions,
    /// invalid branch lengths or a non-binary topology.
    pub fn reroot_at_leaf(&self, leaf: &str, fraction_from_leaf: f64) -> Result<Tree, RerootError> {
        if !fraction_from_leaf.is_finite() || !(0.0..=1.0).contains(&fraction_from_leaf) {
            return Err(RerootError::InvalidFraction);
        }
        let mut graph = Graph::from_tree(self)?;
        let matches = graph
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(index, node)| (node.name.as_deref() == Some(leaf)).then_some(index))
            .collect::<Vec<_>>();
        let target = match matches.as_slice() {
            [] => return Err(RerootError::LeafNotFound),
            [target] => *target,
            _ => return Err(RerootError::AmbiguousLeaf),
        };
        let [(neighbour, length)] = graph.nodes[target].edges.as_slice() else {
            return Err(RerootError::InvalidTopology);
        };
        let neighbour = *neighbour;
        let length = *length;
        graph.remove_edge(target, neighbour);
        let left = Tree::Leaf {
            name: leaf.to_owned(),
        };
        let right = graph.directed_tree(neighbour, target)?;
        Ok(Tree::Clade {
            left: Box::new(left),
            left_length: length * fraction_from_leaf,
            right: Box::new(right),
            right_length: length * (1.0 - fraction_from_leaf),
        })
    }

    fn write(&self, out: &mut String) {
        match self {
            Tree::Leaf { name } => out.push_str(name),
            Tree::Clade {
                left,
                left_length,
                right,
                right_length,
            } => {
                out.push('(');
                left.write(out);
                out.push(':');
                out.push_str(&left_length.to_string());
                out.push(',');
                right.write(out);
                out.push(':');
                out.push_str(&right_length.to_string());
                out.push(')');
            }
        }
    }
}

fn collect_nodes<'a>(tree: &'a Tree, order: TraversalOrder, nodes: &mut Vec<&'a Tree>) {
    if order == TraversalOrder::Preorder {
        nodes.push(tree);
    }
    if let Tree::Clade { left, right, .. } = tree {
        collect_nodes(left, order, nodes);
        collect_nodes(right, order, nodes);
    }
    if order == TraversalOrder::Postorder {
        nodes.push(tree);
    }
}

fn minimum_leaf_name(tree: &Tree) -> &str {
    match tree {
        Tree::Leaf { name } => name,
        Tree::Clade { left, right, .. } => minimum_leaf_name(left).min(minimum_leaf_name(right)),
    }
}

#[derive(Clone)]
struct GraphNode {
    name: Option<String>,
    edges: Vec<(usize, f64)>,
}

struct Graph {
    nodes: Vec<GraphNode>,
}

impl Graph {
    fn from_tree(tree: &Tree) -> Result<Self, RerootError> {
        let mut graph = Self { nodes: Vec::new() };
        let root = graph.add(tree)?;
        if graph.nodes[root].name.is_none() && graph.nodes[root].edges.len() == 2 {
            let first = graph.nodes[root].edges[0];
            let second = graph.nodes[root].edges[1];
            graph.remove_edge(root, first.0);
            graph.remove_edge(root, second.0);
            graph.connect(first.0, second.0, first.1 + second.1);
        }
        Ok(graph)
    }

    fn add(&mut self, tree: &Tree) -> Result<usize, RerootError> {
        let index = self.nodes.len();
        self.nodes.push(GraphNode {
            name: None,
            edges: Vec::new(),
        });
        match tree {
            Tree::Leaf { name } => self.nodes[index].name = Some(name.clone()),
            Tree::Clade {
                left,
                left_length,
                right,
                right_length,
            } => {
                for (child, length) in [
                    (left.as_ref(), *left_length),
                    (right.as_ref(), *right_length),
                ] {
                    if !length.is_finite() || length < 0.0 {
                        return Err(RerootError::InvalidBranchLength);
                    }
                    let child = self.add(child)?;
                    self.connect(index, child, length);
                }
            }
        }
        Ok(index)
    }

    fn connect(&mut self, first: usize, second: usize, length: f64) {
        self.nodes[first].edges.push((second, length));
        self.nodes[second].edges.push((first, length));
    }

    fn remove_edge(&mut self, first: usize, second: usize) {
        self.nodes[first].edges.retain(|(node, _)| *node != second);
        self.nodes[second].edges.retain(|(node, _)| *node != first);
    }

    fn directed_tree(&self, node: usize, parent: usize) -> Result<Tree, RerootError> {
        if let Some(name) = &self.nodes[node].name {
            return Ok(Tree::Leaf { name: name.clone() });
        }
        let children = self.nodes[node]
            .edges
            .iter()
            .filter(|(child, _)| *child != parent)
            .copied()
            .collect::<Vec<_>>();
        let [(left, left_length), (right, right_length)] = children.as_slice() else {
            return Err(RerootError::InvalidTopology);
        };
        Ok(Tree::Clade {
            left: Box::new(self.directed_tree(*left, node)?),
            left_length: *left_length,
            right: Box::new(self.directed_tree(*right, node)?),
            right_length: *right_length,
        })
    }
}

/// Why a Newick string could not be parsed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NewickError {
    /// A character was not what the grammar expected at that position.
    Unexpected(usize),
    /// The string ended in the middle of a tree.
    UnexpectedEnd,
    /// A branch length was not a number.
    BadLength(usize),
    /// A leaf had no name.
    EmptyName(usize),
}

impl Tree {
    /// Parses a Newick string into a strictly binary tree.
    ///
    /// Accepts the shape this module writes: named leaves joined by internal
    /// nodes that each carry two children with branch lengths, ending in an
    /// optional semicolon. Multifurcations and internal-node names are not
    /// accepted.
    ///
    /// # Errors
    ///
    /// Returns [`NewickError`] describing where and how the text departed from
    /// that grammar.
    pub fn from_newick(text: &str) -> Result<Tree, NewickError> {
        let mut parser = Parser {
            bytes: text.as_bytes(),
            position: 0,
        };
        let tree = parser.node()?;
        if parser.peek() == Some(b';') {
            parser.position += 1;
        }
        if parser.position != parser.bytes.len() {
            return Err(NewickError::Unexpected(parser.position));
        }
        Ok(tree)
    }
}

/// A cursor over the bytes of a Newick string.
struct Parser<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.position).copied()
    }

    fn consume(&mut self, byte: u8) -> Result<(), NewickError> {
        match self.peek() {
            Some(found) if found == byte => {
                self.position += 1;
                Ok(())
            }
            Some(_) => Err(NewickError::Unexpected(self.position)),
            None => Err(NewickError::UnexpectedEnd),
        }
    }

    fn node(&mut self) -> Result<Tree, NewickError> {
        if self.peek() == Some(b'(') {
            self.clade()
        } else {
            self.leaf()
        }
    }

    fn clade(&mut self) -> Result<Tree, NewickError> {
        self.consume(b'(')?;
        let left = self.node()?;
        self.consume(b':')?;
        let left_length = self.length()?;
        self.consume(b',')?;
        let right = self.node()?;
        self.consume(b':')?;
        let right_length = self.length()?;
        self.consume(b')')?;
        Ok(Tree::Clade {
            left: Box::new(left),
            left_length,
            right: Box::new(right),
            right_length,
        })
    }

    fn leaf(&mut self) -> Result<Tree, NewickError> {
        let start = self.position;
        while let Some(byte) = self.peek() {
            if matches!(byte, b':' | b',' | b')' | b'(' | b';') {
                break;
            }
            self.position += 1;
        }
        if self.position == start {
            return Err(NewickError::EmptyName(start));
        }
        let name = String::from_utf8_lossy(&self.bytes[start..self.position]).into_owned();
        Ok(Tree::Leaf { name })
    }

    fn length(&mut self) -> Result<f64, NewickError> {
        let start = self.position;
        while let Some(byte) = self.peek() {
            if matches!(byte, b'0'..=b'9' | b'.' | b'-' | b'+' | b'e' | b'E') {
                self.position += 1;
            } else {
                break;
            }
        }
        let text = String::from_utf8_lossy(&self.bytes[start..self.position]);
        match text.parse::<f64>() {
            Ok(value) => Ok(value),
            Err(_) => Err(NewickError::BadLength(start)),
        }
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
