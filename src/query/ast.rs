pub enum Expr {
    SampleUnion(Vec<String>),
    Diff(Box<Expr>, Box<Expr>),
    Union(Box<Expr>, Box<Expr>),
    RelaxedIntersection(Vec<Box<Expr>>, usize)
}
