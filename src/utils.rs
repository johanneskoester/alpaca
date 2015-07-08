
pub fn matrix<T: Copy>(v: T, n: usize, m: usize) -> Vec<Vec<T>> {
    let mut matrix = vec![];
    for _ in 0..n {
        matrix.push(vec![v; m]);
    }

    matrix
}
