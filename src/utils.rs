pub fn log_prob_sum(values: &[f64]) -> f64 {
    let p0 = values.max();
    p0 + (values.map(|p| (p - p0).exp()).sum() - 1).ln_1p()
}


pub fn matrix<T: Copy>(v: T, n: usize, m: usize) -> Vec<Vec<T>> {
    let mut matrix = vec![];
    for _ in 0..n {
        matrix.push(vec![f64::MIN; m]);
    }

    matrix
}
