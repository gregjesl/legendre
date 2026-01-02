use triangle_container::{TriangleContainer, TriangleVec, iterators::DiagonalIterator};

/// Holmes and Featherstone Eq 12
pub fn a_coefficent(l: usize, m: usize) -> f64
{
    let numer = (2 * l - 1) * (2 * l + 1);
    let denom = (l - m) * (l + m);
    (numer as f64 / denom as f64).sqrt()
}

/// Holmes and Featherstone Eq 12
pub fn b_coefficient(l: usize, m: usize) -> f64
{
    let numer = (2 * l + 1) * (l + m - 1) * (l - m - 1);
    let denom = (l - m) * (l + m) * (2 * l - 3);
    (numer as f64 / denom as f64).sqrt()
}

/// Holmes and Featherstone Eq 12
pub fn coefficents(l: usize, m: usize) -> (f64, f64)
{
    let a = a_coefficent(l, m);
    let b = b_coefficient(l, m);
    (a, b)
}

/// Holmes and Featherstone Eq 13
pub fn sectoral_coefficent(l: usize) -> f64
{
    let numer = 2 * l + 1;
    let denom = 2 * l;
    (numer as f64 / denom as f64).sqrt()
}

/// Cache of scalars used to compute the associated Legendre polynomials
#[derive(Clone, Debug)]
pub struct ScalarCache(TriangleVec<(f64, Option<f64>)>);

impl ScalarCache {
    pub fn new(max_l: usize) -> Self {
        // Initialize the container
        let mut container = TriangleVec::new(max_l);

        // Store the first item
        container.set(0, 0, (1.0, None));

        // Compute the diagonal elements
        for l in 1..=max_l {
            container.set(l, l, (sectoral_coefficent(l), None));
        }

        for (l, m) in DiagonalIterator::offset(max_l, 1) {
            let a = a_coefficent(l, m);
            if m == l - 1 {
                container.set(l, m, (a, None));
            } else {
                let b = b_coefficient(l, m);
                container.set(l, m, (a, Some(b)));
            }
        }

        Self(container)
    }

    /// Given a triangle container and an input `x`, fill the container with the values of associated normalized Legendre polynomials
    pub fn fill<C>(&self, container: &mut C, x: f64)
    where C: TriangleContainer<f64>
    {
        assert_eq!(container.max_l(), self.0.max_l(), "Size mismatch");
        let t = x.sin();
        let u = x.cos();

        // Compute the diagonal first
        container.set(0, 0, 1.0);
        let mut prev = u * 3_f64.sqrt();
        container.set(1, 1, prev);
        for l in 2..=container.max_l() {
            let (scalar, _) = self.0.get(l, l);
            let value = u * scalar * prev;
            container.set(l, l, value);
            prev = value;
        }

        // Compute the off-diagonal
        for l in 1..=container.max_l() {
            let m = l - 1;
            let (a, _) = self.0.get(l, m);
            let value = t * a * container.get(l - 1, m);
            container.set(l, m, value);
        }

        for (l, m) in DiagonalIterator::offset(container.max_l(), 2) {
            let (a, Some(b)) = self.0.get(l, m) else { unreachable!() };
            let value = t * a * container.get(l - 1, m) - b * container.get(l - 2, m);
            container.set(l, m, value);
        }
    }

    /// Computes the normalized associated Legendre polynomials
    pub fn compute(&self, u: f64) -> TriangleVec<f64>
    {
        let mut container = TriangleVec::new(self.0.max_l());
        self.fill(&mut container, u);
        container
    }
}