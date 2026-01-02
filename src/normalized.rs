use triangle_container::{TriangleContainer, TriangleVec, iterators::DiagonalIterator};

/// Holmes and Featherstone Eq 12
fn a_coefficent(l: usize, m: usize) -> f64
{
    let numer = (2 * l - 1) * (2 * l + 1);
    let denom = (l - m) * (l + m);
    (numer as f64 / denom as f64).sqrt()
}

/// Holmes and Featherstone Eq 12
fn b_coefficient(l: usize, m: usize) -> f64
{
    let numer = (2 * l + 1) * (l + m - 1) * (l - m - 1);
    let denom = (l - m) * (l + m) * (2 * l - 3);
    (numer as f64 / denom as f64).sqrt()
}

/// Holmes and Featherstone Eq 13
fn sectoral_coefficent(l: usize) -> f64
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

    /// Evaluates the normalized associated Legendre polynomials and stores the results in the provided container
    pub fn fill<C>(&self, container: &mut C, angle: f64)
    where C: TriangleContainer<f64>
    {
        let u = angle.cos();
        let t = angle.sin();

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
    pub fn compute(&self, angle: f64) -> TriangleVec<f64>
    {
        let mut container = TriangleVec::new(self.0.max_l());
        self.fill(&mut container, angle);
        container
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use triangle_container::iterators::RowIterator;
    use factorial::Factorial;
    use num_rational::BigRational;
    use num_bigint::{BigInt, BigUint};
    use num_traits::ToPrimitive;
    use crate::Phase;

    const fn kronecker(m: usize) -> usize
    {
        match m {
            0 => 1,
            _ => 0
        }
    }

    fn normalization_coefficent(l: usize, m: usize) -> f64
    {
        let frac1 = BigRational::new(
            BigInt::from(1), 
            BigInt::from((2 - kronecker(m)) * (2 * l + 1))
        );
        let numer = Factorial::factorial(&BigUint::from(l + m));
        let denom = Factorial::factorial(&BigUint::from(l - m));
        let frac2 = BigRational::new(numer.into(), denom.into());
        let frac = frac1 * frac2;

        // Return the result
        frac.to_f64().unwrap().sqrt()
    }

    #[test]
    fn test_normalization() {
        let lat = 0.5_f64;
        let scalar_cache = ScalarCache::new(10);
        let alp = scalar_cache.compute(lat);
        let expected = crate::batch_associated_legendre_polynomials(10, lat.sin(), &Phase::None).unwrap();
        for (l, m) in RowIterator::new(10) {
            let lhs = alp.get(l, m) * normalization_coefficent(l, m);
            let rhs = expected.get(l, m);
            let error = (lhs - rhs).abs() / rhs.abs();
            assert!(error < 1e-10);
        }
    }
}