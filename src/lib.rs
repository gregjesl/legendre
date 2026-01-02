use std::{collections::HashMap, ops::Neg};
use factorial::Factorial;
use num_rational::BigRational;
use num_bigint::{BigInt, BigUint};
use num_traits::ToPrimitive;
use compensated_summation::KahanBabuskaNeumaier;
use triangle_container::{TriangleVec, TriangleContainer, iterators::DiagonalIterator};

#[cfg(feature = "normalized")]
pub mod normalized;

/// Constructs a BigRational for {integer} values
macro_rules! ratio {
    ($numer:expr, $denom:expr) => {
        BigRational::new(BigInt::from($numer), BigInt::from($denom))
    };
}

/// Computes the leading coefficent of the Legendre polynomal
/// 
/// See Lemma 7.2 of [A Second Course in Ordinary Differential Equations](https://math.libretexts.org/Bookshelves/Differential_Equations/A_Second_Course_in_Ordinary_Differential_Equations%3A_Dynamical_Systems_and_Boundary_Value_Problems_(Herman)/07%3A_Special_Functions/7.02%3A_Legendre_Polynomials)
fn leading(n: u16) -> BigRational
{
    let nfact = BigInt::from(Factorial::factorial(&BigUint::from(n)));
    let twonfact = BigInt::from(Factorial::factorial(&BigUint::from(2 * (n as u32))));
    let twopown = BigInt::from(2).pow(n as u32);
    BigRational::new(1.into(), twopown * nfact.clone()) * BigRational::new(twonfact, nfact)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    CondonShortley,
    None,
}

/// Structure representing a polynomial
#[derive(Clone, Debug)]
pub struct Polynomial<V>
{
    coefficents: HashMap<u16, V>,
}

impl<C> Polynomial<C>
where C: Default + Clone
{
    /// Returns the order of the polynomial
    pub fn order(&self) -> u16
    {
        *self.coefficents.keys().max().unwrap_or(&0)
    }
}

impl<C> Polynomial<C>
where C: Neg<Output = C> + Clone
{
    pub fn negate(&self) -> Self
    {
        Self {
            coefficents: self.coefficents.iter()
                .map(|(k, v)| (*k, v.clone().neg()))
                .collect()
        }
    }
}

impl<C> Neg for Polynomial<C>
where C: Neg<Output = C> + Clone
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.negate()
    }
}

impl Polynomial<BigRational>
{
    /// Generates the Legendre polynomial of order `n`
    pub fn new(n: u16) -> Self
    {
        let mut coefficents = HashMap::new();

        // Compute the leading coefficent
        coefficents.insert(n, leading(n));

        let mut k = n;
        while k > 1 {
            let mut value = coefficents.get(&k).unwrap().clone();
            let numer = BigInt::from((k as i64) * (k - 1) as i64);
            let denom = BigInt::from((n as i64 - k as i64 + 2_i64) * (n as i64 + k as i64 - 1));
            value *= -BigRational::new(numer, denom);
            k -= 2;
            coefficents.insert(k, value);
        }
        Self { coefficents }
    }

    /// Returns the coefficent of x^n
    pub fn coefficent(&self, n: u16) -> BigRational
    {
        self.coefficents.get(&n).unwrap_or(&BigRational::default()).clone()
    }

    /// Differentiates the polynomial
    pub fn differentiate(&self) -> Self
    {
        let mut result = HashMap::new();
        for (power, value) in &self.coefficents {
            if *power == 0 {
                continue;
            }
            result.insert(power - 1, value.clone() * ratio!(*power, 1));
        }
        return Self { coefficents: result }
    }

    /// Evaluates the polynomial
    /// 
    /// # Returns
    /// - `None` if `x > 1.0` or `x < -1.0`
    /// - `None` if `x` cannot be converted to a rational fraction (see [`Polynomial::approx_evaluate`])
    /// - `Some(result)` upon success
    pub fn evaluate(&self, x: f64) -> Option<BigRational>
    {
        if x > 1.0 || x < -1.0 { return None }
        let xfrac = BigRational::from_float(x)?;
        let mut sum = ratio!(0, 1);
        for (power, coefficent) in &self.coefficents {
            sum += coefficent * xfrac.pow(*power as i32);
        }
        Some(sum)
    }

    /// Attemps to evaluate the polynomial, moving to the next possible value if the evaluation fails
    /// 
    /// This method will panic if `x > 1.0` or `x < -1.0`
    pub fn approx_evaluate(&self, x: f64) -> BigRational
    {
        assert!(x <= 1.0 && x >= -1.0);
        match self.evaluate(x) {
            Some(value) => value,
            None => self.approx_evaluate(x.next_up()),
        }
    }
}

impl<V> IntoIterator for Polynomial<V>
{
    type Item = (u16, V);
    type IntoIter = std::collections::hash_map::IntoIter<u16, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.coefficents.into_iter()
    }
}

impl Polynomial<f64>
{
    /// Evaluates the polynomial
    /// 
    /// # Returns
    /// - `None` if `x > 1.0` or `x < -1.0`
    /// - `None` if `x` cannot be converted to a rational fraction (see [`Polynomial::approx_evaluate`])
    /// - `Some(result)` upon success
    pub fn evaluate(&self, x: f64) -> Option<f64>
    {
        if x > 1.0 || x < -1.0 { return None }
        let mut sum = KahanBabuskaNeumaier::new();
        for (power, coefficent) in &self.coefficents {
            sum += coefficent * x.powi(*power as i32);
        }
        Some(sum.total())
    }

    /// Attemps to evaluate the polynomial, moving to the next possible value if the evaluation fails
    /// 
    /// This method will panic if `x > 1.0` or `x < -1.0`
    pub fn approx_evaluate(&self, x: f64) -> f64
    {
        assert!(x <= 1.0 && x >= -1.0);
        match self.evaluate(x) {
            Some(value) => value,
            None => self.approx_evaluate(x.next_up()),
        }
    }
}

impl TryFrom<Polynomial<BigRational>> for Polynomial<f64>
{
    type Error = ();
    fn try_from(value: Polynomial<BigRational>) -> Result<Self, Self::Error> {
        let mut result = HashMap::new();
        for (n, coefficent) in value.coefficents {
            let coeff = coefficent.to_f64().ok_or(())?;
            result.insert(n, coeff);
        };
        Ok(Self { coefficents: result })
    }
}

#[derive(Clone, Debug)]
pub struct AssociatedPolynomial<C>
{
    poly: Polynomial<C>,
    m: u16
}

impl AssociatedPolynomial<BigRational>
{
    /// Generates the associated Legendre polynomial of degree `l` and order `m`
    /// 
    /// - Condon–Shortley phase is typically used in quantom mechanics (set `true`)
    /// - Condon–Shortley phase is not typically used in geophysics (set `false`)
    pub fn new(l: u16, m: u16, phase: &Phase) -> Result<Self, ()>
    {
        if m > l { return Err(()); }
        let mut poly = Polynomial::new(l);
        for _ in 0..m {
            poly = poly.differentiate();
        }
        if phase == &Phase::CondonShortley && m % 2 == 1 {
            poly = -poly;
        }
        Ok(AssociatedPolynomial { poly, m })
    }

    /// Evaluates the polynomial
    /// 
    /// # Returns
    /// - `None` if `x > 1.0` or `x < -1.0`
    /// - `None` if `x` cannot be converted to a rational fraction (see [`AssociatedPolynomial::approx_evaluate`])
    /// - `Some(result)` upon success
    pub fn evaluate(&self, x: f64) -> Option<BigRational>
    {
        let term1 = BigRational::from_float((1.0 - x.powi(2)).powf(self.m as f64 / 2.0))?;
        let term2 = self.poly.evaluate(x)?;
        return Some(term1 * term2);
    }

    /// Attemps to evaluate the polynomial, moving to the next possible value if the evaluation fails
    /// 
    /// This method will panic if `x > 1.0` or `x < -1.0`
    pub fn approx_evaluate(&self, x: f64) -> BigRational
    {
        assert!(x <= 1.0 && x >= -1.0);
        match self.evaluate(x) {
            Some(value) => value,
            None => self.approx_evaluate(x.next_up()),
        }
    }
}

pub fn batch_associated_legendre_polynomials(lmax: u16, x: f64, phase: &Phase) -> Option<TriangleVec<f64>>
{
    if x > 1.0 || x < -1.0 { return None; }
    let mut rational = TriangleVec::new(lmax as usize + 1);
    rational.set(0, 0, ratio!(1, 1));
    let x_ratio = BigRational::from_float(x)?;

    // Compute the diagonal
    {
        let scalar = BigRational::from_float((1.0 - x * x).sqrt())?;
        let mut last = ratio!(1, 1);
        for l in 0..lmax as usize {
            let mut value = ratio!(2 * l as i32 + 1, 1) * &scalar * &last;
            match phase {
                &Phase::CondonShortley => value = -value,
                &Phase::None => { }
            }
            rational.set(l + 1, l + 1, value.clone());
            last = value;
        }
    }
    
    // Compute the off-diagonal
    {
        for l in 0..=lmax as usize {
            let value = &x_ratio * ratio!(2 * l as i32 + 1, 1) * &rational.get(l, l);
            debug_assert!(&rational.get(l, l) != &ratio!(0, 1));
            rational.set(l + 1, l, value);
        }
    }

    // Compute the rest
    {
        for (l, m) in DiagonalIterator::offset((lmax - 1) as usize, 1) {
            let term1 = ratio!(2 * l as i32 + 1, 1) * &x_ratio * &rational.get(l, m);
            let term2 = ratio!(l as i32 + m as i32, 1) * &rational.get(l - 1, m);
            let denom = ratio!(l as i32 - m as i32 + 1, 1);
            let value = (term1 - term2) / denom;
            rational.set(l + 1, m, value);
        }
    }

    let mut result = TriangleVec::new(lmax as usize);
    for (l, m) in DiagonalIterator::new(lmax as usize) {
        result.set(l, m, rational.get(l, m).to_f64()?);
    }

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::ToPrimitive;

    #[test]
    fn test_p0() {
        let poly = Polynomial::new(0);
        assert_eq!(poly.order(), 0);
        assert_eq!(poly.coefficent(0), ratio!(1, 1));

        assert_eq!(poly.evaluate(1.0).unwrap().to_f64().unwrap(), 1.0);
        assert_eq!(poly.evaluate(0.5).unwrap().to_f64().unwrap(), 1.0);
        assert_eq!(poly.evaluate(0.0).unwrap().to_f64().unwrap(), 1.0);
        assert_eq!(poly.evaluate(-0.5).unwrap().to_f64().unwrap(), 1.0);
        assert_eq!(poly.evaluate(-1.0).unwrap().to_f64().unwrap(), 1.0);
    }

    #[test]
    fn test_p1() {
        let poly = Polynomial::new(1);
        assert_eq!(poly.order(), 1);
        assert_eq!(poly.coefficent(1), ratio!(1, 1));
        assert_eq!(poly.coefficent(0), ratio!(0, 1));

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(0.5).unwrap(), ratio!(1, 2));
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(0, 1));
        assert_eq!(poly.evaluate(-0.5).unwrap(), ratio!(-1, 2));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(-1, 1));
    }

    #[test]
    fn test_p2() {
        let poly = Polynomial::new(2);
        assert_eq!(poly.order(), 2);
        assert_eq!(poly.coefficent(2), ratio!(3, 2));
        assert_eq!(poly.coefficent(1), ratio!(0, 1));
        assert_eq!(poly.coefficent(0), ratio!(-1, 2));

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(0.5).unwrap(), ratio!(-1, 8));
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(-1, 2));
        assert_eq!(poly.evaluate(-0.5).unwrap(), ratio!(-1, 8));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(1, 1));
    }

    #[test]
    fn test_p3() {
        let poly = Polynomial::new(3);
        assert_eq!(poly.order(), 3);
        assert_eq!(poly.coefficent(3), ratio!(5, 2));
        assert_eq!(poly.coefficent(2), ratio!(0, 1));
        assert_eq!(poly.coefficent(1), ratio!(-3, 2));
        assert_eq!(poly.coefficent(0), ratio!(0, 1));
    }

    #[test]
    fn test_p4() {
        let poly = Polynomial::new(4);
        assert_eq!(poly.order(), 4);
        assert_eq!(poly.coefficent(4), ratio!(35, 8));
        assert_eq!(poly.coefficent(3), ratio!(0, 1));
        assert_eq!(poly.coefficent(2), ratio!(-30, 8));
        assert_eq!(poly.coefficent(1), ratio!(0, 1));
        assert_eq!(poly.coefficent(0), ratio!(3, 8));
    }

    #[test]
    fn test_p00() {
        let poly = AssociatedPolynomial::new(0, 0, &Phase::CondonShortley).unwrap();

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(0.5).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(-0.5).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(1, 1));
    }

    #[test]
    fn test_p10() {
        let poly = AssociatedPolynomial::new(1, 0, &Phase::CondonShortley).unwrap();

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(0.5).unwrap(), ratio!(1, 2));
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(0, 1));
        assert_eq!(poly.evaluate(-0.5).unwrap(), ratio!(-1, 2));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(-1, 1));
    }

    #[test]
    fn test_p11() {
        let poly = AssociatedPolynomial::new(1, 1, &Phase::CondonShortley).unwrap();

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(0, 1));
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(-1, 1));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(0, 1));
    }

    #[test]
    fn test_p20() {
        let poly = AssociatedPolynomial::new(2, 0, &Phase::CondonShortley).unwrap();

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(1.0).unwrap().to_f64().unwrap(), 1.0);
        assert_eq!(poly.evaluate(0.5).unwrap(), ratio!(-1, 8));
        assert_eq!(poly.evaluate(0.5).unwrap().to_f64().unwrap(), -0.125);
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(-1, 2));
        assert_eq!(poly.evaluate(-0.5).unwrap(), ratio!(-1, 8));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(1, 1));

        let poly = AssociatedPolynomial::new(2, 0, &Phase::None).unwrap();
        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(1, 1));
        assert_eq!(poly.evaluate(1.0).unwrap().to_f64().unwrap(), 1.0);
        assert_eq!(poly.evaluate(0.5).unwrap(), ratio!(-1, 8));
        assert_eq!(poly.evaluate(0.5).unwrap().to_f64().unwrap(), -0.125);
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(-1, 2));
        assert_eq!(poly.evaluate(-0.5).unwrap(), ratio!(-1, 8));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(1, 1));
    }

    #[test]
    fn test_p21() {
        let poly = AssociatedPolynomial::new(2, 1, &Phase::CondonShortley).unwrap();

        assert_eq!(poly.evaluate(1.0).unwrap(), ratio!(0, 1));
        assert_eq!(poly.evaluate(0.0).unwrap(), ratio!(0, 1));
        assert_eq!(poly.evaluate(-1.0).unwrap(), ratio!(0, 1));
    }

    #[test]
    fn test_p22() {
        let poly = AssociatedPolynomial::new(2, 2, &Phase::CondonShortley).unwrap();

        assert_eq!(poly.evaluate(1.0).unwrap().to_f64().unwrap(), 0.0);
        assert_eq!(poly.evaluate(0.5).unwrap().to_f64().unwrap(), 9.0 / 4.0);
        assert_eq!(poly.evaluate(0.0).unwrap().to_f64().unwrap(), 3.0);
        assert_eq!(poly.evaluate(-0.5).unwrap().to_f64().unwrap(), 9.0 / 4.0);
        assert_eq!(poly.evaluate(-1.0).unwrap().to_f64().unwrap(), 0.0);
    }

    #[test]
    fn test_associated_error() {
        assert!(AssociatedPolynomial::new(2, 3, &Phase::CondonShortley).is_err());
    }

    #[test]
    fn test_batch_legendre() {
        let lmax = 4;
        let x = 0.5;
        let result = batch_associated_legendre_polynomials(lmax, x, &Phase::CondonShortley).unwrap();

        for (l, m) in DiagonalIterator::new(lmax as usize) {
            let p = AssociatedPolynomial::new(l as u16, m as u16, &Phase::CondonShortley).unwrap();
            let expected = p.approx_evaluate(x).to_f64().unwrap();
            let computed = result.get(l, m);
            assert!((expected - computed).abs() < 1e-10);
        }
    }
}