# legendre
Rust library for computing Legendre polynomials and associated Legendre polynomials

## Legendre Polynomials
[Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials) can be generated using arbitrary-precision (via `num_bigint` and `num_rational`) or `f64`. 
Arbitrary precision calculations are slower but `f64` calculations begin to accumulate appreciable error as `l` increases. For example, at `l=20` the error is approximately `1e-11` and at `l=46` the error is approximately `1.0`. 

## Associated Legendre Polynomials
[Associated Legendre polynomials](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials) can be generated individually using `AssociatedPolynomial<BigRational>::new(l,m,phase)`. If all associated Legendre polynomials are needed up to `l` and `m`, use `batch_associated_legendre_polynomials(l_max, x, phase)`, which uses recursion to increase the speed of computation. 

## Normalized Associated Legendre Polynomials
If unnormalized, arbitrary-precision associated Legendre polynomials are too computationally-intensive, normalized associated Legendre polynomials can be used. Use the `normalized` feature to enable the `legendre::normalized` module. 

The normaliztion in this crate is based on [a paper by Holmes and Featherstone](https://link.springer.com/article/10.1007/s00190-002-0216-2). 

To generate normalized associated Legendre polynomials, first create a `ScalarCache`, which computes and stores constants required for normalization. Then, use `fill` or `compute` to evalute the normalized associated Legendre polynomials for a given value. 

Note that the normalized associated Legendre polynomial generally uses an angle as an input, which is different than the unnormalized associated Legendre polynomial. 