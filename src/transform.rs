use std::ops::{Add, Mul};

/// Provides the multiplicative scale from a unit into its reference unit.
///
/// If `x` is expressed in `Self`, the corresponding value in the reference
/// unit is
///
/// ```text
/// x_reference = x * scale.
/// ```
///
/// For [`LengthUnit`], the reference unit is the centimetre.
///
/// [`LengthUnit`]: crate::LengthUnit
pub(crate) trait UnitScale<T> {
    /// Returns the multiplicative scale into the reference unit.
    fn scale(self) -> T;
}

/// The value and first two derivatives of a scalar transformation.
///
/// For a backend coordinate `q = q(x)`:
///
/// - `value = q(x)`
/// - `first = dq/dx`
/// - `second = d²q/dx²`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransformJet<T> {
    pub value: T,
    pub first: T,
    pub second: T,
}

impl<T> TransformJet<T> {
    #[must_use]
    pub const fn new(value: T, first: T, second: T) -> Self {
        Self {
            value,
            first,
            second,
        }
    }

    #[must_use]
    pub fn map<U>(self, mut f: impl FnMut(T) -> U) -> TransformJet<U> {
        TransformJet {
            value: f(self.value),
            first: f(self.first),
            second: f(self.second),
        }
    }

    pub fn square(self) -> Self
    where
        T: Clone + Add<Output = T> + Mul<Output = T>,
    {
        let two = |x: T| x.clone() + x;

        Self {
            value: self.value.clone() * self.value.clone(),
            first: two(self.value.clone() * self.first.clone()),
            second: two(self.first.clone() * self.first + self.value * self.second),
        }
    }

    /// Composes `outer(inner(x))`.
    #[must_use]
    pub fn compose(self, inner: Self) -> Self
    where
        T: Clone + Add<Output = T> + Mul<Output = T>,
    {
        let first = self.first.clone() * inner.first.clone();

        let second = self.second * inner.first.clone() * inner.first + self.first * inner.second;

        Self {
            value: self.value,
            first,
            second,
        }
    }
}

/// Transforms a first derivative from `q` coordinates into `x` coordinates.
///
/// Given `q = q(x)` and `df/dq`, returns:
///
/// `df/dx = df/dq * dq/dx`.
#[must_use]
pub fn transform_first<F, S>(df_dq: F, transform: TransformJet<S>) -> F
where
    F: Mul<S, Output = F>,
{
    df_dq * transform.first
}

/// Transforms first and second derivatives from `q` coordinates into `x`
/// coordinates.
///
/// Given `q = q(x)`, this returns:
///
/// `df/dx = df/dq * dq/dx`
///
/// `d²f/dx² = d²f/dq² * (dq/dx)² + df/dq * d²q/dx²`
#[must_use]
pub fn transform_second<F, S>(df_dq: F, d2f_dq2: F, transform: TransformJet<S>) -> (F, F)
where
    F: Clone + Add<Output = F> + Mul<S, Output = F>,
    S: Clone + Mul<Output = S>,
{
    let first = df_dq.clone() * transform.first.clone();

    let first_squared = transform.first.clone() * transform.first;

    let second = d2f_dq2 * first_squared + df_dq * transform.second;

    (first, second)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-14_f64.max(expected.abs() * 1.0e-13);

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}",
        );
    }

    #[test]
    fn transform_jet_constructor_sets_fields() {
        let jet = TransformJet::new(1.0, 2.0, 3.0);

        assert_eq!(jet.value, 1.0);
        assert_eq!(jet.first, 2.0);
        assert_eq!(jet.second, 3.0);
    }

    #[test]
    fn first_derivative_uses_chain_rule() {
        // f(q) = q³ at q = 2
        //
        // df/dq = 3q² = 12
        //
        // q(x) = x² at x = 3
        // dq/dx = 2x = 6
        //
        // df/dx = 12 × 6 = 72
        let coordinate = TransformJet::new(2.0_f64, 6.0, 2.0);

        let result = transform_first(12.0, coordinate);

        assert_close(result, 72.0);
    }

    #[test]
    fn second_derivative_uses_full_chain_rule() {
        // f(q) = q³, evaluated at q = 9:
        //
        // df/dq   = 3q² = 243
        // d²f/dq² = 6q  = 54
        //
        // q(x) = x², evaluated at x = 3:
        //
        // q       = 9
        // dq/dx   = 6
        // d²q/dx² = 2
        //
        // d²f/dx²
        // = 54 × 6² + 243 × 2
        // = 2430
        let coordinate = TransformJet::new(9.0_f64, 6.0, 2.0);

        let result = transform_second(243.0, 54.0, coordinate).1;

        assert_close(result, 2430.0);
    }

    #[test]
    fn linear_coordinate_transform_has_no_curvature_term() {
        let coordinate = TransformJet::new(10.0_f64, 4.0, 0.0);

        let result = transform_second(7.0, 3.0, coordinate).1;

        assert_close(result, 3.0 * 4.0_f64.powi(2));
    }

    #[test]
    fn linear_canonical_function_retains_coordinate_curvature() {
        let coordinate = TransformJet::new(9.0_f64, 6.0, 2.0);

        let result = transform_second(5.0, 0.0, coordinate).1;

        assert_close(result, 10.0);
    }

    #[test]
    fn constant_coordinate_has_zero_transformed_derivatives() {
        let coordinate = TransformJet::new(4.0_f64, 0.0, 0.0);

        assert_eq!(transform_first(12.0, coordinate), 0.0,);

        assert_eq!(transform_second(12.0, 7.0, coordinate,).1, 0.0,);
    }

    #[test]
    fn inverse_jets_compose_to_identity_to_second_order() {
        // q(x) = 1/x at x = 2.
        let forward = TransformJet::new(0.5_f64, -0.25, 0.25);

        // x(q) = 1/q at q = 0.5.
        let inverse = TransformJet::new(2.0_f64, -4.0, 16.0);

        let composed_first = inverse.first * forward.first;

        let composed_second =
            inverse.second * forward.first.powi(2) + inverse.first * forward.second;

        assert_close(composed_first, 1.0);
        assert_close(composed_second, 0.0);
    }

    #[test]
    fn derivative_helpers_support_f32() {
        let coordinate = TransformJet::new(2.0_f32, 3.0, 4.0);

        assert_eq!(transform_first(5.0_f32, coordinate), 15.0,);

        assert_eq!(
            transform_second(5.0_f32, 6.0, coordinate,).1,
            6.0 * 9.0 + 5.0 * 4.0,
        );
    }
}
