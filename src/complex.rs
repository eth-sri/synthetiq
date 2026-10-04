//! Small, dependency-free complex arithmetic. Deliberately avoids fused operations
//! so ordinary multiplication follows the reference implementation's rounding.
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use std::str::FromStr;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}
impl Complex {
    pub const ZERO: Self = Self::new(0.0, 0.0);
    pub const ONE: Self = Self::new(1.0, 0.0);
    pub const I: Self = Self::new(0.0, 1.0);
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    #[inline]
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    #[inline]
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    #[inline]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
}
impl From<f64> for Complex {
    fn from(x: f64) -> Self {
        Self::new(x, 0.0)
    }
}
impl Add for Complex {
    type Output = Self;
    #[inline]
    fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
}
impl Sub for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, b: Self) -> Self {
        Self::new(self.re - b.re, self.im - b.im)
    }
}
impl Neg for Complex {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
impl Mul for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
}
impl Mul<f64> for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, b: f64) -> Self {
        Self::new(self.re * b, self.im * b)
    }
}
impl Mul<Complex> for f64 {
    type Output = Complex;
    #[inline]
    fn mul(self, b: Complex) -> Complex {
        b * self
    }
}
impl Div<f64> for Complex {
    type Output = Self;
    #[inline]
    fn div(self, b: f64) -> Self {
        Self::new(self.re / b, self.im / b)
    }
}
impl Div for Complex {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        let d = b.norm_sqr();
        (self * b.conj()) / d
    }
}
impl AddAssign for Complex {
    #[inline]
    fn add_assign(&mut self, b: Self) {
        self.re += b.re;
        self.im += b.im;
    }
}
impl SubAssign for Complex {
    #[inline]
    fn sub_assign(&mut self, b: Self) {
        self.re -= b.re;
        self.im -= b.im;
    }
}
impl MulAssign for Complex {
    #[inline]
    fn mul_assign(&mut self, b: Self) {
        *self = *self * b;
    }
}
impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({},{})", self.re, self.im)
    }
}
impl FromStr for Complex {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if let Some(inner) = s.strip_prefix('(').and_then(|x| x.strip_suffix(')')) {
            if let Some((r, i)) = inner.split_once(',') {
                return Ok(Self::new(
                    r.trim()
                        .parse()
                        .map_err(|_| format!("Invalid complex number: {s}"))?,
                    i.trim()
                        .parse()
                        .map_err(|_| format!("Invalid complex number: {s}"))?,
                ));
            }
            return inner
                .parse::<f64>()
                .map(Self::from)
                .map_err(|_| format!("Invalid complex number: {s}"));
        }
        s.parse::<f64>()
            .map(Self::from)
            .map_err(|_| format!("Invalid complex number: {s}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn arithmetic() {
        let a = Complex::new(2.0, 3.0);
        let b = Complex::new(-4.0, 5.0);
        assert_eq!(a + b, Complex::new(-2.0, 8.0));
        assert_eq!(a * b, Complex::new(-23.0, -2.0));
        assert!(((a / b) * b - a).abs() < 1e-14);
        assert_eq!(a * a.conj(), Complex::new(13.0, 0.0));
        assert_eq!(a.norm_sqr(), 13.0);
    }
    #[test]
    fn parse_reference_format() {
        for (s, want) in [
            ("(1.25,-2e-3)", Complex::new(1.25, -0.002)),
            ("-0.5", Complex::new(-0.5, 0.0)),
            ("(2)", Complex::new(2.0, 0.0)),
        ] {
            assert_eq!(s.parse::<Complex>().unwrap(), want);
        }
        assert!("(nan,invalid)".parse::<Complex>().is_err());
    }
}
