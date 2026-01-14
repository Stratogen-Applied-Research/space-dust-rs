//! Mathematical operations for astrodynamics calculations.
//!
//! This module provides 3D vector and 3x3 matrix types optimized for
//! coordinate transformations and orbital mechanics calculations.

use std::ops::{Add, Div, Mul, Neg, Sub};

// ============================================================================
// Vector3 - 3D Vector
// ============================================================================

/// A 3D vector for position, velocity, and other spatial quantities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3 {
    /// X component
    pub x: f64,
    /// Y component
    pub y: f64,
    /// Z component
    pub z: f64,
}

impl Vector3 {
    /// Create a new Vector3.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Create a zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Create a unit vector along the X axis.
    #[inline]
    pub fn unit_x() -> Self {
        Self::new(1.0, 0.0, 0.0)
    }

    /// Create a unit vector along the Y axis.
    #[inline]
    pub fn unit_y() -> Self {
        Self::new(0.0, 1.0, 0.0)
    }

    /// Create a unit vector along the Z axis.
    #[inline]
    pub fn unit_z() -> Self {
        Self::new(0.0, 0.0, 1.0)
    }

    /// Create a Vector3 from an array.
    #[inline]
    pub fn from_array(arr: [f64; 3]) -> Self {
        Self::new(arr[0], arr[1], arr[2])
    }

    /// Create a Vector3 from a tuple.
    #[inline]
    pub fn from_tuple(t: (f64, f64, f64)) -> Self {
        Self::new(t.0, t.1, t.2)
    }

    /// Convert to an array.
    #[inline]
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Convert to a tuple.
    #[inline]
    pub fn to_tuple(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }

    /// Compute the magnitude (length) of the vector.
    #[inline]
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Compute the squared magnitude (avoids sqrt).
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Normalize the vector to unit length.
    /// Returns a zero vector if the magnitude is zero.
    #[inline]
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();
        if mag > 0.0 {
            Self::new(self.x / mag, self.y / mag, self.z / mag)
        } else {
            Self::zero()
        }
    }

    /// Compute the dot product with another vector.
    #[inline]
    pub fn dot(&self, other: &Vector3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Compute the cross product with another vector.
    #[inline]
    pub fn cross(&self, other: &Vector3) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Compute the angle between two vectors in radians.
    #[inline]
    pub fn angle(&self, other: &Vector3) -> f64 {
        let dot = self.dot(other);
        let mag_product = self.magnitude() * other.magnitude();
        if mag_product > 0.0 {
            (dot / mag_product).clamp(-1.0, 1.0).acos()
        } else {
            0.0
        }
    }

    /// Scale the vector by a scalar.
    #[inline]
    pub fn scale(&self, scalar: f64) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }

    /// Rotate the vector around an arbitrary axis by an angle (radians).
    pub fn rotate_around_axis(&self, axis: &Vector3, angle: f64) -> Self {
        let axis = axis.normalize();
        let c = angle.cos();
        let s = angle.sin();
        let t = 1.0 - c;

        let x = (t * axis.x * axis.x + c) * self.x
            + (t * axis.x * axis.y - s * axis.z) * self.y
            + (t * axis.x * axis.z + s * axis.y) * self.z;

        let y = (t * axis.x * axis.y + s * axis.z) * self.x
            + (t * axis.y * axis.y + c) * self.y
            + (t * axis.y * axis.z - s * axis.x) * self.z;

        let z = (t * axis.x * axis.z - s * axis.y) * self.x
            + (t * axis.y * axis.z + s * axis.x) * self.y
            + (t * axis.z * axis.z + c) * self.z;

        Self::new(x, y, z)
    }

    /// Rotate the vector around the X axis.
    #[inline]
    pub fn rotate_x(&self, angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(self.x, c * self.y - s * self.z, s * self.y + c * self.z)
    }

    /// Rotate the vector around the Y axis.
    #[inline]
    pub fn rotate_y(&self, angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(c * self.x + s * self.z, self.y, -s * self.x + c * self.z)
    }

    /// Rotate the vector around the Z axis.
    #[inline]
    pub fn rotate_z(&self, angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(c * self.x - s * self.y, s * self.x + c * self.y, self.z)
    }

    /// Linear interpolation between two vectors.
    #[inline]
    pub fn lerp(&self, other: &Vector3, t: f64) -> Self {
        Self::new(
            self.x + t * (other.x - self.x),
            self.y + t * (other.y - self.y),
            self.z + t * (other.z - self.z),
        )
    }

    /// Check if the vector is approximately zero.
    #[inline]
    pub fn is_zero(&self, epsilon: f64) -> bool {
        self.magnitude_squared() < epsilon * epsilon
    }

    /// Check if two vectors are approximately equal.
    #[inline]
    pub fn approx_eq(&self, other: &Vector3, epsilon: f64) -> bool {
        (self.x - other.x).abs() < epsilon
            && (self.y - other.y).abs() < epsilon
            && (self.z - other.z).abs() < epsilon
    }
}

impl Default for Vector3 {
    fn default() -> Self {
        Self::zero()
    }
}

impl Add for Vector3 {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl Sub for Vector3 {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Neg for Vector3 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vector3 {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: f64) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Mul<Vector3> for f64 {
    type Output = Vector3;

    #[inline]
    fn mul(self, vec: Vector3) -> Vector3 {
        Vector3::new(self * vec.x, self * vec.y, self * vec.z)
    }
}

impl Div<f64> for Vector3 {
    type Output = Self;

    #[inline]
    fn div(self, scalar: f64) -> Self {
        Self::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

impl From<[f64; 3]> for Vector3 {
    fn from(arr: [f64; 3]) -> Self {
        Self::from_array(arr)
    }
}

impl From<(f64, f64, f64)> for Vector3 {
    fn from(t: (f64, f64, f64)) -> Self {
        Self::from_tuple(t)
    }
}

impl From<Vector3> for [f64; 3] {
    fn from(v: Vector3) -> Self {
        v.to_array()
    }
}

impl From<Vector3> for (f64, f64, f64) {
    fn from(v: Vector3) -> Self {
        v.to_tuple()
    }
}

// ============================================================================
// Matrix3 - 3x3 Matrix
// ============================================================================

/// A 3x3 matrix for rotations and coordinate transformations.
///
/// The matrix is stored in row-major order:
/// ```text
/// | m11 m12 m13 |
/// | m21 m22 m23 |
/// | m31 m32 m33 |
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3 {
    pub m11: f64,
    pub m12: f64,
    pub m13: f64,
    pub m21: f64,
    pub m22: f64,
    pub m23: f64,
    pub m31: f64,
    pub m32: f64,
    pub m33: f64,
}

impl Matrix3 {
    /// Create a new Matrix3 from individual elements.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        m11: f64,
        m12: f64,
        m13: f64,
        m21: f64,
        m22: f64,
        m23: f64,
        m31: f64,
        m32: f64,
        m33: f64,
    ) -> Self {
        Self {
            m11,
            m12,
            m13,
            m21,
            m22,
            m23,
            m31,
            m32,
            m33,
        }
    }

    /// Create an identity matrix.
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    }

    /// Create a zero matrix.
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }

    /// Create a matrix from row vectors.
    pub fn from_rows(row1: Vector3, row2: Vector3, row3: Vector3) -> Self {
        Self::new(
            row1.x, row1.y, row1.z, row2.x, row2.y, row2.z, row3.x, row3.y, row3.z,
        )
    }

    /// Create a matrix from column vectors.
    pub fn from_cols(col1: Vector3, col2: Vector3, col3: Vector3) -> Self {
        Self::new(
            col1.x, col2.x, col3.x, col1.y, col2.y, col3.y, col1.z, col2.z, col3.z,
        )
    }

    /// Create a rotation matrix about the X axis.
    pub fn rotation_x(angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c)
    }

    /// Create a rotation matrix about the Y axis.
    pub fn rotation_y(angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c)
    }

    /// Create a rotation matrix about the Z axis.
    /// This rotates vectors counter-clockwise when viewed from +Z axis.
    pub fn rotation_z(angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0)
    }

    /// Get a row as a Vector3.
    pub fn row(&self, index: usize) -> Vector3 {
        match index {
            0 => Vector3::new(self.m11, self.m12, self.m13),
            1 => Vector3::new(self.m21, self.m22, self.m23),
            2 => Vector3::new(self.m31, self.m32, self.m33),
            _ => panic!("Row index out of bounds"),
        }
    }

    /// Get a column as a Vector3.
    pub fn col(&self, index: usize) -> Vector3 {
        match index {
            0 => Vector3::new(self.m11, self.m21, self.m31),
            1 => Vector3::new(self.m12, self.m22, self.m32),
            2 => Vector3::new(self.m13, self.m23, self.m33),
            _ => panic!("Column index out of bounds"),
        }
    }

    /// Transpose the matrix.
    pub fn transpose(&self) -> Self {
        Self::new(
            self.m11, self.m21, self.m31, self.m12, self.m22, self.m32, self.m13, self.m23,
            self.m33,
        )
    }

    /// Compute the determinant.
    pub fn determinant(&self) -> f64 {
        self.m11 * (self.m22 * self.m33 - self.m23 * self.m32)
            - self.m12 * (self.m21 * self.m33 - self.m23 * self.m31)
            + self.m13 * (self.m21 * self.m32 - self.m22 * self.m31)
    }

    /// Compute the inverse of the matrix.
    /// Returns None if the matrix is singular.
    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det.abs() < 1e-15 {
            return None;
        }

        let inv_det = 1.0 / det;

        Some(Self::new(
            (self.m22 * self.m33 - self.m23 * self.m32) * inv_det,
            (self.m13 * self.m32 - self.m12 * self.m33) * inv_det,
            (self.m12 * self.m23 - self.m13 * self.m22) * inv_det,
            (self.m23 * self.m31 - self.m21 * self.m33) * inv_det,
            (self.m11 * self.m33 - self.m13 * self.m31) * inv_det,
            (self.m13 * self.m21 - self.m11 * self.m23) * inv_det,
            (self.m21 * self.m32 - self.m22 * self.m31) * inv_det,
            (self.m12 * self.m31 - self.m11 * self.m32) * inv_det,
            (self.m11 * self.m22 - self.m12 * self.m21) * inv_det,
        ))
    }

    /// Multiply the matrix by a vector.
    pub fn mul_vec(&self, v: &Vector3) -> Vector3 {
        Vector3::new(
            self.m11 * v.x + self.m12 * v.y + self.m13 * v.z,
            self.m21 * v.x + self.m22 * v.y + self.m23 * v.z,
            self.m31 * v.x + self.m32 * v.y + self.m33 * v.z,
        )
    }

    /// Multiply two matrices.
    pub fn mul_mat(&self, other: &Matrix3) -> Self {
        Self::new(
            self.m11 * other.m11 + self.m12 * other.m21 + self.m13 * other.m31,
            self.m11 * other.m12 + self.m12 * other.m22 + self.m13 * other.m32,
            self.m11 * other.m13 + self.m12 * other.m23 + self.m13 * other.m33,
            self.m21 * other.m11 + self.m22 * other.m21 + self.m23 * other.m31,
            self.m21 * other.m12 + self.m22 * other.m22 + self.m23 * other.m32,
            self.m21 * other.m13 + self.m22 * other.m23 + self.m23 * other.m33,
            self.m31 * other.m11 + self.m32 * other.m21 + self.m33 * other.m31,
            self.m31 * other.m12 + self.m32 * other.m22 + self.m33 * other.m32,
            self.m31 * other.m13 + self.m32 * other.m23 + self.m33 * other.m33,
        )
    }

    /// Scale the matrix by a scalar.
    pub fn scale(&self, scalar: f64) -> Self {
        Self::new(
            self.m11 * scalar,
            self.m12 * scalar,
            self.m13 * scalar,
            self.m21 * scalar,
            self.m22 * scalar,
            self.m23 * scalar,
            self.m31 * scalar,
            self.m32 * scalar,
            self.m33 * scalar,
        )
    }

    /// Check if this is approximately an identity matrix.
    pub fn is_identity(&self, epsilon: f64) -> bool {
        (self.m11 - 1.0).abs() < epsilon
            && (self.m22 - 1.0).abs() < epsilon
            && (self.m33 - 1.0).abs() < epsilon
            && self.m12.abs() < epsilon
            && self.m13.abs() < epsilon
            && self.m21.abs() < epsilon
            && self.m23.abs() < epsilon
            && self.m31.abs() < epsilon
            && self.m32.abs() < epsilon
    }

    /// Check if this is approximately a rotation matrix (orthogonal with det = 1).
    pub fn is_rotation(&self, epsilon: f64) -> bool {
        // Check determinant is approximately 1
        if (self.determinant() - 1.0).abs() > epsilon {
            return false;
        }

        // Check orthogonality: M * M^T should be identity
        let product = self.mul_mat(&self.transpose());
        product.is_identity(epsilon)
    }
}

impl Default for Matrix3 {
    fn default() -> Self {
        Self::identity()
    }
}

impl Mul<Matrix3> for Matrix3 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mul_mat(&other)
    }
}

impl Mul<Vector3> for Matrix3 {
    type Output = Vector3;

    fn mul(self, v: Vector3) -> Vector3 {
        self.mul_vec(&v)
    }
}

impl Mul<&Vector3> for Matrix3 {
    type Output = Vector3;

    fn mul(self, v: &Vector3) -> Vector3 {
        self.mul_vec(v)
    }
}

impl Mul<f64> for Matrix3 {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self {
        self.scale(scalar)
    }
}

// ============================================================================
// Mathematical Functions
// ============================================================================

/// Evaluate a polynomial using Horner's method.
/// Coefficients are ordered from highest degree to constant term.
/// e.g., `[a, b, c, d]` represents `a*x³ + b*x² + c*x + d`
pub fn poly_eval(coefficients: &[f64], x: f64) -> f64 {
    coefficients.iter().fold(0.0, |acc, &coef| acc * x + coef)
}

/// Linear interpolation between two values.
#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// Clamp a value to a range.
#[inline]
pub fn clamp(value: f64, min: f64, max: f64) -> f64 {
    value.max(min).min(max)
}

/// Sign function: returns -1, 0, or 1 based on the sign of the input.
#[inline]
pub fn sign(value: f64) -> f64 {
    if value > 0.0 {
        1.0
    } else if value < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Safe arcsin that clamps input to [-1, 1].
#[inline]
pub fn safe_asin(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).asin()
}

/// Safe arccos that clamps input to [-1, 1].
#[inline]
pub fn safe_acos(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).acos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPSILON: f64 = 1e-12;

    #[test]
    fn test_vector3_basic() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v.x, 1.0);
        assert_eq!(v.y, 2.0);
        assert_eq!(v.z, 3.0);
    }

    #[test]
    fn test_vector3_magnitude() {
        let v = Vector3::new(3.0, 4.0, 0.0);
        assert!((v.magnitude() - 5.0).abs() < EPSILON);
    }

    #[test]
    fn test_vector3_normalize() {
        let v = Vector3::new(3.0, 4.0, 0.0);
        let n = v.normalize();
        assert!((n.magnitude() - 1.0).abs() < EPSILON);
        assert!((n.x - 0.6).abs() < EPSILON);
        assert!((n.y - 0.8).abs() < EPSILON);
    }

    #[test]
    fn test_vector3_dot() {
        let v1 = Vector3::new(1.0, 2.0, 3.0);
        let v2 = Vector3::new(4.0, 5.0, 6.0);
        assert!((v1.dot(&v2) - 32.0).abs() < EPSILON);
    }

    #[test]
    fn test_vector3_cross() {
        let x = Vector3::unit_x();
        let y = Vector3::unit_y();
        let z = x.cross(&y);
        assert!(z.approx_eq(&Vector3::unit_z(), EPSILON));
    }

    #[test]
    fn test_vector3_rotate_z() {
        let v = Vector3::unit_x();
        let rotated = v.rotate_z(PI / 2.0);
        assert!(rotated.approx_eq(&Vector3::unit_y(), EPSILON));
    }

    #[test]
    fn test_matrix3_identity() {
        let m = Matrix3::identity();
        let v = Vector3::new(1.0, 2.0, 3.0);
        let result = m.mul_vec(&v);
        assert!(result.approx_eq(&v, EPSILON));
    }

    #[test]
    fn test_matrix3_rotation_z() {
        let m = Matrix3::rotation_z(PI / 2.0);
        let v = Vector3::unit_x();
        let result = m.mul_vec(&v);
        assert!(result.approx_eq(&Vector3::unit_y(), EPSILON));
    }

    #[test]
    fn test_matrix3_transpose() {
        let m = Matrix3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        let t = m.transpose();
        assert!((t.m12 - 4.0).abs() < EPSILON);
        assert!((t.m21 - 2.0).abs() < EPSILON);
    }

    #[test]
    fn test_matrix3_determinant() {
        let m = Matrix3::identity();
        assert!((m.determinant() - 1.0).abs() < EPSILON);

        let m2 = Matrix3::rotation_z(0.5);
        assert!((m2.determinant() - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_matrix3_inverse() {
        let m = Matrix3::rotation_z(0.5);
        let inv = m.inverse().unwrap();
        let product = m.mul_mat(&inv);
        assert!(product.is_identity(EPSILON));
    }

    #[test]
    fn test_matrix3_is_rotation() {
        let m = Matrix3::rotation_z(0.5);
        assert!(m.is_rotation(EPSILON));

        let m2 = Matrix3::identity().scale(2.0);
        assert!(!m2.is_rotation(EPSILON));
    }

    #[test]
    fn test_poly_eval() {
        // 2x² + 3x + 1 at x = 2 should be 2*4 + 3*2 + 1 = 15
        let coeffs = [2.0, 3.0, 1.0];
        assert!((poly_eval(&coeffs, 2.0) - 15.0).abs() < EPSILON);
    }

    #[test]
    fn test_lerp() {
        assert!((lerp(0.0, 10.0, 0.5) - 5.0).abs() < EPSILON);
        assert!((lerp(0.0, 10.0, 0.0) - 0.0).abs() < EPSILON);
        assert!((lerp(0.0, 10.0, 1.0) - 10.0).abs() < EPSILON);
    }
}
