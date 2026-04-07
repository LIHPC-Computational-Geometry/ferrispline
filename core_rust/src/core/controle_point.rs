use std::ops::{Add, Mul};

// core/src/geometry.rs
use nalgebra::{Point3, Point4};

pub trait Point3Ext {
    fn add_points(&self, rhs: &Self) -> Self;
}

impl Point3Ext for Point3<f64> {
    fn add_points(&self, rhs: &Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

#[allow(dead_code)]
pub trait ControlPoint: Copy {
    type Homogeneous;

    fn new(x: f64, y: f64, z: f64) -> Self;

    fn to_homogeneous(&self) -> Self::Homogeneous;
    fn from_homogeneous(h: &Self::Homogeneous) -> Self;
}

impl ControlPoint for Point3<f64> {
    type Homogeneous = Point4<f64>;
    fn new(x: f64, y: f64, z: f64) -> Self {
        Point3::new(x, y, z)
    }

    fn to_homogeneous(&self) -> Self::Homogeneous {
        Point4::new(self.x, self.y, self.z, 1.0)
    }

    fn from_homogeneous(h: &Self::Homogeneous) -> Self {
        Point3::new(h.x / h.w, h.y / h.w, h.z / h.w)
    }
}

#[derive(Copy, Clone)]
pub struct RationalPoint3 {
    pub point: Point3<f64>,
    pub weight: f64,
}

impl RationalPoint3 {
    #[allow(dead_code)]
    fn weight(&mut self, weight: f64) -> f64 {
        self.weight = weight;
        self.weight
    }
}

impl ControlPoint for RationalPoint3 {
    type Homogeneous = Point4<f64>;

    fn new(x: f64, y: f64, z: f64) -> Self {
        RationalPoint3 {
            point: Point3::new(x, y, z),
            weight: 1.0,
        }
    }

    fn to_homogeneous(&self) -> Self::Homogeneous {
        Point4::new(
            self.point.x * self.weight,
            self.point.y * self.weight,
            self.point.z * self.weight,
            self.weight,
        )
    }

    fn from_homogeneous(h: &Self::Homogeneous) -> Self {
        RationalPoint3 {
            point: Point3::new(h.x / h.w, h.y / h.w, h.z / h.w),
            weight: h.w,
        }
    }
}

impl Add for RationalPoint3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let h1 = self.to_homogeneous();
        let h2 = rhs.to_homogeneous();
        let res = h1 + h2.coords;
        Self::from_homogeneous(&res)
    }
}

impl Mul<f64> for RationalPoint3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        RationalPoint3 {
            point: self.point * rhs,
            weight: self.weight,
        }
    }
}
