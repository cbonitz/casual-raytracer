extern crate image;

extern crate nalgebra as na;
use na::Vector3;
use std::f64;

pub struct Sphere {
    pub origin: Vector3<f64>,
    pub radius: f64,
}

impl Sphere {
    pub fn new(x: f64, y: f64, z: f64, radius: f64) -> Sphere {
        Sphere {
            origin: Vector3::new(x, y, z),
            radius,
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        let dot_product = (self.origin - ray.origin).dot(&ray.direction);
        if dot_product < 0.0 {
            // ray pointing away from sphere
            return None;
        }
        // for the following calculations, we don't take roots if not necessary

        // use orthogonal decomposition to find how close to the origin the ray comes
        let closest_to_origin_squared =
            (self.origin - ray.origin).norm_squared() - dot_product.powi(2);
        let square_radius = self.radius.powi(2);
        // if this is further than the radius, there is no intersection
        if closest_to_origin_squared > square_radius {
            return None;
        }

        // we have an intersection. calculate values for Intersection struct
        let intersection_distance =
            dot_product - (square_radius - closest_to_origin_squared).sqrt();
        Some(Intersection {
            distance: intersection_distance
        })
    }

    pub fn surface_normal(&self, intersection_point: Vector3<f64>) -> Vector3<f64> {
        (intersection_point - self.origin).normalize()
    }
}

pub struct Intersection {
    pub distance: f64
}

pub struct Ray {
    pub origin: Vector3<f64>,
    pub direction: Vector3<f64>,
}

impl Ray {
    pub fn from_to(origin: Vector3<f64>, target: &Vector3<f64>) -> Ray {
        Ray {
            origin,
            direction: (target - origin).normalize(),
        }
    }
}

fn main() {
    let origin = Vector3::<f64>::new(0.0, 0.0, 0.0);
    let film_distance: f64 = 1000.0;
    let dimensions = (1920 as usize, 1200 as usize);
    let spheres = vec![
        Sphere::new(-200.0, -200.0, 7000.0, 500.0),
        Sphere::new(200.0, 200.0, 5000.0, 500.0)
    ];
    let lights = vec![
        Vector3::<f64>::new(1000.0, 1000.0, 1000.0),
        Vector3::<f64>::new(0.0, 0.0, 0.0),
    ];
    let upper_left = (-(dimensions.0 as f64) / 2.0, -(dimensions.1 as f64) / 2.0);

    let number_of_pixels = dimensions.0 * dimensions.1;
    let mut pixels = vec![0.0; number_of_pixels];
    for (i, pixel) in pixels.iter_mut().enumerate() {
        let x = (i % dimensions.0) as f64;
        let y = (i / dimensions.0) as f64;
        let mut closest_match = f64::MAX;
        let ray = Ray::from_to(
            origin,
            &Vector3::<f64>::new(upper_left.0 + x, upper_left.1 + y, film_distance),
        );
        *pixel = 0.0;
        for sphere in &spheres {
            match sphere.intersect(&ray) {
                None => {}
                Some(Intersection { distance }) => {
                    if distance < closest_match {
                        let intersection_point = ray.origin + ray.direction * distance;
                        let mut illumination = 0.0;
                        for light in &lights {
                            let shadow_ray = Ray::from_to(intersection_point, light);
                            let mut occluded = false;
                            for other_sphere in &spheres {
                                if other_sphere as *const _ != sphere as *const _ {
                                    if let Some(_) = other_sphere.intersect(&shadow_ray) {  
                                        occluded = true;
                                        break;
                                    }
                                }
                            }
                            if occluded {
                                continue;
                            }
                            let surface_normal = sphere.surface_normal(intersection_point);
                            let light_intensity_dot = surface_normal.dot(&shadow_ray.direction);
                            if  light_intensity_dot > 0.0 {
                                illumination += light_intensity_dot / (lights.len() as f64);
                            }
                        }
                        *pixel = illumination;
                        closest_match = distance;
                    }
                }
            }
        }
    }
    let mut imgbuf = image::GrayImage::new(dimensions.0 as u32, dimensions.1 as u32);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let array_index = x as usize + (y as usize * dimensions.0);
        *pixel = image::Luma([(pixels[array_index] * 255.0) as u8]);
    }
    imgbuf.save("render.png").unwrap();
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {}
}
