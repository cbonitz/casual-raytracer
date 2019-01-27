extern crate image;
extern crate rand;

extern crate nalgebra as na;
use na::Vector3;
#[allow(unused_imports)]
use rand::{thread_rng, Rng};
use std::f32;

const MAX_DEPTH: u32 = 5;

pub enum Surface {
    Diffuse { color: Vector3<f32>},
    Specular,
}

pub struct Light {
    pub coordinates: Vector3<f32>,
    pub intensity: f32,
}

impl Light {
    pub fn new(x: f32, y: f32, z: f32, intensity: f32) -> Light {
        Light {
            coordinates: Vector3::new(x, y, z),
            intensity,
        }
    }
}

pub struct Sphere {
    pub origin: Vector3<f32>,
    pub radius: f32,
    pub surface: Surface,
}

impl Sphere {
    pub fn new(x: f32, y: f32, z: f32, radius: f32, surface: Surface) -> Sphere {
        Sphere {
            origin: Vector3::new(x, y, z),
            radius,
            surface,
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
            distance: intersection_distance,
        })
    }

    pub fn surface_normal(&self, intersection_point: Vector3<f32>) -> Vector3<f32> {
        (intersection_point - self.origin).normalize()
    }
}

pub struct Intersection {
    pub distance: f32,
}

pub struct Ray {
    pub origin: Vector3<f32>,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn from_to(origin: Vector3<f32>, target: &Vector3<f32>) -> Ray {
        Ray {
            origin,
            direction: (target - origin).normalize(),
        }
    }

    pub fn origin_direction(origin: Vector3<f32>, normalized_direction: Vector3<f32>) -> Ray {
        Ray {
            origin,
            direction: normalized_direction,
        }
    }
}

fn infinite_color(ray: &Ray) -> Vector3<f32> {
    let background_color = 0.5 * (ray.direction[1] + 1.0) as f32;
    let result = Vector3::new(1.0 - background_color, 1.0 - background_color, 1.0);
    result
}

fn trace_ray(ray: &Ray, spheres: &Vec<Sphere>, lights: &Vec<Light>, depth: u32) -> Vector3<f32> {
    let mut closest_match = f32::MAX;
    let mut pixel = infinite_color(&ray);
    if depth >= MAX_DEPTH {
        return pixel;
    }
    for sphere in spheres {
        match sphere.intersect(&ray) {
            None => {}
            Some(Intersection { distance }) => {
                if distance > 0.0 && distance < closest_match {
                    closest_match = distance;
                    let intersection_point = ray.origin + ray.direction * distance;
                    match sphere.surface {
                        Surface::Diffuse { color } => {
                            let mut light_intensity: f32 = 0.0;
                            for light in lights {
                                let shadow_ray =
                                    Ray::from_to(intersection_point, &light.coordinates);
                                let light_distance =
                                    (light.coordinates - intersection_point).norm();
                                let mut occluded = false;
                                for other_sphere in spheres {
                                    if other_sphere as *const _ != sphere as *const _ {
                                        if let Some(distance) = other_sphere.intersect(&shadow_ray)
                                        {
                                            if distance.distance < light_distance {
                                                occluded = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if occluded {
                                    continue;
                                }
                                let surface_normal = sphere.surface_normal(intersection_point);
                                let light_intensity_scaling_factor = surface_normal.dot(&shadow_ray.direction);
                                if light_intensity_scaling_factor > 0.0 {
                                    light_intensity += light.intensity * light_intensity_scaling_factor;
                                }
                            }
                            pixel = color * light_intensity;
                        }
                        Surface::Specular => {
                            let surface_normal = sphere.surface_normal(intersection_point);
                            let to_source = - ray.direction;
                            let out_direction = 2.0 * to_source.dot(&surface_normal) * surface_normal - to_source;
                            pixel = trace_ray(&Ray::origin_direction(intersection_point, out_direction), spheres, lights, depth + 1);
                        }
                    }
                }
            }
        }
    }
    pixel
}

fn render(
    dimensions: (usize, usize),
    spheres: &Vec<Sphere>,
    lights: &Vec<Light>,
) -> Vec<Vector3<f32>> {
    let film_distance: f32 = 1000.0;
    let origin = Vector3::<f32>::new(0.0, 0.0, 0.0);
    let upper_left = (-(dimensions.0 as f32) / 2.0, (dimensions.1 as f32) / 2.0);
    let number_of_pixels = dimensions.0 * dimensions.1;
    let mut pixels = vec![Vector3::<f32>::new(0.0, 0.0, 0.0); number_of_pixels];
    for (i, pixel) in pixels.iter_mut().enumerate() {
        let offset_right = (i % dimensions.0) as f32;
        let offset_down = (i / dimensions.0) as f32;

        let ray = Ray::from_to(
            origin,
            &Vector3::<f32>::new(
                upper_left.0 + offset_right,
                upper_left.1 - offset_down - 1.0,
                film_distance,
            ),
        );

        *pixel = trace_ray(&ray, spheres, lights, 1);
    }
    pixels
}

fn expose(pixels: &Vec<Vector3<f32>>) -> Vec<Vector3<f32>> {
    let scaling_factor = 1.0 / pixels.iter().map(|x| x[x.imax()]).fold(0.0, f32::max);
    pixels.iter().map(|x| x * scaling_factor).collect()
}

fn main() {
    let dimensions = (1920 as usize, 1200 as usize);
    let spheres = vec![
        Sphere::new(-200.0, 200.0, 7000.0, 500.0, Surface::Diffuse { color: Vector3::<f32>::new(1.0, 1.0, 0.0) }),
        Sphere::new(200.0, -200.0, 5000.0, 500.0, Surface::Diffuse { color: Vector3::<f32>::new(1.0, 0.0, 0.0) }),
        Sphere::new(0.0, -51000.0, 5000.0, 50000.0, Surface::Diffuse { color: Vector3::<f32>::new(0.0, 1.0, 0.0) }),

        Sphere::new(2000.0, 0.0, 10000.0, 1000.0, Surface::Specular),
        Sphere::new(0.0, 0.0, -10001.0, 10000.0, Surface::Specular),
    ];
    let lights = vec![
        Light::new(1000.0, -1000.0, 1000.0, 0.5),
        Light::new(0.0, 0.0, 0.0, 0.5),
    ];

    let pixels = render(dimensions, &spheres, &lights);

    let pixels = expose(&pixels);

    let mut imgbuf = image::RgbImage::new(dimensions.0 as u32, dimensions.1 as u32);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let array_index = x as usize + (y as usize * dimensions.0);
        let channels = pixels[array_index].map(|x| (x * 255.9) as u8);
        for (i, channel) in channels.iter().enumerate() {
            (*pixel)[i] = *channel;
        }
    }
    imgbuf.save("render.png").unwrap();
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {}
}
