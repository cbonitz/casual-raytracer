extern crate image;
extern crate nalgebra as na;
extern crate rand;
extern crate rayon;

use na::Vector3;
#[allow(unused_imports)]
use rand::{thread_rng, Rng};
use rayon::prelude::*;
use std::f32;

const MAX_DEPTH: u32 = 5;
const ANTIALIASING_SQRT: usize = 5;

pub enum Surface {
    Diffuse { color: Vector3<f32> },
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
        let origin_to_ray_origin = ray.origin - self.origin;
        let b = 2.0 * origin_to_ray_origin.dot(&ray.direction);
        let c = origin_to_ray_origin.norm_squared() - self.radius * self.radius;
        // assuming unit direction
        let discriminant = b * b - 4.0 * c;
        if discriminant < 0.0 {
            return None;
        }
        let distance = 0.5 * (-b - discriminant.sqrt());
        if distance > 0.0 {
            Some(Intersection { distance })
        } else {
            None
        }
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
    let mut closest_sphere: Option<&Sphere> = None;

    if depth >= MAX_DEPTH {
        return infinite_color(&ray);
    }

    for sphere in spheres {
        match sphere.intersect(&ray) {
            None => {}
            Some(Intersection { distance }) => {
                if distance < closest_match {
                    closest_match = distance;
                    closest_sphere = Some(&sphere);
                }
            }
        }
    }
    match closest_sphere {
        Some(sphere) => {
            let intersection_point = ray.origin + ray.direction * closest_match;
            match sphere.surface {
                Surface::Diffuse { color } => {
                    let mut light_intensity: f32 = 0.0;
                    for light in lights {
                        let shadow_ray = Ray::from_to(intersection_point, &light.coordinates);
                        let light_distance = (light.coordinates - intersection_point).norm();
                        let mut occluded = false;
                        for other_sphere in spheres {
                            if other_sphere as *const _ != sphere as *const _ {
                                if let Some(distance) = other_sphere.intersect(&shadow_ray) {
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
                        let light_intensity_scaling_factor =
                            surface_normal.dot(&shadow_ray.direction);
                        if light_intensity_scaling_factor > 0.0 {
                            light_intensity += light.intensity * light_intensity_scaling_factor;
                        }
                    }
                    color * light_intensity
                }
                Surface::Specular => {
                    let surface_normal = sphere.surface_normal(intersection_point);
                    let to_source = -ray.direction;
                    let out_direction =
                        2.0 * to_source.dot(&surface_normal) * surface_normal - to_source;
                    trace_ray(
                        &Ray::origin_direction(intersection_point, out_direction),
                        spheres,
                        lights,
                        depth + 1,
                    )
                }
            }
        }
        None => infinite_color(&ray),
    }
}

fn render(
    dimensions: (usize, usize),
    spheres: &Vec<Sphere>,
    lights: &Vec<Light>,
) -> Vec<Vector3<f32>> {
    let film_distance: f32 = 1.0;
    let origin = Vector3::<f32>::new(0.0, 0.0, 0.0);
    let inverse_origin = dimensions.1 as f32 / dimensions.0 as f32;
    let delta_x = 2.0 / dimensions.0 as f32;
    let delta_y = 2.0 * inverse_origin / dimensions.1 as f32;
    let upper_left = (-1.0, inverse_origin);
    let number_of_pixels = dimensions.0 * dimensions.1;
    let indices: Vec<usize> = (0..number_of_pixels).collect();
    let antialiasing_scale_factor = 1.0 / ANTIALIASING_SQRT as f32;
    let antiliasing_step_x = delta_x / ANTIALIASING_SQRT as f32;
    let antiliasing_step_y = delta_y / ANTIALIASING_SQRT as f32;
    let subpixel_offsets_x: Vec<f32> = (0..ANTIALIASING_SQRT)
        .map(|i| 0.5 * antiliasing_step_x + i as f32 * antiliasing_step_x)
        .collect();
    let subpixel_offsets_y: Vec<f32> = (0..ANTIALIASING_SQRT)
        .map(|i| 0.5 * antiliasing_step_y + i as f32 * antiliasing_step_y)
        .collect();
    indices
        .par_iter()
        .map(|i| -> Vector3<f32> {
            let mut pixel = Vector3::<f32>::new(0.0, 0.0, 0.0);
            let offset_right = (i % dimensions.0) as f32 * delta_x;
            let offset_down = (i / dimensions.0) as f32 * delta_y;
            let pixel_origin_x = upper_left.0 + offset_right;
            let pixel_origin_y = upper_left.1 - offset_down;
            for x in 0..ANTIALIASING_SQRT {
                for y in 0..ANTIALIASING_SQRT {
                    let ray = Ray::from_to(
                        origin,
                        &Vector3::<f32>::new(
                            pixel_origin_x + subpixel_offsets_x[x],
                            pixel_origin_y - subpixel_offsets_y[y],
                            film_distance,
                        ),
                    );
                    pixel += trace_ray(&ray, spheres, lights, 1);
                }
            }
            pixel * antialiasing_scale_factor
        })
        .collect()
}

fn expose(pixels: &Vec<Vector3<f32>>) -> Vec<Vector3<f32>> {
    let scaling_factor = 1.0 / pixels.iter().map(|x| x[x.imax()]).fold(0.0, f32::max);
    pixels.iter().map(|x| x * scaling_factor).collect()
}

fn main() {
    let dimensions = (1920 as usize, 1080 as usize);
    let spheres = vec![
        Sphere::new(
            -0.2,
            0.2,
            7.0,
            0.5,
            Surface::Diffuse {
                color: Vector3::<f32>::new(1.0, 1.0, 0.0),
            },
        ),
        Sphere::new(
            0.2,
            -0.2,
            5.0,
            0.5,
            Surface::Diffuse {
                color: Vector3::<f32>::new(1.0, 0.0, 0.0),
            },
        ),
        Sphere::new(
            0.0,
            -51.0,
            5.0,
            50.0,
            Surface::Diffuse {
                color: Vector3::<f32>::new(0.0, 1.0, 0.0),
            },
        ),
        Sphere::new(2.0, 0.0, 10.0, 1.0, Surface::Specular),
        Sphere::new(0.0, 0.0, -10.1, 10.0, Surface::Specular),
    ];
    let lights = vec![
        Light::new(1.0, -1.0, 1.0, 0.5),
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
