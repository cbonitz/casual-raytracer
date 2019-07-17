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

/// A point light source
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

/// Types of surfaces.
/// Currently diffuse and perfectly specular are implemented.
pub enum Surface {
    Diffuse { color: Vector3<f32> },
    Specular,
}

pub struct Sphere {
    pub origin: Vector3<f32>,
    pub radius: f32,
    pub surface: Surface,
}

/// A sphere, defined by an origin, its radius and its surface.
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

    pub fn surface_normal(&self, intersection_point: &Vector3<f32>) -> Vector3<f32> {
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

/// The colour to draw if a ray goes to infinity.
///
/// A very useful debugging tool from Peter Shirley's "Ray Tracing in One Weekend".
fn infinite_color(ray: &Ray) -> Vector3<f32> {
    let background_color = 0.5 * (ray.direction[1] + 1.0) as f32;
    let result = Vector3::new(1.0 - background_color, 1.0 - background_color, 1.0);
    result
}

/// Recursive function to compute the color to draw for a given ray and a given scenery,
/// at recursion depth `depth`.
///
/// If `depth` is greater than or equal to `MAX_DEPTH`, we return the color at infinity.
fn trace_ray(ray: &Ray, spheres: &Vec<Sphere>, lights: &Vec<Light>, depth: u32) -> Vector3<f32> {
    let mut closest_match = f32::MAX;
    let mut closest_sphere: Option<&Sphere> = None;

    if depth >= MAX_DEPTH {
        return infinite_color(&ray);
    }
    // find the closest intersection, if any
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
            // Compute an intersection point. Doing this here avoids unnecessary alloctions.
            let intersection_point = ray.origin + ray.direction * closest_match;
            match sphere.surface {
                // Diffuse surface
                // Find all unoccluded light sources and add up their light contributions.
                Surface::Diffuse { color } => {
                    let mut light_intensity: f32 = 0.0;
                    for light in lights {
                        // Use a shadow ray to check whether a light source is visible:
                        // If the ray between the intersection point and a light source does not
                        // intersect other scene object before the light source, the light source is
                        // unoccluded.
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
                        let surface_normal = sphere.surface_normal(&intersection_point);
                        let light_intensity_scaling_factor =
                            surface_normal.dot(&shadow_ray.direction);
                        if light_intensity_scaling_factor > 0.0 {
                            light_intensity += light.intensity * light_intensity_scaling_factor;
                        }
                    }
                    color * light_intensity
                }
                // Specular surface: Recurse using a reflected ray.
                Surface::Specular => {
                    let surface_normal = sphere.surface_normal(&intersection_point);
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
        // no intersection, return color at infinity
        None => infinite_color(&ray),
    }
}

/// Render an image of size `dimensions` depicting the spheres and
/// lights, projecting rays from the origin of the coordinate frame
/// through an image plane at x=y=0, z=1, with the left edge being at
/// (-1, 0, 1), the right one at (1, 0, 1), and height determined by
/// the aspect ratio of `dimensions`.
/// 
/// The result is a 2-D Vec of dimension (number of pixels, 3),
/// with the first dimension being the pixels in top-left to 
/// bottom-right scanline order, and the second dimensions being
/// the colors (RGB)
pub fn render(
    dimensions: (usize, usize),
    spheres: &Vec<Sphere>,
    lights: &Vec<Light>,
) -> Vec<Vector3<f32>> {
    // Hardcoded camero parameters
    let film_distance: f32 = 1.0;
    let origin = Vector3::<f32>::new(0.0, 0.0, 0.0);
    // Using the aspect ratio of the image, compute the distance between
    // pixel centers as well as the upper and lower edges of the image
    // in world coordinates.
    let inverse_aspect_ratio = dimensions.1 as f32 / dimensions.0 as f32;
    let delta_x = 2.0 / dimensions.0 as f32;
    let delta_y = 2.0 * inverse_aspect_ratio / dimensions.1 as f32;
    let upper_left = (-1.0, inverse_aspect_ratio);

    // Create indices for the pixels as well as the subpixel offsets
    // to feed into the actual ray tracing loop. We use fixed subpixel
    // offsets as randomization didn't seem to provide any benefits but
    // was slower
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
    
    // Compute pixels in parallel...
    indices
        .par_iter()
        .map(|i| -> Vector3<f32> {
            // ... and subpixels within a thread.
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

/// Expose the image, currently using scaling to [0, 1]
pub fn expose(pixels: &Vec<Vector3<f32>>) -> Vec<Vector3<f32>> {
    let scaling_factor = 1.0 / pixels.iter().map(|x| x[x.imax()]).fold(0.0, f32::max);
    pixels.iter().map(|x| x * scaling_factor).collect()
}

/// Save an image to a `filename`.
pub fn save(
    filename: &str,
    dimensions: (usize, usize),
    pixels: &Vec<Vector3<f32>>,
) -> std::io::Result<()> {
    let mut imgbuf = image::RgbImage::new(dimensions.0 as u32, dimensions.1 as u32);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let array_index = x as usize + (y as usize * dimensions.0);
        let channels = pixels[array_index].map(|x| (x * 255.9) as u8);
        for (i, channel) in channels.iter().enumerate() {
            (*pixel)[i] = *channel;
        }
    }
    imgbuf.save(filename)
}

#[cfg(test)]
mod tests {
    use super::*;

    // let's be generous for now
    const EPSILON: f32 = 0.0001;

    /// Ray intersecting a sphere in an easy to compute position.
    #[test]
    fn ray_intersects_sphere() {
        let origin = Vector3::<f32>::new(0.0, 0.0, 0.0);
        let sphere_origin = Vector3::<f32>::new(0.0, 0.0, 1.0);
        let sphere_radius = 0.5;
        let ray = Ray::from_to(origin, &sphere_origin);
        let sphere = Sphere::new(
            sphere_origin[0],
            sphere_origin[1],
            sphere_origin[2],
            sphere_radius,
            Surface::Specular,
        );
        let intersection = sphere.intersect(&ray).unwrap();
        assert!((intersection.distance - 0.5).abs() < EPSILON);
    }
    /// Trivial case of a surface normal calculation.
    #[test]
    fn sphere_surface_normal() {
        let sphere = Sphere::new(0.0, 0.0, 1.0, 0.1, Surface::Specular);
        let intersection_point = Vector3::<f32>::new(0.0, 0.0, 0.0);
        let surface_normal = sphere.surface_normal(&intersection_point);
        let expected_surface_normal = Vector3::<f32>::new(0.0, 0.0, -1.0);
        for i in 0..3 {
            assert!((surface_normal[i] - expected_surface_normal[i]).abs() < EPSILON);
        }
    }
}
