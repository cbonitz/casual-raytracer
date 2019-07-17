extern crate image;
extern crate nalgebra as na;

use na::Vector3;
use tracelib::{expose, render, save, Light, Sphere, Surface};

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

    save("render.png", dimensions, &pixels).unwrap();
}
