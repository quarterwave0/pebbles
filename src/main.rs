mod trapezoidal;
mod finitedifference;
mod ode;
mod pde;

use std::f32::consts::{E, PI};
use std::fs::File;
use std::io::Write;

fn main() {

    let integral_solution_1d: f32 = trapezoidal::trapz(&fint, -1.0, 1.0, 0.01, false);
    let integral_solution_2d: f32 = trapezoidal::trapz_2(&fint2, -4.0, -4.0, 4.0, 4.0, 0.01, 0.01);
    let differentiated_1d: Vec<f32> = finitedifference::differentiate_1d(&fint, -1.0, 1.0, 0.001, 3);
    let differentiated_2d: Vec<Vec<[f32; 2]>> = finitedifference::differentiate_2d(&fint2, -2.0, 2.0, -2.0, 2.0, 0.001, 0.001, 3);
    let euler_ode: Vec<f32> = ode::euler(&ode, 1.0, 0.001, 0.0, 3.0);
    let rk5_ode: Vec<f32> = ode::rk5(&ode, 1.0, 0.001, 0.0, 3.0);

    let euler_system_ode: Vec<Vec<f32>> = ode::euler_system(vec![&dho_ode_a, &dho_ode_b], vec![0.0, 1.0], 0.001, 0.0, 2.0);
    let rk5_system_ode: Vec<Vec<f32>> = ode::rk5_system(vec![&dho_ode_a, &dho_ode_b], vec![0.0, 1.0], 0.001, 0.0, 2.0);

    let solution_pde: Vec<Vec<f32>> = pde::ctcs2_1d(&wave_pde, &pde_icond, vec![10.0], 1.5, 0.001, 1.0, 0.0, 0.01);

    println!("Integral 1D: {0}", integral_solution_1d);
    println!("Integral 2D: {0}", integral_solution_2d);
    println!();
    println!("Derivative 1D: {0}", differentiated_1d[((2.0/0.001)-1.0) as usize]);
    println!("Derivative 2D: ∂f/∂x={0}, ∂f/∂y={1}", differentiated_2d[3999][3999][0], differentiated_2d[3999][3999][1]);
    println!();
    println!("ODE Euler 1st: {0}", euler_ode[1500]);
    println!("ODE RK5 1st: {0}", rk5_ode[1500]);
    println!("ODE Euler System: y={0}, dy/dt={1}", euler_system_ode[1000][0], euler_system_ode[1000][1]);
    println!("ODE RK5 System: y={0}, dy/dt={1}", rk5_system_ode[1000][0], rk5_system_ode[1000][1]);

    let mut file = File::create("pde_output.txt").unwrap();

    for a in solution_pde {
        for b in 0..a.len() {
            write!(file, "{:?}", a[b]).unwrap();
            if b != a.len()-1 {
                write!(file, ":").unwrap();
            }
        }
        writeln!(file).unwrap();
    }
}

fn fint(x: f32) -> f32 {
    x.powf(2.0)
}

fn fint2(x: f32, y: f32) -> f32 {
    (x.powf(2.0))*(y.powf(2.0))
}

fn ode(y: f32, t: f32) -> f32 {
    -t*y // y' + y(t)*t
}

//m*y'' + c*y' + ky = 0 for unforced dampened harmonic oscillator
// let a = y and b = y'
// a' = y' and b' = y''
// thus:
// b' = (-k*a - c*b)/m
// a' = b
// such that:
// a = y
// b = dy/dt

fn dho_ode_a(ab: &Vec<f32>, _t: f32) -> f32 {
    ab[1] //da/dt = b; dy/dt = dy/dt
}

fn dho_ode_b(ab: &Vec<f32>, _t: f32) -> f32 {
    let (m, c, k) = (1.0, 0.2, 1.0);
    (-1.0*k*ab[0] - c*ab[1]) / m //db/dt = (-ka - cb)/m; d^2y/dt^2 = (-ky - c*dy/dt)/m
}

fn wave_pde(dx: f32, constants: &Vec<f32>) -> f32 {
    return dx * constants[0].powi(2) // d^2u/dt^2 = c^2 * d^2u/dx^2
}

fn pde_icond(x: f32) -> f32 {
    //let stddev: f32 = 0.1;
    //let mean: f32 = 0.25;
    //(1.0 / (stddev * (2.0*PI).powf(0.5))) * E.powf(-0.5*((x-mean)/stddev).powi(2))

    if (x>=0.25) & (x<=0.75) {
        return E.powf(-1.0 / (1.0 - (4.0*x-2.0).powi(2))) //bump function
    }
    return 0.0
}