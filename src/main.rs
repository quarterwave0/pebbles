mod trapezoidal;
mod finitedifference;
mod ode;

fn main() {

    let integral_soltuion_1d: f32 = trapezoidal::trapz(&fint, -1.0, 1.0, 0.01, false);
    let integral_solution_2d: f32 = trapezoidal::trapz_2(&fint2, -4.0, -4.0, 4.0, 4.0, 0.01, 0.01);
    let differentiated_1d: Vec<f32> = finitedifference::differentiate_1d(&fint, -1.0, 1.0, 0.001, 3);
    let differentiated_2d: Vec<Vec<[f32; 2]>> = finitedifference::differentiate_2d(&fint2, -2.0, 2.0, -2.0, 2.0, 0.001, 0.001, 3);
    let euler_ode: Vec<f32> = ode::euler(&ode, 1.0, 0.001, 0.0, 3.0);
    let rk5_ode: Vec<f32> = ode::rk5(&ode, 1.0, 0.001, 0.0, 3.0);

    let euler_system_ode: Vec<Vec<f32>> = ode::euler_system(vec![&dho_ode_a, &dho_ode_b], vec![0.0, 1.0], 0.001, 0.0, 2.0);
    let rk5_system_ode: Vec<Vec<f32>> = ode::rk5_system(vec![&dho_ode_a, &dho_ode_b], vec![0.0, 1.0], 0.001, 0.0, 2.0);

    println!("Integral 1D: {0}", integral_soltuion_1d);
    println!("Integral 2D: {0}", integral_solution_2d);
    println!();
    println!("Derivative 1D: {0}", differentiated_1d[((2.0/0.001)-1.0) as usize]);
    println!("Derivative 2D: ∂f/dx={0}, ∂f/dy={1}", differentiated_2d[3999][3999][0], differentiated_2d[3999][3999][1]);
    println!();
    println!("ODE Euler 1st: {0}", euler_ode[1500]);
    println!("ODE RK5 1st: {0}", rk5_ode[1500]);
    println!("ODE Euler System: y={0}, dy/dt={1}", euler_system_ode[1000][0], euler_system_ode[1000][1]);
    println!("ODE RK5 System: y={0}, dy/dt={1}", rk5_system_ode[1000][0], rk5_system_ode[1000][1]);
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

//m*y'' + c*y' + ky = 0 for unforced dampened harmonic oscilator
// let a = y and b = y'
// a' = y' and b' = y''
// thus:
// b' = (-k*a - c*b)/m
// a' = b
// such that:
// a = x
// b = dx/dt

fn dho_ode_a(ab: &Vec<f32>, _t: f32) -> f32 {
    ab[1] //da/dt = b; dy/dt = dy/dt
}

fn dho_ode_b(ab: &Vec<f32>, _t: f32) -> f32 {
    let (m, c, k) = (1.0, 0.2, 1.0);
    (-1.0*k*ab[0] - c*ab[1]) / m //db/dt = (-ka - cb)/m; d^2y/dt^2 = (-ky - c*dy/dt)/m
}

