mod trapezoidal;

fn main() {

    let integral_soltuion_1d: f32 = trapezoidal::trapz(&fint, &-1.0, &1.0, &0.01, false);
    let integral_solution_2d: f32 = trapezoidal::trapz_2(&fint2, &-4.0, &-4.0, &4.0, &4.0, &0.01, &0.01);

    println!("Evaluated 1D: {0}", integral_soltuion_1d);
    println!("Evaluated 2D: {0}", integral_solution_2d);
}

fn fint(x: &f32) -> f32 {
    x.powf(2.0)
}

fn fint2(x: &f32, y: &f32) -> f32 {
    (x.powf(2.0))*(y.powf(2.0))
}

