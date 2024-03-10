mod trapezoidal;
mod finitedifference;

fn main() {

    let integral_soltuion_1d: f32 = trapezoidal::trapz(&fint, &-1.0, &1.0, &0.01, false);
    let integral_solution_2d: f32 = trapezoidal::trapz_2(&fint2, &-4.0, &-4.0, &4.0, &4.0, &0.01, &0.01);
    let differentiated_1d: Vec<f32> = finitedifference::differentiate_1d(&fint, -1.0, 1.0, 0.001, 3);
    let differentiated_2d: Vec<Vec<[f32; 2]>> = finitedifference::differentiate_2d(&fint2, -2.0, 2.0, -2.0, 2.0, 0.001, 0.001, 3);

    println!("Integral 1D: {0}", integral_soltuion_1d);
    println!("Integral 2D: {0}", integral_solution_2d);

    println!("Derivative 1D: {0}", differentiated_1d[((2.0/0.001)-1.0) as usize]);
    println!("Derivative 2D: [{0}, {1}]", differentiated_2d[3999][3999][0], differentiated_2d[3999][3999][1]);
}

fn fint(x: &f32) -> f32 {
    x.powf(2.0)
}

fn fint2(x: &f32, y: &f32) -> f32 {
    (x.powf(2.0))*(y.powf(2.0))
}