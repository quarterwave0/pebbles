fn main() {

    let integral_soltuion_1d: f32 = trapz(&fint, -1.0, 1.0, 0.01);
    let integral_solution_2d: f32 = trapz_2(&fint2, -4.0, -4.0, 4.0, 4.0, 0.01, 0.01);

    println!("Evaluated 1D: {0}", integral_soltuion_1d);
    println!("Evaluated 2D: {0}", integral_solution_2d);
}

fn fint(x: f32) -> f32 {
    x.powf(2.0)
}

fn fint2(x: f32, y:f32) -> f32 {
    (x.powf(2.0))*(y.powf(2.0))
}

fn trapz(eval_function: &dyn Fn(f32) -> f32, lower_bound: f32, upper_bound: f32, step_size: f32) -> f32 {
    let n_steps: i32 = ((upper_bound-lower_bound)/step_size) as i32;
    let mut intg: f32 = 0.0;

    for i in 0..n_steps {

        let trapz_left: f32 = (i as f32 * step_size) + lower_bound;
        let trapz_right: f32 = ((i as f32 + 1.0) * step_size) + lower_bound;

        intg += ((eval_function(trapz_left)+eval_function(trapz_right)) / 2.0) * step_size;
    }

    intg
}

fn trapz_2(eval_function: &dyn Fn(f32, f32) -> f32, lower_bound_a: f32, lower_bound_b: f32, upper_bound_a: f32, upper_bound_b:f32, step_size_a: f32, step_size_b:f32) -> f32 {
    let mut intg: f32 = 0.0;
    let n_steps_b: i32 = ((upper_bound_b-lower_bound_b)/step_size_b) as i32;

    for i in 0..n_steps_b { // like \int\int f(x,y) dydx, we are doing dy with trapz and dx right here

        let trapz_left_b: f32 = (i as f32 * step_size_b) + lower_bound_b;
        let trapz_right_b: f32 = ((i as f32 + 1.0) * step_size_b) + lower_bound_b;
        let dim_reduction_left = |a: f32| {eval_function(a, trapz_left_b)};
        let dim_reduction_right = |a: f32| {eval_function(a, trapz_right_b)};


        intg += ((trapz(&dim_reduction_left, lower_bound_a, upper_bound_a, step_size_a) + trapz(&dim_reduction_right, lower_bound_a, upper_bound_a, step_size_a)) / 2.0) * step_size_a;
    }
    intg
}

