pub fn trapz(eval_function: &dyn Fn(f32) -> f32, lower_bound: f32, upper_bound: f32, step_size: f32, simpson: bool) -> f32 {
    let n_steps: i32 = ((upper_bound-lower_bound)/step_size) as i32;
    let mut intg: f32 = 0.0;

    for i in 0..n_steps {

        let (trapz_bound_l, trapz_bound_r) = trapezoid_lr(step_size, i, lower_bound);
        if !simpson {
            intg += ((eval_function(trapz_bound_l)+eval_function(trapz_bound_r)) / 2.0) * step_size;
        }
        else {
            intg += ((trapz_bound_r-trapz_bound_l) / 6.0) * (eval_function(trapz_bound_l) + 4.0*eval_function((trapz_bound_l+trapz_bound_r)/2.0) + eval_function(trapz_bound_r));
        }
    }

    intg
}

pub fn trapz_2(eval_function: &dyn Fn(f32, f32) -> f32, lower_bound_a: f32, lower_bound_b: f32, upper_bound_a: f32, upper_bound_b: f32, step_size_a: f32, step_size_b: f32) -> f32 {
    let mut intg: f32 = 0.0;
    let n_steps_b: i32 = ((upper_bound_b-lower_bound_b)/step_size_b) as i32;

    for i in 0..n_steps_b { // like \int\int f(x,y) dydx, we are doing dy with trapz and dx right here

        let (trapz_bound_l, trapz_bound_r) = trapezoid_lr(step_size_b, i, lower_bound_b);
        let dim_reduction_left = |a: f32| {eval_function(a, trapz_bound_l)};
        let dim_reduction_right = |a: f32| {eval_function(a, trapz_bound_r)};

        let trapezoid_l: f32 = trapz(&dim_reduction_left, lower_bound_a, upper_bound_a, step_size_a, false);
        let trapezoid_r: f32 = trapz(&dim_reduction_right, lower_bound_a, upper_bound_a, step_size_a, false);

        intg += ((trapezoid_l + trapezoid_r) / 2.0) * step_size_a;
    }
    intg
}

fn trapezoid_lr(step_size: f32, incrementor: i32, bound: f32) -> (f32, f32){
    ((incrementor as f32 * step_size) + bound, ((incrementor as f32 + 1.0) * step_size) + bound)
}