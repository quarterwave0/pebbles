pub fn euler(ode: &dyn Fn(&f32, &f32) -> f32, initial_value: f32, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<f32> {

    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<f32> = vec![initial_value];

    for i in 1..point_number{
        let derivative: f32 = ode(&(solution[(i-1) as usize]), &((i as f32)*h)); //first order ODEs only
        solution.push(solution[(i-1) as usize] + derivative*h)
    }

    solution
}

pub fn rk5(ode: &dyn Fn(&f32, &f32) -> f32, initial_value: f32, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<f32> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<f32> = vec![initial_value];

    for i in 1..point_number{
        let prior: f32 = solution[(i-1) as usize];

        let k1: f32 = ode(&(prior), &((i as f32)*h));
        let k2: f32 = ode(&(prior), &((i as f32)*h + (h/2.0))) + h*(k1/2.0);
        let k3: f32 = ode(&(prior), &((i as f32)*h + (h/2.0))) + h*(k2/2.0);
        let k4: f32 = ode(&(prior), &((i as f32)*h + h)) + h*k3;

        solution.push(prior + ((h/6.0) *(k1 + (2.0*k2) + (2.0*k3) + k4)));
    }

    solution
}