use core::f32;
use std::vec;

pub fn euler(ode: &dyn Fn(f32, f32) -> f32, initial_value: f32, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<f32> {

    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<f32> = vec![initial_value];

    for ii in 1..point_number{
        let i = ii as f32;
        let iu: usize = ii as usize;

        let derivative: f32 = ode(solution[iu-1], i*h); //first order ODEs only
        solution.push(solution[iu-1] + derivative*h)
    }

    solution
}

pub fn euler_system(ode_system: Vec<&dyn Fn(&Vec<f32>, f32) -> f32>, initial_values: Vec<f32>, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<Vec<f32>> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<Vec<f32>> = vec![initial_values];
    let order: usize = ode_system.len();

    for ii in 1..point_number{
        let i = ii as f32;
        let iu: usize = ii as usize;

        let prior_state: Vec<f32> = solution[iu-1].clone(); //prior estimates for all orders; initial state if this is i=1
        solution.push(vec![0.0; order]); //Add a new vector for this timestep's orders
        for o in 0..order {
            let derivative: f32 = ode_system[o](&prior_state, i*h); //for each timestep, generate the derivative for each order from all prior values
            solution[iu][o] = prior_state[o] + derivative * h; //generate the next step from the derivative for each order
        }
    }

    solution
}

pub fn rk5(ode: &dyn Fn(f32, f32) -> f32, initial_value: f32, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<f32> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<f32> = vec![initial_value];

    for ii in 1..point_number{
        let i = ii as f32;
        let iu = ii as usize;

        let prior: f32 = solution[iu-1];
        let k1: f32 = ode(prior, i*h);
        let k2: f32 = ode(prior, i*h + (h/2.0)) + h*(k1/2.0);
        let k3: f32 = ode(prior, i*h + (h/2.0)) + h*(k2/2.0);
        let k4: f32 = ode(prior, i*h + h) + h*k3;

        solution.push(prior + ((h/6.0) *(k1 + (2.0*k2) + (2.0*k3) + k4)));
    }

    solution
}

pub fn rk5_system(ode_system: Vec<&dyn Fn(&Vec<f32>, f32) -> f32>, initial_values: Vec<f32>, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<Vec<f32>> {
    let point_number:i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<Vec<f32>> = vec![initial_values];
    let order: usize = ode_system.len();

    for ii in 1..point_number{
        let i: f32 = ii as f32;
        let iu: usize = ii as usize;

        let prior_state: Vec<f32> = solution[iu-1].clone(); //prior estimates for all orders; initial state if this is i=1
        solution.push(vec![0.0; order]); //Add a new vector for this timestep's orders
        for o in 0..order {
            let k1: f32 = ode_system[o](&prior_state, i*h);
            let k2: f32 = ode_system[o](&prior_state, i*h + (h/2.0)) + h*(k1/2.0);
            let k3: f32 = ode_system[o](&prior_state, i*h + (h/2.0)) + h*(k2/2.0);
            let k4: f32 = ode_system[o](&prior_state, i*h + h) + h*k3;
            
            solution[iu][o] = prior_state[o] + ((h/6.0) *(k1 + (2.0*k2) + (2.0*k3) + k4));
        }
    }

    solution
}