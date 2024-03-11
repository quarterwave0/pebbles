use core::f32;
use std::vec;

pub fn euler(ode: &dyn Fn(f32, f32) -> f32, initial_value: f32, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<f32> {

    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<f32> = vec![initial_value];

    for i in 1..point_number{
        let derivative: f32 = ode(solution[(i-1) as usize], (i as f32)*h); //first order ODEs only
        solution.push(solution[(i-1) as usize] + derivative*h)
    }

    solution
}

pub fn euler_system(ode_system: Vec<&dyn Fn(&Vec<f32>, f32) -> f32>, initial_values: Vec<f32>, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<Vec<f32>> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<Vec<f32>> = vec![initial_values];
    let order: i32 = ode_system.len() as i32;

    for i in 1..point_number{
        let prior_state: Vec<f32> = solution[(i-1) as usize].clone(); //prior estimates for all orders; initial state if this is i=1
        solution.push(vec![0.0; order as usize]); //Add a new vector for this timestep's orders
        for o in 0..order {
            let derivative: f32 = ode_system[o as usize](&prior_state, (i as f32)*h); //for each timestep, generate the derivative for each order from all prior values 
            solution[i as usize][o as usize] = prior_state[o as usize] + derivative * h; //generate the next step from the derivative for each order
        }
    }

    solution
}

pub fn rk5(ode: &dyn Fn(f32, f32) -> f32, initial_value: f32, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<f32> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<f32> = vec![initial_value];

    for i in 1..point_number{
        let prior: f32 = solution[(i-1) as usize];

        let k1: f32 = ode(prior, (i as f32)*h);
        let k2: f32 = ode(prior, (i as f32)*h + (h/2.0)) + h*(k1/2.0);
        let k3: f32 = ode(prior, (i as f32)*h + (h/2.0)) + h*(k2/2.0);
        let k4: f32 = ode(prior, (i as f32)*h + h) + h*k3;

        solution.push(prior + ((h/6.0) *(k1 + (2.0*k2) + (2.0*k3) + k4)));
    }

    solution
}

pub fn rk5_system(ode_system: Vec<&dyn Fn(&Vec<f32>, f32) -> f32>, initial_values: Vec<f32>, h:f32, lower_bound: f32, upper_bound: f32) -> Vec<Vec<f32>> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    let mut solution: Vec<Vec<f32>> = vec![initial_values];
    let order: i32 = ode_system.len() as i32;

    for i in 1..point_number{
        let prior_state: Vec<f32> = solution[(i-1) as usize].clone(); //prior estimates for all orders; initial state if this is i=1
        solution.push(vec![0.0; order as usize]); //Add a new vector for this timestep's orders
        for o in 0..order {
            let k1: f32 = ode_system[o as usize](&prior_state, (i as f32)*h);
            let k2: f32 = ode_system[o as usize](&prior_state, (i as f32)*h + (h/2.0)) + h*(k1/2.0);
            let k3: f32 = ode_system[o as usize](&prior_state, (i as f32)*h + (h/2.0)) + h*(k2/2.0);
            let k4: f32 = ode_system[o as usize](&prior_state, (i as f32)*h + h) + h*k3;
            
            solution[i as usize][o as usize] = prior_state[o as usize] + ((h/6.0) *(k1 + (2.0*k2) + (2.0*k3) + k4));
        }
    }

    solution
}