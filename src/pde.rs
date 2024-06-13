pub fn ctcs22_1d(pde: &dyn Fn(f32, &Vec<f32>) -> f32, init_con: &dyn Fn(f32) -> f32, constants: Vec<f32>, timespan: f32, time_step: f32, rx: f32, lx: f32, space_step: f32) -> Vec<Vec<f32>> { //pde: (d^2u/dx^2, const_args) -> d^2u/dt^2
    let x_width: i32 = ((rx - lx + (2.0 * space_step)) / space_step) as i32;
    let time_width: i32 = (timespan / time_step) as i32;

    let mut solution: Vec<Vec<f32>> = vec![vec![0.0; x_width as usize]; time_width as usize];

    for x in 1..x_width-1 {
        solution[0][x as usize] = init_con(x as f32 * space_step); //set up initial condition
        solution[1][x as usize] = init_con(x as f32 * space_step);
    }

    for i in 1..time_width-1 {
        for x in 1..x_width-1 {

            let iu = i as usize;
            let xu = x as usize;

            let d2udx2 = (solution[iu][xu+1] - (2.0 * solution[iu][xu]) + solution[iu][xu-1]) / space_step.powi(2);
            solution[iu+1][xu] = (pde(d2udx2, &constants) * time_step.powi(2))  + (2.0 * solution[iu][xu]) - solution[iu-1][xu];
        }
    }

    return solution
}

pub fn ctcs22_2d(pde: &dyn Fn(Vec<f32>, &Vec<f32>) -> f32, init_con: &dyn Fn(f32, f32) -> f32, constants: Vec<f32>, timespan: f32, time_step: f32, rx: f32, lx: f32, space_step: f32) -> Vec<Vec<Vec<f32>>> { //pde: ([d^2u/dx^2, d^2u/dy^2], const_args) -> d^2u/dt^2
    let xy_width: i32 = ((rx - lx + (2.0 * space_step)) / space_step) as i32;
    let time_width: i32 = (timespan / time_step) as i32;

    let mut solution: Vec<Vec<Vec<f32>>> = vec![vec![vec![0.0; xy_width as usize]; xy_width as usize]; time_width as usize]; //time[x[y]]

    for x in 1..xy_width-1 {
        for y in 1..xy_width - 1 {
            solution[0][x as usize][y as usize] = init_con(x as f32 * space_step, y as f32 * space_step); //set up initial condition
            solution[1][x as usize][y as usize] = init_con(x as f32 * space_step, y as f32 * space_step); //set up initial condition
        }
    }

    for i in 1..time_width-1 {
        for x in 1..xy_width-1 {
            for y in 1..xy_width-1 {

                let iu = i as usize;
                let xu = x as usize;
                let yu = y as usize;

                let d2udx2 = (solution[iu][xu+1][yu] - (2.0 * solution[iu][xu][yu]) + solution[iu][xu-1][yu]) / space_step.powi(2);
                let d2udy2 = (solution[iu][xu][yu+1] - (2.0 * solution[iu][xu][yu]) + solution[iu][xu][yu-1]) / space_step.powi(2);

                solution[iu+1][xu][yu] = (pde(vec![d2udx2, d2udy2], &constants) * time_step.powi(2))  + (2.0 * solution[iu][xu][yu]) - solution[iu-1][xu][yu];
                //Take the rest of the second order central difference and move it over to the right hand, so *dt^2 and +2f(x)-f(x-h) leaving f(x+h
            }
        }
    }

    return solution
}