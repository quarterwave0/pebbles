pub fn ctcs22_1d(pde: &dyn Fn(f32, &Vec<f32>) -> f32, init_con: &dyn Fn(f32) -> f32, constants: Vec<f32>, timespan: f32, time_step: f32, rx: f32, lx: f32, space_step: f32) -> Vec<Vec<f32>> { //pde: (d^2u/dx^2, const_args) -> d^2u/dt^2
    let x_width: i32 = ((rx - lx + (2.0 * space_step)) / space_step) as i32;
    let time_width: f32 = timespan / time_step;

    let x_width_u: usize = x_width as usize;
    let time_width_u: usize = time_width as usize;

    let mut solution: Vec<Vec<f32>> = vec![vec![0.0; x_width_u]; time_width_u];

    for xi in 1..x_width-1 {
        let x: f32 = xi as f32;
        let xu: usize = xi as usize;

        solution[0][xu] = init_con(x * space_step); //set up initial condition
        solution[1][xu] = init_con(x * space_step);
    }

    for i in 1..time_width_u-1 {
        for x in 1..x_width_u-1 {

            let d2udx2 = (solution[i][x+1] - (2.0 * solution[i][x]) + solution[i][x-1]) / space_step.powi(2);
            solution[i+1][x] = (pde(d2udx2, &constants) * time_step.powi(2))  + (2.0 * solution[i][x]) - solution[i-1][x];
        }
    }

    solution
}

pub fn ctcs22_2d(pde: &dyn Fn(Vec<f32>, &Vec<f32>) -> f32, init_con: &dyn Fn(f32, f32) -> f32, constants: Vec<f32>, timespan: f32, time_step: f32, rx: f32, lx: f32, space_step: f32) -> Vec<Vec<Vec<f32>>> { //pde: ([d^2u/dx^2, d^2u/dy^2], const_args) -> d^2u/dt^2
    let xy_width: i32 = ((rx - lx + (2.0 * space_step)) / space_step) as i32;
    let time_width: f32 = timespan / time_step;

    let xy_width_u: usize =  xy_width as usize;
    let time_width_u: usize = time_width as usize;

    let mut solution: Vec<Vec<Vec<f32>>> = vec![vec![vec![0.0; xy_width_u]; xy_width_u]; time_width_u]; //time[x[y]]

    for xi in 1..xy_width-1 {
        let x: f32 = xi as f32;
        for yi in 1..xy_width-1 {
            let y: f32 = yi as f32;

            solution[0][x as usize][y as usize] = init_con(x * space_step, y * space_step); //set up initial condition
            solution[1][x as usize][y as usize] = init_con(x * space_step, y * space_step); //set up initial condition
        }
    }

    for i in 1..time_width_u-1 {
        for x in 1..xy_width_u-1 {
            for y in 1..xy_width_u-1 {

                let d2udx2 = (solution[i][x+1][y] - (2.0 * solution[i][x][y]) + solution[i][x-1][y]) / space_step.powi(2);
                let d2udy2 = (solution[i][x][y+1] - (2.0 * solution[i][x][y]) + solution[i][x][y-1]) / space_step.powi(2);

                solution[i+1][x][y] = (pde(vec![d2udx2, d2udy2], &constants) * time_step.powi(2))  + (2.0 * solution[i][x][y]) - solution[i-1][x][y];
                //Take the rest of the second order central difference and move it over to the right hand, so *dt^2 and +2f(x)-f(x-h) leaving f(x+h)
            }
        }
    }

    solution
}