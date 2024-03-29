pub fn ctcs2_1d(pde: &dyn Fn(f32, &Vec<f32>) -> f32, init_con: &dyn Fn(f32) -> f32, constants: Vec<f32>, timespan: f32, time_step: f32, rx: f32, lx: f32, space_step: f32) -> Vec<Vec<f32>> { //pde: (d^2u/dx^2, const_args) -> d^2u/dt^2
    let x_width: i32 = ((rx - lx + (2.0 * space_step)) / space_step) as i32;
    let time_width: i32 = (timespan / time_step) as i32;

    let mut solution: Vec<Vec<f32>> = vec![vec![0.0; x_width as usize]; time_width as usize];

    for i in 0..time_width-1 {
        for x in 1..x_width-1 {

            let iu = i as usize;
            let xu = x as usize;

            if i == 0 {
                solution[iu][xu] = init_con(x as f32 * space_step); //set up initial condition
            }

            else {

                if i==1 {
                    solution[iu][xu] = init_con(x as f32 * space_step); //since we don't actually set the first position
                }

                let d2udx2 = (solution[iu][xu+1] - (2.0 * solution[iu][xu]) + solution[iu][xu-1]) / space_step.powi(2);
                solution[iu+1][xu] = (pde(d2udx2, &constants) * time_step.powi(2))  + (2.0 * solution[iu][xu]) - solution[iu-1][xu];
                //Take the rest of the second order central difference and move it over to the right hand, so *dt^2 and +2f(x)-f(x-h) leaving f(x+h)
            }
        }
    }

    solution.remove(0); //trim off the extraneous state used for the first step

    return solution
}