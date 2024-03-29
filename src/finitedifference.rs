pub fn df1(eval_function: &dyn Fn(f32) -> f32, x: f32, h: f32, mode: i8) -> f32 {
    return match mode {
        0 => (eval_function(x + h) - eval_function(x)) / h, //forward
        1 => (eval_function(x) - eval_function(x - h)) / h, //backward
        2 | _ => (eval_function(x + (h / 2.0)) - eval_function(x - (h / 2.0))) / h, //central
    }
}

pub fn df2(eval_function: &dyn Fn(f32, f32) -> f32, x: f32, y: f32, h: f32, k: f32, mode: i8) -> [f32; 2] { 

    let partial_x = |xp: f32| {eval_function(xp, y)}; //∂f/dx
    let partial_y = |yp: f32| {eval_function(x, yp)}; //∂f/dy

    let pder_x: f32 = df1(&partial_x, x, h, mode);
    let pder_y: f32 = df1(&partial_y, y, k, mode);

    let derivative: [f32; 2] = [pder_x, pder_y];

    derivative
}

pub fn differentiate_1d(eval_function: &dyn Fn(f32) -> f32, lower_bound:f32, upper_bound:f32, h:f32, mode: i8) -> Vec<f32> {
    let point_number: i32 = ((upper_bound - lower_bound) / h).round() as i32;
    
    let mut derivative: Vec<f32> = vec![0.0; point_number as usize];

    for i in 0..point_number {
        let x: f32 = lower_bound + ((i as f32) * h);
        derivative[i as usize] = df1(eval_function, x, h, mode);
    }

    derivative
}

pub fn differentiate_2d(eval_function: &dyn Fn(f32, f32) -> f32, lower_bound_x:f32, upper_bound_x:f32, lower_bound_y:f32, upper_bound_y:f32, h: f32, k: f32, mode: i8) -> Vec<Vec<[f32; 2]>> {

    let point_number_x: i32 = ((upper_bound_x - lower_bound_x) / h).round() as i32;
    let point_number_y: i32 = ((upper_bound_y - lower_bound_y) / k).round() as i32;


    let mut gradient: Vec<Vec<[f32; 2]>> = vec![vec![[0.0; 2]; point_number_y as usize]; point_number_x as usize]; // [x coordinate][y coordinate][(∂f/dx, ∂f/dy)]

    for i in 0..point_number_x{
        for j in 0..point_number_y{
            let x: f32 = lower_bound_x + ((i as f32) * h);
            let y: f32 = lower_bound_y + ((j as f32) * k);
            gradient[i as usize][j as usize] = df2(eval_function, x, y, h, k, mode);
        }
    }

    gradient
}