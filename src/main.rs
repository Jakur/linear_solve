extern crate num;
extern crate num_traits;
extern crate num_rational;
extern crate nalgebra;
extern crate rand;

use std::fmt;

use num_traits::identities::Zero;
use num_traits::identities::One;
use num_rational::Rational64;
use num_rational::Ratio;
use nalgebra::DMatrix;
use rand::Rng;

mod tableau;
mod parse;
use parse::{create_table, FileFormat};
use tableau::Tableau;

struct ParameterLP {
    matrix: DMatrix<Rational64>,
    lambda_count: usize, //Index of the first slack column, also number of lambda
    fixed_lambda: Vec<Rational64>, //Fixed coefficients of lambda for a single tableaux solution
    optim_function: Vec<Rational64>,
    initial_condition: DMatrix<Rational64>,
    artificial: Vec<bool>,
    x_count: usize,
    rng: rand::ThreadRng,
}

impl Tableau for ParameterLP {
    fn matrix(&self) -> &DMatrix<Rational64>{
        &self.matrix
    }

    fn matrix_mut(&mut self) -> &mut DMatrix<Rational64> {
        &mut self.matrix
    }
    fn num_variables(&self) -> usize {
        self.lambda_count
    }
    fn solve(&mut self, phase_one: bool) -> bool {
        let num_solutions = 100; //Todo compute this, or find other halting condition
        let mut solutions = Vec::new();
        let mut counter = 0;
        while solutions.len() < num_solutions && counter < 10000 {
            counter += 1;
            if phase_one {
                let artificial_cols = self.artificial.iter().enumerate()
                    .filter_map(|(index, b)| {
                        if *b {
                            Some(index)
                        } else {
                            None
                        }
                    }).collect();
                if !self.artificial_solve(artificial_cols) {
                    println!("Unable to find a feasible solution under these X.");
                    self.matrix = self.initial_condition.clone();
                    self.update_objective_row();
                    continue;
                }
            }
            while self.pivot() {}
            let lambda = self.read_solution();
            let mut vec = vec![Ratio::zero(); self.x_count + 1];
            for (sl_index, slice) in self.optim_function
                .chunks(self.x_count + 1).enumerate() {
                for (var_index, val) in slice.iter().enumerate() {
                    vec[var_index] += *val * lambda[sl_index];
                }
            }
            let solution = Inequality::new(vec);

            if !solutions.contains(&solution) {
                println!("Solution {}, Attempt {}", solutions.len() + 1, counter);
                print!("Fixed λ coefficients: [");
                for lambda in self.fixed_lambda.iter() {
                    print!("{} ", lambda);
                }
                println!("]");

                print!("Solution vector (in X): ");
                for x in solution.data.iter() {
                    print!("{} ", x);
                }
                print!("\nSolution vector (in Lambda): ");
                for y in lambda.iter() {
                    print!("{} ", y);
                }
                println!("\nResulting Inequality: {}\n", solution.integer_form());

                solutions.push(solution);

            }
            self.matrix = self.initial_condition.clone();
            self.update_objective_row();
        }
        println!("Counter: {}", counter);
        return true
    }
    fn artificial_cols(&self) -> &Vec<bool> {
        &self.artificial
    }
}

impl ParameterLP {
    fn new_standard(x_count: usize, optim: Vec<Rational64>, constraints: Vec<Vec<Rational64>>,
                    equal_rows: Vec<usize>) -> (ParameterLP, bool) {
        //Fix all x initially to zero
        let variables = constraints[0].len() - 1;
        let fixed_lambda = vec![Ratio::from_integer(0); variables]; //Todo fix this
//        let mut fixed_x = vec![Ratio::from((39, 12)); 6];
//        fixed_x[4] = Ratio::from((9, 2));
//        fixed_x[5] = Ratio::from((9, 2));
        let mut counter = x_count;
        let mut init_optim = Vec::new();
        for _ in 0..variables {
            init_optim.push(optim[counter]);
            counter += x_count + 1;
        }
        let (matrix, artificial, phase_one) = matrix_init(&init_optim,
                                                          &constraints, &equal_rows);
        let mut para = ParameterLP {
            matrix: matrix.clone(),
            lambda_count: variables,
            fixed_lambda,
            optim_function: optim,
            initial_condition: matrix,
            artificial,
            x_count,
            rng: rand::thread_rng(),
        };
        para.update_objective_row();
        (para, phase_one)
    }
    fn update_objective_row(&mut self) {
        let mut col = 0;
        while col < self.lambda_count { //For each x column
            let a: i64 = self.rng.gen_range(-100, 100);
            let b: i64 = self.rng.gen_range(1, 2);
            self.matrix_mut()[(0, col)] = Ratio::from((a, b));
            self.fixed_lambda[col] = Ratio::from((a, b));
            col += 1;
        }
    }
}

struct LP {
    matrix: DMatrix<Rational64>,
    artificial: Vec<bool>,
    variables: usize, //Number of non-slack variables
}

impl Tableau for LP {
    fn matrix(&self) -> &DMatrix<Rational64> {
        &self.matrix
    }

    fn matrix_mut(&mut self) -> &mut DMatrix<Rational64> {
        &mut self.matrix
    }

    fn num_variables(&self) -> usize {
        self.variables
    }

    fn solve(&mut self, phase_one: bool) -> bool {
        if phase_one {
            println!("Initiating Phase I to reach a feasible starting point");
            let artificial_cols = self.artificial.iter().enumerate()
                .filter_map(|(index, b)| {
                if *b {
                    Some(index)
                } else {
                    None
                }
            }).collect();
            if self.artificial_solve(artificial_cols) {
                println!("Phase I successful beginning Simplex method proper");
            } else {
                println!("Unable to find a feasible solution.");
                return false
            }
        }
//        println!("Starting Tableau: {}", self.matrix());

        while self.pivot() {}

//        println!("Finished Tableau: {}", self.matrix());
        let solution = self.read_solution();
        print!("Solution: ");
        for val in solution {
            print!("{} ", val);
        }
        println!("\nObjective Function Value: {}", self.matrix()[(0, self.matrix().ncols() - 1)]);
        return true
    }
    fn artificial_cols(&self) -> &Vec<bool> {
        &self.artificial
    }
}


impl LP {
    ///Creates a tableaux from a vector of inequalities in the standard form
    fn new_standard(optim: Vec<Rational64>, constraints: Vec<Vec<Rational64>>,
                    equals: Vec<usize>) -> (LP, bool) {
        let variables = constraints[0].len() - 1;
        let (matrix, art, phase_flag) = matrix_init(&optim, &constraints, &equals);
        let lp = LP {
            matrix,
            artificial: art,
            variables,
        };
        return (lp, phase_flag);
    }
}

///Does initializations for the general simplex method. Returns the initial matrix, a vector
///representing which columns are artificial, and a boolean representing whether or not a Phase 1
///is needed.
fn matrix_init(optim: &Vec<Rational64>, constraints: &Vec<Vec<Rational64>>,
               equals: &Vec<usize>) -> (DMatrix<Rational64>, Vec<bool>, bool) {
    let variables = constraints[0].len() - 1;
    let rows = constraints.len() + 1;
    let cols = constraints.len() + 2 + variables;
    let mut art = vec![false; cols];
    let mut phase_flag = {
        if equals.len() == 0 {
            false
        } else {
            true
        }
    };
    let mut art_cols = Vec::new();
    for r in equals.iter() {
        let col = variables + *r;
        art[col] = true;
        art_cols.push(col);
    }
    //Penultimate column special a_0, for correcting negative inequalities
    let mut matrix = DMatrix::from_element(rows+1, cols, Ratio::zero());
    //Init optimization row
    for i in 0..optim.len() {
        matrix[(0, i)] = optim[i];
    }
    //Init variable columns
    for i in 1..rows {
        for j in 0..variables {
            matrix[(i, j)] = constraints[i-1][j];
        }
    }
    let mut row = 1;
    //Init slack variables
    for col in variables..cols-2 {
        matrix[(row, col)] = Ratio::one();
        row += 1;
        if art[col] {
            matrix[(rows, col)] = Ratio::one();
        }
    }
    //Init row values
    let mut slack_col = variables;
    for r in 1..rows {
        matrix[(r, cols-1)] = constraints[r-1][variables];
        if !art[slack_col] { //If not an equality constraint
            if matrix[(r, cols-1)] < Ratio::zero() { //If negative rhs
                matrix[(r, cols-2)] = Ratio::from_integer(-1); //init a_0 value
                matrix[(rows, cols-2)] = Ratio::one();
                phase_flag = true;
            }
        }
        slack_col += 1;
    }
    //Init optim row
    for i in 0..optim.len() {
        matrix[(0, i)] = optim[i];
    }
    return (matrix, art, phase_flag);
}

#[derive(PartialEq)]
struct Inequality {
    data: Vec<Rational64>,
    less_or_eq: bool,
}

impl Inequality {
    ///Creates a formal inequality type from an vector of the form x_1 + x_2 + ... + const ≤ 0
    fn new(mut data: Vec<Rational64>) -> Inequality {
        //Normalize to set the leftmost nonzero x coefficient to 1
        let mut less_or_eq = true;
        let index = {
            let mut answer = None;
            for i in 0..data.len() {
                if data[i] != Ratio::zero() {
                    answer = Some(i);
                    break;
                }
            }
            answer
        };
        if let Some(num) = index {
            let norm = data[num];
            if norm < Ratio::zero() { //Dividing by a negative number, so flip the direction
                less_or_eq = false;
            }
            for i in num..data.len() { //We know all previous must be 0
                data[i] /= norm;
            }
        }
        Inequality {
            data,
            less_or_eq,
        }
    }
    ///Returns a new equivalent inequality in integer form
    fn integer_form(&self) -> Inequality {
        let mut lcm = num::integer::lcm(*self.data[0].denom(), *self.data[1].denom());
        for i in 2..self.data.len() {
            lcm = num::integer::lcm(lcm, *self.data[i].denom());
        }
        let mut new_vec = Vec::with_capacity(self.data.len());
        for rat in self.data.iter() {
            new_vec.push(*rat * lcm);
        }
        Inequality {
            data: new_vec,
            less_or_eq: self.less_or_eq,
        }
    }
}

impl fmt::Display for Inequality {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let vec = &self.data;
        let first = vec[0..vec.len() - 1].iter().enumerate()
            .find(|tup| *tup.1 != Ratio::zero()); //First nonzero x value
        let symbol = {
            if self.less_or_eq {
                "≤"
            } else {
                "≥"
            }
        };
        match first {
            Some((index, val)) => {
                write!(f, "{}[x{}]", val, index + 1);
                for i in index + 1..vec.len() - 1 {
                    if vec[i] != Ratio::zero() {
                        write!(f, " + {}[x{}]", vec[i], i + 1);
                    }
                }
                write!(f, " {} {}", symbol, vec[vec.len() - 1] * -1)
            }
            None => {
                write!(f, "Trivial solution 0 {} {}", symbol, vec[vec.len() - 1] * -1)
            }
        }
    }
}

fn read_file(file_string: Option<String>) -> Option<(String, FileFormat)> {
    use std::fs::File;
    use std::io::Read;
    match file_string {
        Some(file_string) => {
            let format = {
                if file_string.ends_with("mps") {
                    FileFormat::MPS
                } else {
                    FileFormat::Custom
                }
            };
            let mut f = File::open(file_string).expect("File not found!");
            let mut contents = String::new();
            f.read_to_string(&mut contents).expect("File reading failed!");
            return Some((contents, format))
        },
        None => {
            println!("No input file given.");
            return None;
        }
    };
}

fn main() {
    use std::env;

    let args: Vec<String> = env::args().collect();
    let (input_data, format) = read_file(args.get(1)
        .map(|s| s.to_string())).expect("Unable to read file!");

    let (mut table, phase_one) = create_table(&input_data, format);
    table.solve(phase_one);
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_solutions() {
        let test_arr = [
            "4 3 0\n2 3 6\n-3 2 3\n0 2 5\n2 1 4",
            "20 10 15 0\n3 2 5 55\n2 1 1 26\n1 1 3 30\n5 2 4 57",
            "1000 1200 0\n10 5 200\n2 3 60\n1 0 34\n0 1 14",
        ];
        let solutions = [
            vec![Ratio::new(3, 2), Ratio::from_integer(1)],
            vec![Ratio::new(18, 10), Ratio::new(208, 10), Ratio::new(16, 10)],
            vec![Ratio::from_integer(15), Ratio::from_integer(10)]
        ];
        for i in 0..test_arr.len() {
            let string = test_arr[i];
            let (mut table, phase_one) = create_table(string, FileFormat::Custom);
            table.solve(phase_one);
            assert_eq!(table.read_solution(), solutions[i]);
        }

    }
    #[test]
    fn test_afiro() {
        let (input_data, format) = read_file(Some("afiro.mps".to_string()))
            .expect("Unable to read file!");
        assert_eq!(format, FileFormat::MPS);
        let (mut table, phase_one) = create_table(&input_data, format);
        table.solve(phase_one);
        let objective_value = table.matrix()[(0, table.matrix().ncols() - 1)];
        let objective_float = *objective_value.numer() as f64 / *objective_value.denom() as f64;
        let answer = 464.753142857142857;
        let loss = (objective_float - answer).abs();
        assert!(loss < 1E-7);
    }
}

