extern crate num_traits;
extern crate num_rational;
extern crate nalgebra;

use num_traits::identities::Zero;
use num_traits::identities::One;
use num_rational::Rational64;
use num_rational::Ratio;
use nalgebra::DMatrix;

mod tableau;
use tableau::Tableau;

enum TabularData {
    LP(Vec<Rational64>, Vec<Vec<Rational64>>, Vec<usize>), //Opt, Constraints, Equality Indices
    //x_count, Opt, Constr, Eq Indices
    ParameterLP(usize, Vec<Rational64>, Vec<Vec<Rational64>>, Vec<usize>),
}

struct ParameterLP {
    matrix: DMatrix<Rational64>,
    lambda_count: usize, //Index of the first slack column, also number of lambda
    fixed_x: Vec<Rational64>, //Fixed values of x for a single tableaux solution
    optim_function: PolyVec,
    initial_condition: DMatrix<Rational64>,
    artificial: Vec<bool>,
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
        let num_solutions = 2; //Todo compute this, or find other halting condition
        let mut solutions = Vec::new();
        while solutions.len() < num_solutions {
            println!("X: {:?}", self.fixed_x);
            let mut print = Vec::new();
            if phase_one {
                print.push(format!("Initial Matrix: {}", self.matrix));
                print.push("Initiating Phase I to reach a feasible starting point".to_string());
                let artificial_cols = self.artificial.iter().enumerate()
                    .filter_map(|(index, b)| {
                        if *b {
                            Some(index)
                        } else {
                            None
                        }
                    }).collect();
                if self.artificial_solve(artificial_cols) {
                    print.push(format!("Phase I successful beginning Simplex method proper with:\n{}",
                             self.matrix));
                } else {
                    println!("Unable to find a feasible solution under these X.");
                    self.matrix = self.initial_condition.clone();
                    self.update_fixed_x();
                    self.update_objective_row();
                    continue;
                }
            }
            while self.pivot() {}
            let lambda = self.read_solution();
            let mut vec = vec![Ratio::zero(); self.fixed_x.len() + 1];
            for (sl_index, slice) in self.optim_function.data
                .chunks(self.fixed_x.len() + 1).enumerate() {
                for (var_index, val) in slice.iter().enumerate() {
                    vec[var_index] += *val * lambda[sl_index];
                }
            }
            if !solutions.contains(&vec) {
                for line in print {
                    println!("{}", line);
                }
                print!("Solution vector (in X): ");
                for x in vec.iter() {
                    print!("{} ", x);
                }
                print!("\nSolution vector (in Lambda): ");
                for y in lambda.iter() {
                    print!("{} ", y);
                }
                let mut art = Vec::new();
                for (index, b) in self.artificial_cols().iter().enumerate() {
                    if *b {
                        art.push(index);
                    }
                }
                println!("\nArtificial columns: {:?}", art);
                println!("{}", self.matrix);
                solutions.push(vec);
            }
            self.matrix = self.initial_condition.clone();
            self.update_fixed_x();
            self.update_objective_row();
        }
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
        let fixed_x = vec![Ratio::from_integer(0); x_count]; //Todo fix this
        let mut counter = x_count;
        let mut init_optim = Vec::new();
        for _ in 0..variables {
            init_optim.push(optim[counter]);
            counter += x_count + 1;
        }
        let (matrix, artificial, phase_one) = matrix_init(&init_optim,
                                                          &constraints, &equal_rows);

        let opt = PolyVec::new(optim, fixed_x.len() + 1, variables);
        let mut para = ParameterLP {
            matrix: matrix.clone(),
            lambda_count: variables,
            fixed_x,
            optim_function: opt,
            initial_condition: matrix,
            artificial,
        };
        para.update_objective_row(); //Necessary if all fixed_x != 0
        (para, phase_one)
    }
    fn update_objective_row(&mut self) {
        let mut col = 0;
        while col < self.lambda_count { //For each x column
            let mut total = Ratio::zero();
            let offset = col * (self.fixed_x.len() + 1);
            for index in 0..self.fixed_x.len() {
                total += self.optim_function.data[offset+index] * self.fixed_x[index];
            }
            total += self.optim_function.data[offset + self.fixed_x.len()]; //Constant term
            self.matrix_mut()[(0, col)] = total;
            col += 1;
        }
    }
    fn update_fixed_x(&mut self) {
        //Todo find heuristic
        for i in 0..self.fixed_x.len() {
            self.fixed_x[i] += 1;
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
        println!("Starting Tableau: {}", self.matrix());

        while self.pivot() {}

        println!("Finished Tableau: {}", self.matrix());
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
    ///Creates a Tableau from a single slice representing a standard LP
    fn new_from_tabular(n: usize, m: usize, slice: &[i64]) -> LP {
        //Debug testing array
        let table = LP {
            matrix: DMatrix::from_iterator(n, m, slice.iter()
                .map(|n| Ratio::from_integer(*n))),
            artificial: Vec::new(),
            variables: 2,
        };
        table
    }
}
///Parses a matrix as written in text as a vector of row vectors
fn parse_inequalities(text: &str) -> TabularData {
    let mut greater = Vec::new();
    let mut equals = Vec::new();
    let mut iter = text.lines().enumerate();
    let (_, opt_line) = iter.next().expect("Optimization line malformed");
    let split: Vec<_> = opt_line.split(|c| c == '[' || c == ']').collect();
    let mut vec: Vec<Vec<Rational64>> = iter
        .map(|(index, line)| line.split_whitespace()
        .filter_map(|x| {
            match x.parse() {
                Ok(num) => {
                    Some(num)
                }
                _ => {
                    match x {
                        ">=" => {greater.push(index-1)},
                        "=" => {equals.push(index-1)},
                        "<=" => {} //Already canonical form, okay
                        _ => {panic!("Unrecognized string supplied.")}
                    }
                    None
                }
            }
        }).collect()).collect();

    //Flip the signs to make greater than or equal to into less than or equal to
    for index in greater.into_iter() {
        for i in 0..vec[index].len() {
            vec[index][i] *= -1;
        }
    }
    //Make the rhs of all equality constraints nonnegative
    for index in equals.iter() {
        let index = *index;
        let length = vec[index].len();
        if vec[index][length-1] < Ratio::zero() {
            let opposite: Vec<_> = vec[index].iter()
                .map(|num| Ratio::from_integer(-1) * num).collect();
            vec[index] = opposite;
        }
    }
    if split.len() == 1 { //Normal LP
        let opt = split[0].split_whitespace().filter_map(|x| {
            //Make optimization values negative
            match x.parse::<Rational64>() {
                Ok(num) => {Some(num * -1)},
                _ => None,
            }
        }).collect();
        return TabularData::LP(opt, vec, equals);
    } else { //Parametrized LP
        let x_count = split[1].split_whitespace().count() - 1; //-1 for constant term
        let xs: Vec<_> = split.into_iter().map(|lambda| {
            lambda.split_whitespace().filter_map(|x| {
                match x.parse::<Rational64>() {
                    Ok(num) => {Some(num * -1)}, //Todo figure out if * -1 is appropriate
                    _ => None,
                }
            })
        }).flatten().collect();
        println!("X count: {}", x_count);
        return TabularData::ParameterLP(x_count, xs, vec, equals);
    }
}

fn create_table(text: &str) -> (Box<Tableau>, bool) {
    let parsed = parse_inequalities(text);
    match parsed {
        TabularData::LP(opt, con, eq) => {
            let (lp, phase_one) = LP::new_standard(opt, con, eq);
            (Box::new(lp), phase_one)
        }
        TabularData::ParameterLP(x_count, var_opt, lambda_con, eq) => {
            let (plp, phase_one) = ParameterLP::new_standard(x_count, var_opt, lambda_con, eq);
            (Box::new(plp), phase_one)
        }
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

struct PolyVec {
    data: Vec<Rational64>,
    poly_len: usize, //Number of terms in each polynomial
    poly_count: usize, //Number of polynomials
}

impl PolyVec {
    fn new(data: Vec<Rational64>, poly_len: usize, poly_count: usize) -> PolyVec {
        PolyVec {
            data,
            poly_len,
            poly_count
        }
    }
    fn operate(&mut self, row: &DMatrix<Rational64>, multipliers: Vec<Rational64>) {
        for (col, m) in multipliers.into_iter().enumerate() {
            if m == Ratio::zero() {continue;}
            let r = row * m;
            for val in r.iter() {
               for p_index in 0..self.poly_count {
                   let index = self.index(p_index, col);
                   self.data[index] += *val;
               }
            }
        }
    }
    fn index(&self, poly_index: usize, col: usize) -> usize {
        return poly_index * self.poly_len + col
    }
}

fn main() {
    use std::env;
    use std::fs::File;
    use std::io::Read;

    let args: Vec<String> = env::args().collect();
    let input_data = match args.get(1) {
        Some(file_string) => {
            let mut f = File::open(file_string).expect("File not found!");
            let mut contents = String::new();
            f.read_to_string(&mut contents).expect("File reading failed!");
            contents
        },
        _ => {
            println!("No input file given.");
            return;
        }
    };

    let (mut table, phase_one) = create_table(&input_data);
    table.solve(phase_one);
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_standard_parsing() {
        let test_str1 = "4 3 0\n2 3 6\n-3 2 3\n0 2 5\n2 1 4";
        let test_str2 = "20 10 15 0\n3 2 5 55\n2 1 1 26\n1 1 3 30\n5 2 4 57";
        let test_str3 = "1000 1200 0\n10 5 200\n2 3 60\n1 0 34\n0 1 14";
        let (tab1, _) = create_table(test_str1);
        let (tab2, _) = create_table(test_str2);
        let (tab3, _) = create_table(test_str3);
        let test1 = [
            -4, 2, -3, 0, 2,
            -3, 3, 2, 2, 1,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1,
            0, 6, 3, 5, 4]; //Column major
        let test2 = [
            -20, 3, 2, 1, 5,
            -10, 2, 1, 1, 2,
            -15, 5, 1, 3, 4,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1,
            0, 55, 26, 30, 57];
        let test3 = [
            -1000, 10, 2, 1, 0,
            -1200, 5, 3, 0, 1,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1,
            0, 200, 60, 34, 14];
        let tab1 = (tab1.matrix(),
                    LP::new_from_tabular(5, 7, &test1));
        let tab2 = (tab2.matrix(),
                    LP::new_from_tabular(5, 8, &test2));
        let tab3 = (tab3.matrix(),
                    LP::new_from_tabular(5, 7, &test3));
        assert_eq!(tab1.0, &tab1.1.matrix);
        assert_eq!(tab2.0, &tab2.1.matrix);
        assert_eq!(tab3.0, &tab3.1.matrix);
    }
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
            let (mut table, phase_one) = create_table(string);
            table.solve(phase_one);
            assert_eq!(table.read_solution(), solutions[i]);
        }

    }
}

