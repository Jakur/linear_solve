extern crate num_traits;
extern crate num_rational;
extern crate nalgebra;

use num_traits::identities::Zero;
use num_traits::identities::One;
use num_rational::Rational64;
use num_rational::Ratio;
use nalgebra::DMatrix;

enum TabularData {
    LP(Vec<Rational64>, Vec<Vec<Rational64>>, Vec<usize>), //Opt, Constraints, Equality Indices
    //x_count, Opt, Constr, Eq Indices
    ParameterLP(usize, Vec<Rational64>, Vec<Vec<Rational64>>, Vec<usize>),
}

trait Tableau {
    fn matrix(&self) -> &DMatrix<Rational64>;
    fn matrix_mut(&mut self) -> &mut DMatrix<Rational64>;
    fn num_variables(&self) -> usize;
    fn solve(&mut self, phase_one: bool) -> bool;
    ///Eliminate on this column, making its value 1 and all other values in the column 0
    fn eliminate(&mut self, row_index: usize, col_index: usize) {
        //Set the value of (pivot_row, pivot_col) to 1
        let value = Ratio::one() / self.matrix()[(row_index, col_index)];
        {
            let mut row = self.matrix_mut().row_mut(row_index);
            row *= value;
        }
        let mut row_vec = Vec::new();
        for value in self.matrix().row(row_index).iter() {
            row_vec.push(value.clone());
        }
        let mat = self.matrix_mut();
        for other_row in 0..mat.nrows() {
            let col_val = mat[(other_row, col_index)];
            if other_row != row_index && col_val != Ratio::zero() {
                let mut row_copy = DMatrix::from_row_slice(1, mat.ncols(),
                                                           &row_vec[..]);
                let mut row = mat.row_mut(other_row);
                let scaling = Ratio::from_integer(-1) * col_val;
                row_copy *= scaling;
                row += row_copy;
            }
        }
    }
    fn pivot(&mut self) -> bool {
        //Check optimal
        if self.is_optimal() {
            return false;
        }
        let pivot_col = self.choose_var();
        let pivot_row = self.choose_row(pivot_col);
        match pivot_row {
            Some(r) => {
                //Zero out the rest of the pivot_column with row operations using pivot_row
                self.eliminate(r, pivot_col);
                return true
            }
            None => {
                return false
            }
        }
    }
    ///Choose the nonbasic variable that will have the best effect for optimization
    fn choose_var(&self) -> usize {
        let mut best_col = 0;
        let matrix = self.matrix();
        //println!("Art cols: {:?}", self.artificial_cols());
        for col in 0..matrix.ncols()-2 {
            if matrix[(0, col)] < matrix[(0, best_col)] && !self.artificial_cols()[col] {
                best_col = col;
            }
        }
        return best_col;
    }
    fn choose_row(&self, pivot_col: usize) -> Option<usize> {
        let matrix = self.matrix();
        let mut pivot_row = None;
        let mut candidate_ratio: Rational64 = Ratio::zero(); //Placeholder until set
        //Find the best ratio for the pivot column
        for index in 1..matrix.nrows()-1 {
            let pivot_value = matrix[(index, pivot_col)];
            if pivot_value > Ratio::zero() {
                let ratio = matrix[(index, matrix.ncols()-1)] / pivot_value;
                match pivot_row {
                    Some(_r) => {
                        if ratio < candidate_ratio {
                            pivot_row = Some(index);
                            candidate_ratio = ratio;
                        }
                    }
                    None => {
                        pivot_row = Some(index);
                        candidate_ratio = ratio;
                    }
                }
            }
        }
        return pivot_row;
    }
    fn is_optimal(&self) -> bool {
        let matrix = self.matrix();
        for col in 0..matrix.ncols() - 2 {
            if matrix[(0, col)] < Ratio::zero() && !self.artificial_cols()[col] {
                return false;
            }
        }
        return true;
    }
    ///Try to eliminate artificial variables to get a feasible initial tableau.
    fn artificial_solve(&mut self, artificial_cols: Vec<usize>) -> bool {
        let nrows = self.matrix().nrows();
        let ncols = self.matrix().ncols();
        for col in artificial_cols {
            //println!("{}", self.matrix());
            let row = col - self.num_variables() + 1; //Todo verify this
            //println!("Row: {} Col: {}", row, col);
            self.eliminate(row, col);
        }
        //Eliminate a_0, isolating it on the row with the most negative slack
        let mut best_row = None;
        for row in 1..nrows - 1 {
            if self.matrix()[(row, ncols - 2)] == Ratio::from_integer(-1) {
                match best_row {
                    Some(b) => {
                        if self.matrix()[(row, ncols-1)] < self.matrix()[(b, ncols-1)] {
                            best_row = Some(row);
                        }
                    },
                    None => {best_row = Some(row);}
                }
            }
        }
        //println!("{}", self.matrix());
        match best_row {
            Some(row) => {
                self.eliminate(row, ncols - 2);
            }
            None => {} //No a_0 to eliminate
        }
        //Find optimality for the bottom row--the w row
        while self.matrix()[(nrows-1, ncols-1)] < Ratio::zero() {
            //println!("{}", self.matrix());
            let mut best_col = None;
            for col in 1..ncols - 1 {
                match best_col {
                    Some(b) => {
                        if self.matrix()[(nrows-1, col)] < self.matrix()[(nrows-1, b)] {
                            best_col = Some(col);
                        }
                    }
                    None => {best_col = Some(col);}
                }
            }
            match best_col {
                Some(best_col) => {
                    if self.matrix()[(nrows-1, best_col)] < Ratio::zero() {
                        let row = self.choose_row(best_col);
                        if let Some(best_row) = row {
                            self.eliminate(best_row, best_col);
                        } else {
                            if self.matrix()[(nrows-1, ncols-1)] < Ratio::zero() {
                                return false; //Infeasible
                            } else {
                                return true;
                            }
                        }
                    }
                }
                None => {
                    if self.matrix()[(nrows-1, ncols-1)] < Ratio::zero() {
                        return false; //Infeasible
                    } else {
                        return true;
                    }
                }
            }
        }
        //println!("{}", self.matrix());
        true
    }
    ///Find and return the optimal values of the true variables
    fn read_solution(&self) -> Vec<Rational64> {
        let mut values = Vec::new();
        'outer: for var_column in 0..self.num_variables() {
            let mut row = None; //Row suspected to hold the value of the variable
            for (index, val) in self.matrix().column(var_column).iter().enumerate() {
                if val != &Ratio::zero() {
                    match row {
                        Some(_) => {
                            values.push(Ratio::zero()); //Confirmed nonbasic
                            continue 'outer;
                        }
                        None => {
                            row = Some(index);
                        }
                    }
                }
            }
            if let Some(r) = row {
                values.push(self.matrix()[(r, self.matrix().ncols()-1)]);
            } else { //The whole column was zero. Unknown if this can happen
                values.push(Ratio::zero());
            }
        }
        return values;
    }
    fn artificial_cols(&self) -> &Vec<bool>;
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

