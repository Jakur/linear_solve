extern crate num_traits;
extern crate num_rational;
extern crate nalgebra;

use num_traits::identities::Zero;
use num_traits::identities::One;
use num_rational::Rational64;
use num_rational::Ratio;
use nalgebra::DMatrix;

enum TabularData {
    LP(Vec<Rational64>, Vec<Vec<Rational64>>),
    ParameterLP(usize, Vec<Rational64>, Vec<Vec<Rational64>>)
}

trait Tableau {
    fn matrix(&self) -> &DMatrix<Rational64>;
    fn matrix_mut(&mut self) -> &mut DMatrix<Rational64>;
    fn pivot(&mut self) -> bool;
    ///Choose the nonbasic variable that will have the best effect for optimization
    fn choose_var(&self) -> usize {
        let mut best_col = 0;
        let matrix = self.matrix();
        for col in 0..matrix.ncols()-1 {
            if matrix[(0, col)] < matrix[(0, best_col)] {
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
        for index in 1..matrix.nrows() {
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
        for col in 0..matrix.ncols()-1 {
            if matrix[(0, col)] < Ratio::zero() {
                return false;
            }
        }
        return true;
    }
    fn read_solution(&self) -> Vec<Rational64>;
}

struct ParameterLP {
    matrix: DMatrix<Rational64>,
    lambda_count: usize, //Index of the first slack column, also number of lambda
    fixed_x: Vec<Rational64>, //Fixed values of x for a single tableaux solution
    optim_function: Vec<Rational64>,
}

impl Tableau for ParameterLP {
    fn matrix(&self) -> &DMatrix<Rational64>{
        &self.matrix
    }

    fn matrix_mut(&mut self) -> &mut DMatrix<Rational64> {
        &mut self.matrix
    }

    fn pivot(&mut self) -> bool {
        //Check optimal
        if self.is_optimal() {
            self.print_objective_fn();
            return false;
        }
        let pivot_col = self.choose_var();
        let pivot_row = self.choose_row(pivot_col);

        match pivot_row {
            Some(r) => {
                //Set the value of (pivot_row, pivot_col) to 1
                let value = Ratio::one() / self.matrix[(r, pivot_col)];
                {
                    let mut row = self.matrix.row_mut(r);
                    row *= value;
                }
                //Zero out the rest of the pivot_column with row operations using pivot_row
                let mut row_vec = Vec::new();
                for value in self.matrix.row(r).iter() {
                    row_vec.push(value.clone());
                }
                //Todo operations on variable objective function?
                for row_index in 0..self.matrix.nrows() {
                    let col_val = self.matrix[(row_index, pivot_col)];
                    if row_index != r && col_val != Ratio::zero() {
                        let mut row_copy = DMatrix::from_row_slice(1, self.matrix.ncols(),
                                                                   &row_vec[..]);
                        let mut row = self.matrix.row_mut(row_index);
                        let scaling = Ratio::from_integer(-1) * col_val;
                        row_copy *= scaling;
                        row += row_copy;
                    }
                }

                return true
            }
            None => {
                return false
            }
        }
    }

    fn read_solution(&self) -> Vec<Ratio<i64>> {
        unimplemented!()
    }
}

impl ParameterLP {
    fn new_standard(x_count: usize, optim: Vec<Rational64>,
                    constraints: Vec<Vec<Rational64>>) -> ParameterLP {
        let variables = constraints[0].len() - 1;
        let rows = constraints.len() + 1;
        let cols = constraints.len() + 1 + variables;
        let mut matrix = DMatrix::from_element(rows, cols, Ratio::zero());
        //Init variable columns
        for i in 1..rows {
            for j in 0..variables {
                matrix[(i, j)] = constraints[i-1][j];
            }
        }
        let mut row = 1;
        //Init slack variables
        for col in variables..cols-1 {
            matrix[(row, col)] = Ratio::one();
            row += 1;
        }
        //Init row values
        for r in 1..rows {
            matrix[(r, cols-1)] = constraints[r-1][variables];
        }
        //Fix all x initially to zero
        let fixed_x = vec![Ratio::from_integer(3); x_count]; //Todo fix this

        let mut counter = x_count;
        println!("{:?}", optim);
        for col in 0..variables {
            matrix[(0, col)] = optim[counter];
            println!("Optim of col: {}", optim[counter]);
            counter += x_count + 1;
        }
        let mut para = ParameterLP {
            matrix,
            lambda_count: variables,
            fixed_x,
            optim_function: optim,
        };
        para.update_objective_row();
        para
    }
    fn update_objective_row(&mut self) {
        let mut col = 0;
        while col < self.lambda_count { //For each x column
            let mut total = Ratio::zero();
            let offset = col * (self.fixed_x.len() + 1);
            for index in 0..self.fixed_x.len() {
                total += self.optim_function[offset+index] * self.fixed_x[index];
            }
            total += self.optim_function[offset + self.fixed_x.len()]; //Constant term
            self.matrix_mut()[(0, col)] = total;
            col += 1;
        }
    }
    fn print_objective_fn(&self) {
        let mut counter = 1;
        for chunk in self.optim_function.chunks(self.fixed_x.len() + 1) {
            println!("Lambda{} {:?}", counter, chunk);
            counter += 1;
        }
    }
}

struct LP {
    matrix: DMatrix<Rational64>,
    real_variables: Vec<usize>, //Not slack variables
}

impl Tableau for LP {
    fn matrix(&self) -> &DMatrix<Rational64> {
        &self.matrix
    }

    fn matrix_mut(&mut self) -> &mut DMatrix<Rational64> {
        &mut self.matrix
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
                //Set the value of (pivot_row, pivot_col) to 1
                let value = Ratio::one() / self.matrix[(r, pivot_col)];
                {
                    let mut row = self.matrix.row_mut(r);
                    row *= value;
                }
                //Zero out the rest of the pivot_column with row operations using pivot_row
                let mut row_vec = Vec::new();
                for value in self.matrix.row(r).iter() {
                    row_vec.push(value.clone());
                }
                for row_index in 0..self.matrix.nrows() {
                    let col_val = self.matrix[(row_index, pivot_col)];
                    if row_index != r && col_val != Ratio::zero() {
                        let mut row_copy = DMatrix::from_row_slice(1, self.matrix.ncols(),
                                                                   &row_vec[..]);
                        let mut row = self.matrix.row_mut(row_index);
                        let scaling = Ratio::from_integer(-1) * col_val;
                        row_copy *= scaling;
                        row += row_copy;
                    }
                }

                return true
            }
            None => {
                return false
            }
        }
    }
    ///Find and return the optimal values of the true variables
    fn read_solution(&self) -> Vec<Rational64> {
        let mut values = Vec::new();
        'outer: for var_column in self.real_variables.iter() {
            let mut row = None; //Row suspected to hold the value of the variable
            for (index, val) in self.matrix.column(*var_column).iter().enumerate() {
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
                values.push(self.matrix[(r, self.matrix.ncols()-1)]);
            } else { //The whole column was zero. Unknown if this can happen
                values.push(Ratio::zero());
            }
        }
        return values;
    }
}


impl LP {
    ///Creates a tableaux from a vector of inequalities in the standard form
    fn new_standard(optim: Vec<Rational64>, constraints: Vec<Vec<Rational64>>) -> LP {
        let variables = constraints[0].len() - 1;
        let rows = constraints.len() + 1;
        let cols = constraints.len() + 1 + variables;
        let mut matrix = DMatrix::from_element(rows, cols, Ratio::zero());
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
        for col in variables..cols-1 {
            matrix[(row, col)] = Ratio::one();
            row += 1;
        }
        //Init row values
        for r in 1..rows {
            matrix[(r, cols-1)] = constraints[r-1][variables];
        }
        let real_variables: Vec<_> = (0..variables).collect();
        LP {
            matrix,
            real_variables,
        }
    }
    ///Creates a Tableau from a single slice representing a standard LP
    fn new_from_tabular(n: usize, m: usize, slice: &[i64]) -> LP {
        //Debug testing array
        let table = LP {
            matrix: DMatrix::from_iterator(n, m, slice.iter()
                .map(|n| Ratio::from_integer(*n))),
            real_variables: vec![0, 1],
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
    //Add the -1 * row to the tableaux to represent row >= x and row <= x
    for index in equals.into_iter() {
        let opposite: Vec<_> = vec[index].iter()
            .map(|num| Ratio::from_integer(-1) * num).collect();
        vec.push(opposite);
    }
    if split.len() == 1 { //Normal LP
        let opt = split[0].split_whitespace().filter_map(|x| {
            //Make optimization values negative
            match x.parse::<Rational64>() {
                Ok(num) => {Some(num * -1)},
                _ => None,
            }
        }).collect();
        return TabularData::LP(opt, vec);
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
        for x in xs.iter() {
            println!("{}", x);
        }
        println!("X count: {}", x_count);
        return TabularData::ParameterLP(x_count, xs, vec);
    }
}

fn create_table(text: &str) -> Box<Tableau> {
    let parsed = parse_inequalities(text);
    match parsed {
        TabularData::LP(opt, con) => {
            Box::new(LP::new_standard(opt, con))
        }
        TabularData::ParameterLP(x_count, var_opt, lambda_con) => {
            Box::new(ParameterLP::new_standard(x_count, var_opt, lambda_con))
        }
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

    let mut table = create_table(&input_data);

    println!("Starting Tableau: {}", table.matrix());

    while table.pivot() {}

    println!("Finished Tableau: {}", table.matrix());
//    let solution = table.read_solution();
//    print!("Solution: ");
//    for val in solution {
//        print!("{} ", val);
//    }
    println!("\nObjective Function Value: {}", table.matrix()[(0, table.matrix().ncols() - 1)]);
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_standard_parsing() {
        let test_str1 = "4 3 0\n2 3 6\n-3 2 3\n0 2 5\n2 1 4";
        let test_str2 = "20 10 15 0\n3 2 5 55\n2 1 1 26\n1 1 3 30\n5 2 4 57";
        let test_str3 = "1000 1200 0\n10 5 200\n2 3 60\n1 0 34\n0 1 14";
        let tab1 = create_table(test_str1);
        let tab2 = create_table(test_str2);
        let tab3 = create_table(test_str3);
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
            let mut table = create_table(string);
            while table.pivot() {}
            assert_eq!(table.read_solution(), solutions[i]);
        }

    }
}

