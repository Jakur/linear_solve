extern crate num_traits;
extern crate num_rational;
extern crate nalgebra;

use num_traits::identities::Zero;
use num_traits::identities::One;
use num_rational::Rational64;
use num_rational::Ratio;
use nalgebra::DMatrix;


struct Tableau {
    matrix: DMatrix<Rational64>,
    real_variables: Vec<usize>, //Not slack variables
}

impl Tableau {
    ///Creates a tableaux from a vector of inequalities in the standard form
    fn new_standard(vector: Vec<Vec<Rational64>>) -> Tableau {
        let variables = vector[1].len() - 1;
        let rows = vector.len(); //constraints + 1
        let cols = vector.len() + variables; //constraints + 1 + vars
        let mut matrix = DMatrix::from_element(rows, cols, Ratio::zero());
        //Init variable columns
        for i in 0..rows {
            for j in 0..variables {
                matrix[(i, j)] = vector[i][j];
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
            matrix[(r, cols-1)] = vector[r][variables];
        }
        let real_variables: Vec<_> = (0..variables).collect();
        Tableau {
            matrix,
            real_variables,
        }
    }
    fn new(n: usize, m: usize, slice: &[i64]) -> Tableau {
        //Debug testing array
        let table = Tableau {
            matrix: DMatrix::from_iterator(n, m, slice.iter()
                .map(|n| Ratio::from_integer(*n))),
            real_variables: vec![0, 1],
        };
        table
    }
    fn pivot(&mut self) -> bool {
        //Check optimal
        if self.is_optimal() {
            return false;
        }
        let pivot_col = self.choose_var();
        let mut pivot_row = None;
        let mut candidate_ratio: Rational64 = Ratio::zero(); //Placeholder until set
        //Find the best ratio of nonbasic columns
        for col in 0..self.matrix.ncols() {
            for (index, val) in self.matrix.column(col).iter().enumerate() {
                let pivot_value = self.matrix[(index, pivot_col)];
                if pivot_value > Ratio::zero() {
                    let ratio = self.matrix[(index, self.matrix.ncols()-1)] / pivot_value;
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
        }

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
                        let mut row_copy = DMatrix::from_row_slice(1, self.matrix.ncols(), &row_vec[..]);
                        let mut row = self.matrix.row_mut(row_index);
//                        row.copy_from(&(row + self.matrix.row(r)));
                        let scaling = Ratio::from_integer(-1) * col_val;
//                        row.add_to(scaling, self.matrix.row(r), 1);
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
    ///Choose the nonbasic variable that will have the best effect for optimization
    fn choose_var(&self) -> usize {
        let mut best_col = 0;
        for col in 0..self.matrix.ncols()-1 {
            if self.matrix[(0, col)] < self.matrix[(0, best_col)] {
                best_col = col;
            }
        }
        return best_col;
    }

    fn is_optimal(&self) -> bool {
        for col in 0..self.matrix.ncols()-1 {
            if self.matrix[(0, col)] < Ratio::zero() {
                return false;
            }
        }
        return true;
    }
    ///Find and return the optimal values of the true variables
    fn read_solution(&self) -> Vec<Rational64> {
        let mut values = Vec::new();
        for var_column in self.real_variables.iter() {
            let mut row = None; //Row suspected to hold the value of the variable
            for (index, val) in self.matrix.column(*var_column).iter().enumerate() {
                if val != &Ratio::zero() {
                    match row {
                        Some(_) => {
                            values.push(Ratio::zero()); //Confirmed nonbasic
                            continue;
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

fn main() {
    let n = 5;
    let m = 7;
    let test1 =
        [-4, 2, -3, 0, 2,
        -3, 3, 2, 2, 1,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1,
        0, 6, 3, 5, 4]; //Column major
    let test2 =
    [-1, 2, 4, 2,
    -2, 1, 2, 5,
    1, 1, 3, 5,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
    0, 14, 28, 30];
    let test3 =
    [-1000, 10, 2, 1, 0,
    -1200, 5, 3, 0, 1,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 200, 60, 34, 14];

    //test_string is equivalent to test3
    let test_string = "1000 1200 0\n10 5 200\n2 3 60\n1 0 34\n0 1 14";
    let mut test_vec: Vec<Vec<Rational64>> = test_string.lines()
        .map(|line| line.split_whitespace()
        .map(|x| x.parse().unwrap()).collect()).collect();
    //Make first row, i.e. objective function values, negative
    test_vec[0] = test_vec[0].iter().map(|num| num * -1).collect();

    //let mut table = Tableau::new(n, m, &test3);
    let mut table = Tableau::new_standard(test_vec);

    println!("Starting Tableau: {}", table.matrix);

    while table.pivot() {}

    let solution = table.read_solution();
    println!("Finished Tableau: {}", table.matrix);
    print!("Solution: ");
    for val in solution {
        print!("{} ", val);
    }
    println!("\nObjective Function Value: {}", table.matrix[(0, table.matrix.ncols() - 1)]);
}

