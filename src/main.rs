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
    nonbasic_columns: Vec<usize>,
}

impl Tableau {
    fn new(n: usize, m: usize) -> Tableau {
        //Debug testing array
        let arr = [1, 0, 0, 0, 0,
        -4, 2, -3, 0, 2,
        -3, 3, 2, 2, 1,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1,
        0, 6, 3, 5, 4]; //Column major

        let table = Tableau {
            matrix: DMatrix::from_iterator(n, m, arr.iter()
                .map(|n| Ratio::from_integer(*n))),
            nonbasic_columns: vec![1, 2],
        };
        table
    }
    fn pivot(&mut self) -> bool {
        let pivot_col = self.choose_var();
        let mut pivot_row = None;
        let mut candidate_ratio: Rational64 = Ratio::zero(); //Placeholder until set
        //Find the best ratio of nonbasic columns
        for col in 1..self.matrix.ncols() {
            if !self.nonbasic_columns.contains(&col) {
                for (index, val) in self.matrix.column(col).iter().enumerate() {
                    if val.is_integer() && val.numer() == &1 {
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
                println!("{}", self.matrix);
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
                println!("{}", self.matrix);
                //Find the nonbasic column we replaced
                for i in 0..self.nonbasic_columns.len() {
                    if self.nonbasic_columns[i] == pivot_col {
                        self.nonbasic_columns[i] = r;
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
        let mut best_col = self.nonbasic_columns[0];
        for col in self.nonbasic_columns.iter() {
            if self.matrix[(0, *col)] < self.matrix[(0, best_col)] {
                best_col = *col;
            }
        }
        return best_col;
    }

}

fn main() {
    use std::env;
    let args: Vec<String> = env::args().collect();
    let n = args[1].parse().unwrap();
    let m = args[2].parse().unwrap();
    let mut table = Tableau::new(n, m);
    //println!("{:?}", table.matrix[[0, 0]]);
    while table.pivot() {}
}

