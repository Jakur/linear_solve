extern crate num_traits;
extern crate num_rational;
extern crate nalgebra;

use num_traits::identities::Zero;
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

        let mut table = Tableau {
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
                                Some(r) => {
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
        println!("Pivot row: {} Pivot Column: {}", pivot_row.unwrap(), pivot_col);
        println!("Candidate ratio: {:?}", candidate_ratio);
        match pivot_row {
            Some(r) => {
                return true //placeholder
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
    table.pivot();
}

