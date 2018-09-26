extern crate num_traits;
extern crate num_rational;
extern crate ndarray;

use num_traits::identities::Zero;
use num_rational::Rational64;
use num_rational::Ratio;
use ndarray::Array2;


struct Tableau {
    matrix: Array2<Rational64>,
    basic_columns: Vec<usize>,
}

impl Tableau {
    fn new(n: usize, m: usize) -> Tableau {
        //Debug testing array
        let arr = [[1, -4, -3, 0, 0, 0, 0, 0],
            [0, 2, 3, 1, 0, 0, 0, 6],
            [0, -3, 2, 0, 1, 0, 0, 3],
            [0, 0, 2, 0, 0, 1, 0, 5],
            [0, 2, 1, 0, 0, 0, 1, 4]];

        let mut table = Tableau {
            matrix: Array2::zeros((n, m)),
            basic_columns: vec![0, 0, 0, 1, 1, 1, 1],
        };
        for x in arr.iter().enumerate() {
            for y in x.1.iter().enumerate() {
                table.matrix[[x.0, y.0]] = Ratio::from_integer(*y.1);
            }
        }
        table
    }
    fn pivot(&mut self) -> bool {
        let pivot_col = self.choose_var();
        let mut pivot_row = None;
        let mut candidate_ratio = Ratio::zero(); //Placeholder until set
        //Find the best ratio of nonbasic columns
        for col in 1..self.matrix.cols() {
            if self.basic_columns[col] == 1 {
                for (index, val) in self.matrix.column(col).iter().enumerate() {
                    if val.is_integer() && val.numer() == &1 {
                        let pivot_value = self.matrix[[index, pivot_col]];
                        if pivot_value > Ratio::zero() {
                            let ratio = self.matrix[[index, self.matrix.cols()-1]] / pivot_value;
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
        let mut best_col = self.basic_columns[0];
        for col in self.basic_columns.iter() {
            if self.matrix[[0, *col]] < self.matrix[[0, best_col]] {
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

