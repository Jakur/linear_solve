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
            basic_columns: vec![1, 2],
        };
        for x in arr.iter().enumerate() {
            for y in x.1.iter().enumerate() {
                table.matrix[[x.0, y.0]] = Ratio::from_integer(*y.1);
            }
        }
        table
    }
    fn pivot(&mut self) {
        let var = self.choose_var();
        //Find the best ratio of nonbasic columns
        for col in 1..self.matrix.cols() {
            if !self.basic_columns.contains(&col) {
                for (index, val) in self.matrix.column(col).iter().enumerate() {
                    if val.is_integer() && val.numer() == &1 {
                        let pivot_value = self.matrix[[index, var]];
                        if pivot_value > Ratio::zero() {
                            let ratio = self.matrix[[index, self.matrix.cols()-1]] / pivot_value;
                            println!("{:?}", ratio);
                        }
                    }
                }
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
    let mut table = Tableau::new(5, 8);
    //println!("{:?}", table.matrix[[0, 0]]);
    table.pivot();
}

