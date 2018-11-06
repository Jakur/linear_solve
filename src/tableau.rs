use num_traits::identities::Zero;
use num_traits::identities::One;
use num_rational::Rational64;
use num_rational::Ratio;
use nalgebra::DMatrix;

pub trait Tableau {
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
        let mut chosen = vec![false; self.matrix().nrows()];
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
                            if !chosen[index] {
                                row = Some(index);
                            }
                        }
                    }
                }
            }
            if let Some(r) = row {
                values.push(self.matrix()[(r, self.matrix().ncols()-1)]);
                chosen[r] = true;
            } else { //The whole column was zero. Unknown if this can happen
                values.push(Ratio::zero());
            }
        }
        return values;
    }
    fn artificial_cols(&self) -> &Vec<bool>;
}