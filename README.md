# linear_solve
A solver for parametric linear programs written in Rust. Specifically, it solves programs of the form maximize<br/>   
λ<sub>1</sub>(c<sub>1</sub>x<sub>1</sub> + c<sub>2</sub>x<sub>2</sub> + ... c<sub>0</sub>) + 
λ<sub>2</sub>(c<sub>1</sub>x<sub>1</sub> + c<sub>2</sub>x<sub>2</sub> + ... c<sub>0</sub>) + ... 
and prints resulting inequalities in x. <br/>
## Input
Input is a file in the root project directory. The file must fit the form described below.
### Optimization Line
Input consists of the optimization function on the first line then constraints on all further lines. 
All terms are separated by whitespace.
Only the coefficients of the terms are supplied by the user.
Parametric input has terms consisting of polynomials in x, so to indicate they are one term 
the polynomials are enclosed in square brackets, e.g. <br/> 
[1 0 -4] [2 3 0] corresponds to λ<sub>1</sub>(x<sub>1</sub> - 4) + λ<sub>2</sub>(2x<sub>1</sub> + 3x<sub>2</sub>)<br/>
### Constraint Lines
In the parametric case, constraints are equalities or inequalities in lambda. 
Constant terms are put on the right hand side of the equality or inequality, which can be either "=", ">=" or "<="
Additionally, there is an implicit constraint that all variables have non-negative values, i.e. λ<sub>i</sub> ≥ 0. <br/>
