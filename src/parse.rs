use super::*;
use num::FromPrimitive;

#[derive(Debug, PartialEq)]
pub enum FileFormat {
    Custom, //Our format
    MPS,
}

pub enum TabularData {
    LP(Vec<Rational64>, Vec<Vec<Rational64>>, Vec<usize>), //Opt, Constraints, Equality Indices
    //x_count, Opt, Constr, Eq Indices
    ParameterLP(usize, Vec<Rational64>, Vec<Vec<Rational64>>, Vec<usize>),
}

///Parses a matrix as written in text as a vector of row vectors
pub fn parse_inequalities(text: &str) -> TabularData {
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

pub fn parse_mps(text: &str) -> TabularData {
    use std::collections::HashMap;

    let mut row_names = HashMap::new();
    let mut cols = HashMap::new();
    let mut lines = text.lines();
    lines.next();
    lines.next();
    let mut greater = Vec::new();
    let mut equals = Vec::new();
    let mut line = lines.next().unwrap();
    let mut row = 0;
    let mut optim_row = 0;
    while !line.starts_with("COLUMNS") { //Read Rows
        let mut split = line.split_whitespace();
        let first = split.next().unwrap();
        match first {
            "E" => {equals.push(row); greater.push(false);},
            "G" => {greater.push(true);},
            "N" => { //Should represent the optimization row
                greater.push(false);
                optim_row = row;
            }
            _ => {greater.push(false);}, //Should only be "L"
        }
        let second = split.next().unwrap();
        row_names.insert(second, row);

        row += 1;
        line = lines.next().unwrap();
    }
    line = lines.next().unwrap();
    while !line.starts_with("RHS") { //Read Cols
        let split: Vec<_> = line.split_whitespace().collect();
        if split.len() >= 3 {
            let num_cols = cols.len();
            let (ref mut entry, _) = cols.entry(split[0])
                .or_insert((Vec::new(), num_cols));
            let value1: f64 = split[2].parse().unwrap();
            entry.push((row_names.get(split[1]).unwrap(), value1));
            if split.len() >= 5 {
                let value2: f64 = split[4].parse().unwrap();
                entry.push((row_names.get(split[3]).unwrap(), value2));
            }
        }
        line = lines.next().unwrap();
    }
    line = lines.next().unwrap();
    let mut last_col = Vec::with_capacity(row_names.len());
    while !line.starts_with("ENDATA") { //Todo bounds, etc.
        let split: Vec<_> = line.split_whitespace().collect();
        if split.len() >= 3 {
            if split.len() % 2 == 1 { //This mps redundantly repeats the name of the z col
                let value1: f64 = split[2].parse().unwrap();
                last_col.push((row_names.get(split[1])
                                   .expect(&format!("Couldn't find row {}", split[1])), value1));
                if split.len() >= 5 {
                    let value2: f64 = split[4].parse().unwrap();
                    last_col.push((row_names.get(split[3]).unwrap(), value2));
                }
            } else { //Format which omits name of z col
                let value1: f64 = split[1].parse().unwrap();
                last_col.push((row_names.get(split[0])
                                   .expect(&format!("Couldn't find row {}", split[0])), value1));
                if split.len() >= 4 {
                    let value2: f64 = split[3].parse().unwrap();
                    last_col.push((row_names.get(split[2]).unwrap(), value2));
                }
            }
        }
        line = lines.next().unwrap();
    }
    let num_cols = cols.len() + 1; //+1 for incoming Z_Col
    let num_rows = row_names.len();
    cols.insert("Z_Col", (last_col, num_cols - 1));
    let mut tuples: Vec<_> = cols.into_iter().flat_map(|(_key, (col, col_index))| {
        col.into_iter().map(move |(row_num, value)| {
            (*row_num, col_index, value)
        })
    }).collect();
    tuples.sort_by(|a, b| a.0.cmp(&b.0));
    println!("Size: {}", tuples.len());
    let mut tup_iter = tuples.into_iter();
    let mut tuple = tup_iter.next();
    let mut opt = None;
    let mut constraints = Vec::with_capacity(num_rows);

    for i in 0..num_rows {
        let mut row_vec = vec![Ratio::zero(); num_cols];
        while let Some(t) = tuple {
            if t.0 != i {
                break;
            }
            if greater[i] {
                row_vec[t.1] = Ratio::from_f64(-1.0 * t.2)
                    .expect("Invalid Rational");
            } else {
                row_vec[t.1] = Ratio::from_f64(t.2)
                    .expect("Invalid Rational");
            }
            tuple = tup_iter.next();
        }
        if i == optim_row {
            opt = Some(row_vec);
        } else {
            constraints.push(row_vec);
        }
    }
    let mut nonzero = 0;
    for row in constraints.iter() {
        print!("[ ");
        for val in row.iter().rev().skip(1) {
            if *val != Ratio::zero() {
                nonzero += 1;
            }
            print!("{} ", val);
        }
        println!("] ");
    }
    let opt = opt.unwrap();
    for val in opt.iter() {
        if *val != Ratio::zero() {
            nonzero += 1;
        }
    }
    let opt = Some(opt);
    println!("Nonzero: {}", nonzero);
    return TabularData::LP(opt.expect("No optimization row found"), constraints, equals);
}

pub fn create_table(text: &str, format: FileFormat) -> (Box<Tableau>, bool) {
    let parsed = match format {
        FileFormat::Custom => parse_inequalities(text),
        FileFormat::MPS => parse_mps(text),
    };
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