COMPILE:

	cargo build --release


## The lp format is
## objective (to maximize)
## Coefficients, direction, rhs

TEST:
	cargo run test0.lp

TEST-A:
	cargo run test-a.lp


TEST-B:
	cargo run test-b.lp

TEST-C:
	cargo run test-c.lp
