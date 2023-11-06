# Linear Hypothesis Testing for High Dimensional Tobit Models

The `TobitTest` package provides methods for computing partial penalized estimators and partial penalized Wald, score, and likelihood ratio test statistics for high dimensional Tobit models.

- `tobitADMM` fits a partial penalized Tobit model with a SCAD penalty using an alternating direction method of multipliers (ADMM) algorithm
- `pick_lambda_tobit` selects the penalty parameter for partial penalized Tobit models using a generalized information criterion (GIC)
- `compute_wald`, `compute_score`, and `compute_lrt` compute the partial penalized Wald, score, and likelihood ratio test statistics

See Jacobson and Zou (*accepted*) for details.

## References

Jacobson, T. and Zou, H. (*accepted*) "Linear Hypothesis Testing for High Dimensional Tobit Models," *Statistica Sinica*.
