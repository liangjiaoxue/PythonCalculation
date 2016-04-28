# PythonCalculation
Python code for calculation

#GiniCorrelationCoefficient.py
Calculate Gini Correlation Coefficient (GCC). This speeded up version of GCC on rsgcc R Package 
(https://cran.r-project.org/web/packages/rsgcc/index.html). The running speed is improve by 
1) parallelization. and 2) precalculation of some statics per gene to avoid repeative calculation in loops.

numpy, multipleprocessing, pandas modules are imported in this code.



