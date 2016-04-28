#PythonCalculation
Python code for calculation

#GiniCorrelationCoefficient.py
Calculate Gini Correlation Coefficient (GCC). This is one speeded up version in Python of GCC from rsgcc Package  
(https://cran.r-project.org/web/packages/rsgcc/index.html). The running speed is improved by  
1) parallelization. and 2) pre-calculation of some statics per gene to avoid repetitive calculation in loops.

numpy, multipleprocessing, pandas modules are imported in this code.  

usage:  
python GiniCorrelationCoefficient.py    matrix   output    threadnum  
matrix : input matrix file with headline(sample IDs) and row names(gene IDs), each row contains data for each gene.  
output: output file  
threadnum: number of cores to be used  

