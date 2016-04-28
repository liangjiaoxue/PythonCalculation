#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Calculate Gini Correlation Coefficient (GCC). This speeded up version of GCC on rsgcc R Package
(https://cran.r-project.org/web/packages/rsgcc/index.html). The running speed is improved by
1) parallelization. and 2) pre-calculation of some statics per gene to avoid repetitive calculation in loops.

numpy, multipleprocessing, pandas modules are imported in this code.
usage:
python GiniCorrelationCoefficient.py matrix  output threadnum
matrix : input matrix file with headline(sample IDs) and row names(gene IDs), each row contains data for each gene.
output: output file
threadnum: number of cores to be used
'''

from sys import argv
import numpy as np
import multiprocessing as mp
import pandas as pd


code, matrix_in,  output_file ,thread_num0 = argv

thread_num = int(thread_num0)
## read data into array
with open(matrix_in,"rU") as textFile:
    lines = [line.rstrip().split("\t") for line in textFile]

array = np.array(lines[1:]) # remove head line
gene = array[:,0]
data = (array[:,1:]).astype(np.float) # remove first column

# sort array
data_sorted = np.apply_along_axis(np.sort,axis=1, arr=data)
data_order  = np.apply_along_axis(np.argsort,axis=1, arr=data)

def self_sum (vector ) :
    len_x = len(vector)
    sum_my = 0
    for x_index in range(0,len_x):
        sum_my += (2*(x_index+1) - len_x -1)* vector[x_index]
    return sum_my

## pre-calculation of sum to avoid repetitive calculation
data_self_sum =  np.apply_along_axis(self_sum,axis=1, arr= data_sorted)



# data partition
shape = data.shape
data_num = shape[1]
gene_num = shape[0]

index2range = {}
part_num = int(gene_num/thread_num)

for index in range(thread_num) :
    start_in = index* part_num
    end_in   = (index+1)* part_num
    if index == thread_num -1 :
        end_in = gene_num
    index2range[index] = range(start_in,end_in)


def worker(index_in, list_in):
    range_in = index2range[index_in]
    #print(range_in)
    array2d =  np.empty((len(range_in),gene_num), float)
    for x_pos in range_in :
        for y_pos in range(gene_num) :
            sum_x_self = data_self_sum[x_pos]
            y_order = data_order[y_pos,:]
            x_data  = data[x_pos,:]
            x_data_by_y = x_data[y_order]
            sum_x_by_y = self_sum(x_data_by_y)
            out = sum_x_by_y/sum_x_self
            array2d[x_pos - part_num*index_in  ,y_pos] = out
    list_in[index_in] = array2d


if __name__ == '__main__':
    manager = mp.Manager()
    return_list = manager.list(range(thread_num))
    jobs = []
    for i in range(thread_num):
        p = mp.Process(target=worker, args=(i,return_list))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()


out =  np.empty((0,gene_num), float)
for i in range(thread_num):
    array_out = return_list[i]
    out = np.append(out,array_out, axis=0)


# maximum value
for i in range(gene_num) :
    for j in range(i):
        if abs(out[i,j]) < abs(out[j,i]) :
            out[i, j] = out[j, i]
        else :
            out[j, i] = out[i, j]


df = pd.DataFrame(out, index=gene, columns=gene)
df.to_csv(output_file, index=True, header=True, sep="\t")
### write out files


