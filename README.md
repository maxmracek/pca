## Usage Instructions
Performs PCA on a given file, which defaults to "input_data.csv" if no argument is given

Input file should be a .csv file, with each row being the data associated with one sample. There needs to be more rows than columns or it will not work. 

After PCA is done, the results are written into two different files. The files are named PREFIX_results.csv and PREFIX_new_vectors.csv. PREFIX defaults to "pca", but otherwise will be the second argument in the command.

In PREFIX_results.csv:
Line 1 is every eigenvalue of the covariance matrix.  
Line 2 is the amount of variance captured by nth principal axis.  
Line 3 is the cumulative amount of variance captured by the nth principal axis.  
The next lines are the principal components of the data, written as row vectors, one for each eigenvalue.  

In PREFIX_new_vectors.csv:
Each line corresponds to the sample at the same line number in the input data, except transformed to be written using the principal components as their basis. The points in the new vectors are written in the order of most to least variance captured by the respective principal component. 

## Notes
Errors may occur without proper information on how to fix them.   
The algorithm itself is inefficient and has room for improvement so it may have a long runtime with larger datasets.
