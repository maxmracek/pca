## Usage Instructions
Performs PCA on a given file, which defaults to "input_data.csv" if no argument is given

Input file should be a .csv file, with each row being the data associated with one sample. There needs to be more rows than columns or it will not work. 

After PCA is done, the results are written into a file called "pca_results.csv". The output file format is as follows:  
Line 1 is every eigenvalue of the covariance matrix.  
Line 2 is the amount of variance captured by nth principal axis.  
Line 3 is the cumulative amount of variance captured by the nth principal axis.  
The next lines are the principal axes of the data, written as row vectors, one for each eigenvalue.  
The final lines are the principal components, written as row vectors, with one for each eigenvalue.

## Notes
Errors may occur without proper information on how to fix them.   
The algorithm itself is inefficient and has room for improvement so it may have a long runtime with larger datasets.