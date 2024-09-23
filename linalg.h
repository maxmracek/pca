#include <vector>
#include <iostream>

//Subtracts the second vector from the first vector
std::vector<double> vecsub(std::vector<double> v1, std::vector<double> v2);

//Scales all elements in a vector by a given value
std::vector<double> scale(std::vector<double> v1, double s);

//Finds the dot product of two vectors
double dotprod(std::vector<double>v1,std::vector<double> v2);

//Finds the magnitude of the vector
double mag(std::vector<double> v1);

//Projects v1 onto v2
std::vector<double> proj(std::vector<double> v1, std::vector<double> v2);

//Multiplies two matrices with column vectors M1*M2
std::vector<std::vector<double>> colmatmult(std::vector<std::vector<double>> mat1,std::vector<std::vector<double>> mat2);

//Prints out Matrices defined by column vectors
void print2dcolmmat(std::vector<std::vector<double>> mat);