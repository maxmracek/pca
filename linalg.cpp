#include "linalg.h"
#include <cmath>

std::vector<double> vecsub(std::vector<double> v1,std::vector<double> v2){
    if(v1.size()!=v2.size()){
        std::cout<< "Vectors Different Lengths";
        return {0};
    }
    std::vector<double> diffmat; 
    for(int i = 0; i < v1.size();i++){
        diffmat.push_back(v1[i]-v2[i]);
    }
    return diffmat;  
}

std::vector<double> scale(std::vector<double> v1, double s){
    std::vector<double> result;
    for(int i = 0;i<v1.size();i++){
        result.push_back(v1[i]*s);
    }
    return result;
}

double dotprod(std::vector<double> v1, std::vector<double> v2){
    if(v1.size()!=v2.size()){
        std::cout<< "Vectors are not the same size";
        return 0;
    }
    double dotprod = 0;
    for(int i = 0; i < v1.size();i++){
        dotprod+= v1[i]*v2[i];
    }
    return dotprod;
}

double mag(std::vector<double> v1){
    double magnitude = dotprod(v1,v1);
    return std::sqrt(magnitude);
}

std::vector<double> proj(std::vector<double> v1, std::vector<double> v2){
    double scalar = dotprod(v1,v2)/dotprod(v2,v2);
    //this checks for same len alrdy
    std::vector<double> projection;
    for(int i = 0;i < v2.size();i++){
        projection.push_back(v2[i]*scalar);
    }
    return projection;
}

std::vector<std::vector<double>> colmatmult(std::vector<std::vector<double>> mat1, std::vector<std::vector<double>> mat2){
    std::vector<std::vector<double>> prodmat;
    for(int i = 0; i < mat1.size(); i++){
        std::vector<double> newcol = {};
        for (int j = 0; j < mat1[i].size();j++){
            double colsum = 0;
            for(int k = 0; k <mat2[0].size();k++){
                colsum+=mat1[k][j] * mat2[i][k];
            }
            newcol.push_back(colsum);
        }
        prodmat.push_back(newcol);
    }
    return prodmat;
}

void print2dcolmmat(std::vector<std::vector<double>> mat){
    for(int i = 0; i < mat[0].size();i++){
        for(int j = 0; j < mat.size();j++){
            std::cout<< mat[j][i] << " ";
        }
        std::cout << std::endl;
    }
}