#include "linalg.h"
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char* argv[]){
    //Input Data
    std::string filename = "input_data.csv";
    if(argc != 1){
        filename = argv[1];
    }
    std::ifstream myfilein{filename};

    std::string inputLine;
    std::vector<std::vector<double>> mat1;
    std::getline(myfilein,inputLine);
    std::stringstream ss0(inputLine);
    std::string val0;
    while(std::getline(ss0,val0,',')){
        mat1.push_back({std::stod(val0)});
    }
    while(std::getline(myfilein,inputLine)){
        int vars = 0;
        std::stringstream ss(inputLine);
        std::string val;
        while(std::getline(ss,val,',')){
            mat1[vars].push_back(std::stod(val));
            vars++;
        }
    }
    if(mat1.size() >= mat1[0].size()){
        std::cout<< "Input needs to have more rows than colums" << std::endl;
        return 0;
    }

    //Find Mean Shifted Matrix
    std::vector<std::vector<double>> meanMat; //each nested vector is one variable (column), len of vector is number of samples
    for (int i = 0; i< mat1.size(); i++){ //for each variable/column
        double mean = 0;
        for (int j = 0; j < mat1[0].size();j++){ //find the mean for each sample for the variable
            mean+=mat1[i][j];
        }
        mean = mean/mat1[0].size();
        std::vector<double> meanvec;
        for (int j = 0; j < mat1[0].size();j++){
            meanvec.push_back(mat1[i][j] - mean);  
        }
        meanMat.push_back(meanvec); //add adjusted column to new matrix
    }

    //Covariance Matrix Calculation
    //Does something like backwards transposed matrix multiplication
    std::vector<std::vector<double>> covMat;
    for(int i = 0; i< meanMat.size();i++){ 
        std::vector<double> sumVec;
        for(int j = 0; j < meanMat.size();j++){ 
            double sum = 0;    
            for(int k = 0; k<meanMat[0].size();k++){ 
                sum += meanMat[i][k] * meanMat[j][k]; 
            }
            sumVec.push_back(sum);
        }
        sumVec = scale(sumVec,(1.0/meanMat[0].size()));
        covMat.push_back(sumVec);
    } 
    
    //QR Algorithm
    std::vector<std::vector<double>> newA = covMat;
    std::vector<std::vector<double>> Q;
    std::vector<std::vector<double>> R;
    std::vector<std::vector<double>> eigenvectors;
    std::vector<std::vector<double>> oldA = {};
    for(int i = 0; i < newA.size();i++){
        oldA.push_back(std::vector<double> (newA[i].size(),0));
    }
    int loopnum = 0;
    int maxloop = 10001; //has errors if theres an eigenvalue equal to 0
    while(loopnum < maxloop){ 
        Q = {};
        R = {};
        loopnum+=1;
        std::vector<double> u1 = {newA[0]};
        std::vector<double> e1;
        std::vector<std::vector<double>> U;
        U.push_back(u1);
        double u1mag = mag(U[0]);
        e1 = scale(u1,1/u1mag);
        Q.push_back(e1);
        std::vector<double> r1 (newA.size(),0.0);
        r1[0] = u1mag;
        R.push_back(r1);

        //calculates Q and R matrices
        for(int k = 1;k<newA.size();k++){ //For each column after the first
            std::vector<double> uvec;
            std::vector<double> avec = newA[k];
            std::vector<double> prevvec = avec;
            std::vector<double> nextvec;
            for (int i = 0; i < k; i++){
                nextvec = vecsub(prevvec,proj(prevvec,U[i]));
                prevvec = nextvec;
            }
            U.push_back(nextvec);
            Q.push_back(scale(U[k],1/mag(U[k]))); // adds e vector to Q
            std::vector<double> rvec(newA[0].size(),0.0); 
            for(int i = 0; i < k+1; i++){
                rvec[i] = dotprod(Q[i],avec);
            }  
            R.push_back(rvec);
        }

        //Calculates Eigenvectors Iteratively
        if(loopnum == 1){
            eigenvectors = Q;
        }
        else{
            eigenvectors = colmatmult(eigenvectors,Q);
        }

        //Finds R*Q for the next loop
        newA = colmatmult(R,Q);

        if(std::isnan(newA[0][0])){
            std::cout << "Correlation Matrix is singular" << std::endl;
            return 0;
        }

        //Early end if it converges
        if(loopnum%100 == 0){
            double cusum = 0;
            for(int i = 0; i < newA.size();i++){
                for(int j = 0; j < newA.size();j++){
                    cusum += (fabs(newA[i][j])-fabs(oldA[i][j]));
                }
            }
            if(cusum <= 1e-100 && maxloop-loopnum > 100){
                loopnum = maxloop - 50;
            }
            else{
                oldA = newA;
            }
        }
    }

    std::vector<double> eigenvalues;
    for(int i = 0; i < newA.size() && i < newA[0].size();i++){
        eigenvalues.push_back(newA[i][i]);
    }

    double eigensum = 0.0;
    for(int i = 0; i < eigenvalues.size();i++){
        eigensum += fabs(eigenvalues[i]);
    }

    std::vector<double> scaledeigenvalues;
    for(int i = 0;i<eigenvalues.size();i++){
        scaledeigenvalues.push_back(fabs(eigenvalues[i])/eigensum);
    }

    std::string outputfile;
    if(argc == 3){
        outputfile = argv[2];
    }
    else{
        outputfile = "pca";
    }

    std::ofstream myfile;
    //file output: line of eigenvalues, line of % variance captured, line of cumulative % captured
    //eigenvectors as row vectors
    //samples as row vectors with principal components as basis in the same order as above
    myfile.open(outputfile + "_results.csv");

    myfile << eigenvalues[0];
    for(int i = 1; i < eigenvalues.size();i++){
        myfile << ", " << eigenvalues[i];
    }
    myfile << std::endl;

    myfile << scaledeigenvalues[0];
    for(int i = 1; i < scaledeigenvalues.size();i++){
        myfile << ", " <<scaledeigenvalues[i];
    }
    myfile << std::endl;

    double cueigensum = scaledeigenvalues[0];
    myfile << cueigensum;
    for(int i = 1; i < scaledeigenvalues.size();i++){
        cueigensum+= scaledeigenvalues[i];
        myfile << ", " << cueigensum; 
    }
    myfile << std::endl;

    for(int i = 0; i < eigenvectors.size();i++){
        myfile << eigenvectors[i][0];
        for(int j = 1; j < eigenvectors[i].size();j++){
            myfile << ", " << eigenvectors[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    std::ofstream myfile2;
    myfile2.open(outputfile + "_new_vectors.csv");
    for(int i = 0; i < meanMat[0].size();i++){
        std::vector<double> sample;
        for(int j = 0; j < meanMat.size();j++){
            sample.push_back(meanMat[j][i]);
        }
        std::vector<double> projection;
        for(int j = 0; j < eigenvectors.size()-1;j++){
            myfile2 << dotprod(sample,eigenvectors[j]) << ',';
        }
        myfile2 << dotprod(sample,eigenvectors[eigenvectors.size()-1]);
        if(i != meanMat[0].size()-1){
            myfile2 << std::endl;
        }
    }
    myfile2.close();

    return 0;
}