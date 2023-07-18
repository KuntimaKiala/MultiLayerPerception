#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include <bits/stdc++.h>

using namespace std ;




int main(){
vector<vector<double>> in{{1,2,8,6}, 
                          {3,7,5,4},
                          {4,2,8,9},
                          {3,5,6,7}};

    
    /*vector<vector<double>> out{{-1,0,1}, 
                               {-1,0,1},
                               {-1,0,1}};


    /*vector<vector<double>> in{{1,2}, 
                              {2,8}};*/

    vector<vector<double>> out{{1,0}, 
                               {0,1}};
                             
    
    Matrix mat(in) ;
    Matrix kernel(out) ;
    cout << kernel.rows() << endl ;
    cout << kernel.columns() << endl ;
    Matrix matrix ;
    int padding =0;
    int stride =1;
    if (stride == 0) stride = 1 ;
    Matrix output = matrix.Convolution2D(mat,kernel,padding,stride) ;


    std::cout << "\n" ;
    
    std::cout << "Ouput :\n";
    for(int i = 0; i < output.rows(); i++){
        for(int j = 0; j < output.columns(); j++){

            std::cout << output.at(i,j) << " ";
            
        }
        std::cout << "\n";

    }

    std::cout << "\n";

    
    return 1 ;
}