#include <iostream>
#include "Matrix.hpp"
#include <random>
#include <chrono>
#include <cmath>


Matrix::Matrix()
{

}



Matrix::Matrix(std::vector<std::vector<double>> MAT) : element(MAT), _MAT(MAT){

 
    _rows = MAT.size();
    _columns = MAT[0].size();
}


Matrix::Matrix(uint32_t rows, uint32_t columns, bool rand_init )   : _rows(rows), _columns(columns), _MAT({ }){
_MAT.resize(_rows, std::vector<double>(_columns, 0.0) );
element.resize(_rows, std::vector<double>(_columns, 0.0) );


if(rand_init == true){
    Matrix m = xavier_init() ;
  
    
    
    
    for (uint32_t row = 0; row<  m.rows(); row++){

        
        for (uint32_t col = 0; col <  m.columns() ; col++){
            
            _MAT[row][col]    =  (double) m.element[row][col];
            element[row][col] =  _MAT[row][col];
         
        }
       
    }

}



}


Matrix Matrix::xavier_init() {

Matrix output(_rows, _columns, false) ;
double xavier_stddev = std::sqrt(2.0/(_rows + _columns));
int seed = std::chrono::system_clock::now().time_since_epoch().count();


std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0,xavier_stddev);


for (uint32_t row = 0; row< this->_MAT.size(); row++){
        for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){

            output._MAT[row][col]    = distribution(generator);
            output.element[row][col] = output._MAT[row][col] ;
            
        }
    }

return output;

}

Matrix Matrix::ApplyFunction(std::function <double(const double & m)> func){

    Matrix output(_rows, _columns, false) ;
    
    for (uint32_t row = 0; row< this->_MAT.size(); row++){
        for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){

            output._MAT[col][row]    = func(at(row,col)) ;
            output.element[col][row] = func(at(row,col)) ;
        }
    }
    return output;
}



Matrix Matrix::DiagonalMatrix(uint32_t rows, uint32_t columns) {
Matrix Diag(rows, columns, false);
for (uint32_t row = 0; row< this->_MAT.size(); row++){
        for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){

            if(row == col) {
                
                Diag._MAT[row][col] = _MAT[row][col];
                Diag.element[row][col] = element[row][col];
                
            }
        }

}

return Diag;

}






Matrix Matrix::operator-(Matrix &m){
    assert( _columns == m._columns && "different of numbers of columns matrix_summation");
    assert( _rows    == m._rows    && "different of numbers of rows matrix_summation");

    Matrix result = Matrix(_rows, m._columns, false) ;

    for (uint32_t row = 0; row< this->_MAT.size(); row++){
        for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){
      
         result._MAT[row][col]    = this->_MAT[row][col] - m._MAT[row][col] ;
         result.element[row][col] = this->_MAT[row][col] - m._MAT[row][col] ;
        
        }
    
    }
    return result;

}


Matrix Matrix::operator+(Matrix &m){
    assert( _columns == m._columns && "different of numbers of columns matrix_summation");
    assert( _rows    == m._rows    && "different of numbers of rows matrix_summation");

    Matrix result = Matrix(_rows, m._columns) ;

    for (uint32_t row = 0; row< this->_MAT.size(); row++){
        for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){
      
         result._MAT[row][col] = this->_MAT[row][col] + m._MAT[row][col] ;
         result.element[row][col] = this->_MAT[row][col] + m._MAT[row][col] ;
        }
    
    }
    return result;

}

Matrix Matrix::operator+(const double &m){
    Matrix result = Matrix(_rows, _columns) ;
    for (uint32_t row = 0; row< this->_MAT.size(); row++){
            for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){
                result._MAT[row][col] =    this->_MAT[row][col] + m ;
                result.element[row][col] = this->_MAT[row][col] + m;
           
            }
    }
    return result;
}

Matrix Matrix::operator-(const double &m){
    Matrix result = Matrix(_rows, _columns) ;
    for (uint32_t row = 0; row< this->_MAT.size(); row++){
            for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){
                result._MAT[row][col] =    this->_MAT[row][col] - m ;
                result.element[row][col] = this->_MAT[row][col] - m;
           
            }
    }
    return result;
}


Matrix Matrix::Hadamard(Matrix & m1){

   
    assert( _columns == m1._columns && "number of row different of numbers of columns Hadamard");
    assert( _rows    == m1._rows && "number of columns different of numbers of rows Hadamard");

    Matrix result = Matrix(_rows, _columns) ;

    for (uint32_t row = 0; row < m1.rows(); row++){
         for (uint32_t col = 0; col < m1.columns(); col++){
           result._MAT[row][col] = m1.at(row,col) * this->at(row,col);
           result.element[row][col] = m1.at(row,col) * this->at(row,col);
           
           
        }


    }

    return result ;
}



Matrix Matrix::operator*(Matrix &m){


    assert( _columns == m._rows && "number of row different of numbers of columns operator*");
    
   
    
    Matrix result = Matrix(_rows, m._columns) ;
    
    for (uint32_t row_1 = 0; row_1 < this->_MAT.size(); row_1++){
        for (uint32_t col_2 = 0; col_2 < m._columns ; col_2++){
            result._MAT[row_1][col_2] = 0.0;
            for (uint32_t row_2 = 0; row_2 < m._rows ; row_2++) {
                result._MAT[row_1][col_2] += this->at(row_1, row_2) * m.at(row_2, col_2) ;
                result.element[row_1][col_2] += this->at(row_1, row_2) * m.at(row_2, col_2) ;

            }
        }

    }

    
    return result ;

}


Matrix Matrix::operator*(const double &m){
    Matrix result = Matrix(_rows, _columns) ;
    for (uint32_t row = 0; row< this->_MAT.size(); row++){
            for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){
                result._MAT[row][col] = this->_MAT[row][col]* m;
                result.element[row][col]  = this->element[row][col] * m;
            }

    }

    return result ;
}


double Matrix::at(uint32_t rows, uint32_t columns){

return _MAT[rows][columns];
}


Matrix Matrix::transpose() { 

    Matrix transposed(_columns, _rows, false);
    for (uint32_t row = 0; row< this->_MAT.size(); row++){
        for (uint32_t col = 0; col < this->_MAT[0].size() ; col++){

            transposed._MAT[col][row] = this->_MAT[row][col] ;
            transposed.element[col][row] = this->_MAT[row][col] ;
        }

}

return transposed;

}


Matrix Matrix::sigmoid () {
    
    /*  Returns the value of the sigmoid function f(x) = 1/(1 + e^-x).
        Input: m1, a vector.
        Output: 1/(1 + e^-x) for every element of the input matrix m1.
    */
   
   
    Matrix output(_rows,_columns,false);
 
 
   
    
    for( uint32_t r = 0; r < _rows; r++ ) {
        for( uint32_t c = 0; c < _columns; c++ ) {
            //m1.element[r][c] = m1.element[r][c] * (-1) ;
            //std::cout << "sig :" <<m1.element[r][c] <<"\n" ;
            output._MAT[r][c] = (double)(1.0 / (1.0 + (double)std::exp(_MAT[r][c]*(-1))));
            output.element[r][c] = (double)(1.0 / (1.0 + (double)std::exp(-element[r][c]*(-1))));
        }
    }
    
 
    return output;
}

Matrix Matrix::sigmoid_deriv (){


    Matrix output(_rows, _columns);
    Matrix output_minus_one(_rows, _columns);
    Matrix result(_rows, _columns);
    output = output.sigmoid() ;
    output_minus_one = output ;
    output_minus_one = output_minus_one*(-1) + 1 ;
   
    result = output.Hadamard(output_minus_one) ;
  
 

    return result;

}

Matrix Matrix::ReLu(Matrix & x){

   
 
    Matrix A(x.rows(), x.columns());


    for( int i = 0; i< x.rows(); ++i )
    for (int j=0; j< x.columns(); ++j) {
        if (x.at(i,j) <= 0){
            A.element[i][j]=0.0;
        }
        else A.element[i][j]=x.at(i,j);
    }
    return std::move(A);

}

Matrix Matrix::ReLu_deriv (Matrix& y) {
    Matrix B(y.rows(), y.columns());
    for( int i = 0; i < y.rows(); ++i )
    for (int j=0; j < y.columns() ; ++j)
    {
        if (y.at(i,j) <= 0.0){
            B.element[i][j]=0.0;
        }
        else B.element[i][j]=1.0;
    }
    return std::move(B);
}



void Matrix::shape() {

    std::cout << "Matrix Shape : (" << _MAT.size() << "," << _MAT[0].size() << ")" << std::endl;

}



int Matrix::columns() {
    return _MAT[0].size() ;
}


int Matrix::rows() {
return _MAT.size() ;

}




Matrix::~Matrix()
{

}