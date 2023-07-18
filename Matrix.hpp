#ifndef MATRIX_H
#define MATRIX_H

#pragma once
#include <vector>
#include <cstdint>
#include<iostream>
#include <cassert>
#include<stdlib.h>
#include<time.h>
#include <cmath>
#include <functional>



class Matrix{
    public:
        
        Matrix();
        Matrix(uint32_t rows, uint32_t columns,bool rand_init = true);
        Matrix(std::vector<std::vector<double>> MAT);
        double at(uint32_t rows, uint32_t columns);
        Matrix operator+(Matrix &m);
        Matrix operator-(Matrix &m);
        Matrix operator+(const double &m);
        Matrix operator-(const double &m);
        Matrix operator*(Matrix &m) ;
        Matrix operator*(const double &m) ;
        Matrix DiagonalMatrix() ;
        Matrix  Convolution2D(Matrix &matrix, Matrix &kernel, int padding =0, int stride =1, bool vis =false ) ;
        Matrix  block(int row, int column, int w, int h, int kernel_size, bool vis = false) ;
        Matrix  padding(int padding) ;
        double  sum() ;
        Matrix Hadamard(Matrix &m1);
        Matrix transpose() ;
        Matrix sigmoid () ;
        Matrix sigmoid_deriv ();
        Matrix ReLu(Matrix & x);
        Matrix ReLu_deriv (Matrix& y);
        //template<typename T>
        Matrix ApplyFunction(std::function <double(const double & m)> func);
        Matrix xavier_init() ;
        void shape() ;
        int columns() ;
        int rows() ;
        ~Matrix();
        std::vector<std::vector<double>> element;

    private :
    uint32_t _rows ;
    uint32_t _columns ;
    std::vector<std::vector<double>> _MAT;
    std::vector<std::vector<double>> _diag;
    

};


#endif