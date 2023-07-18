#ifndef NEURONS_HPP
#define NEURONS_HPP

#pragma once


#include <cmath>
#include "Matrix.hpp"
class Neurons
{
public:
    Neurons();
    Neurons(uint32_t rows, uint32_t columns);
    ~Neurons();
    int size() ;



    Matrix neurons ;
    
    
};

#endif