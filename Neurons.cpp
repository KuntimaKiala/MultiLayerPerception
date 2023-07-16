#include "Neurons.hpp"

Neurons::Neurons()
{

}


Neurons::Neurons(uint32_t rows, uint32_t columns){
    Matrix _neurons(rows, columns, false) ;
    neurons = _neurons;

}


Neurons::~Neurons()
{

  
}

int Neurons::size() {
return neurons.columns();

}


