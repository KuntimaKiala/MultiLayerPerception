#ifndef SIMPLENEURALNETWORK_HPP
#define SIMPLENEURALNETWORK_HPP

#pragma once
#include "Topology.hpp"
#include "Matrix.hpp"
#include <cassert>
#include <string>
//typedef std::vector<std::vector<Neurons>> Layers ;
//typedef std::vector<Neurons*> Layers ;
class SimpleNeuralNetWork
{
public:
  
    SimpleNeuralNetWork();
    SimpleNeuralNetWork(Topology& toposlogy);

    void train(std::vector<Matrix>& input_data, std::vector<Matrix>& output_data);
    void FeedForward();
    static int predit(const double &y_hat) ;
    void backPropagation();
    void calcErrors(Matrix& output);
    void updateWeights();
    void prediction(std::vector<Matrix>& input_data, std::vector<Matrix>& output_data) ;

    double sigmoid (const double& m1);
    Matrix cost(Matrix& y, Matrix& y_hat) ;
    inline double sigmoid (double& m1);
    inline double sigmoid_deriv (double& m1);
    inline double ReLu(double & x);
    inline double ReLu_deriv (double& y);

    ~SimpleNeuralNetWork();

    std::vector<Matrix > weights ; 
    std::vector<Matrix > bias ; 
    
    Matrix X ;
    Matrix Y ;
    std::vector<Matrix> errors;
    int epochs = 0 ;
    double learningRate = 0.2;
    Layers layers ;
  
    
    
private:
Topology * _topology ;
};

#endif