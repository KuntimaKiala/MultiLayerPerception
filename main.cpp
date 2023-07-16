#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include "SimpleNeuralNetWork.hpp"
#include "Topology.hpp"
#include <eigen3/Eigen/Eigen>
#include <bits/stdc++.h>


typedef std::vector<std::vector<Matrix>> Data ;




int main(){
    
    std::srand(time(0));
    Eigen::MatrixXf MatrixF ;
    
    
    std::vector< std::vector<std::vector<double>>> in{{{0,0}}, 
                                                      {{0,1}}, 
                                                      {{1,0}}, 
                                                      {{1,1}}};

    std::vector< std::vector<std::vector<double>>> out{{{0}},
                                                       {{1}},
                                                       {{1}},
                                                       {{0}}};

    Data inputs ;
    Data outputs ;
      
    for(int i = 0; i < in.size(); i++) {

        std::vector<Matrix> i1;
        i1.push_back(Matrix(in[i]));
        std::vector<Matrix> i2;
        i2.push_back(i1[0]) ;
        inputs.push_back(i2);


        std::vector<Matrix> o1;
        o1.push_back(Matrix(out[i]));
        std::vector<Matrix> o2;
        o2.push_back(o1[0]) ;
        outputs.push_back(o2);

        //delete i1[0];
        //delete o1[0];

    }


    
  

    
    Topology topology({2,3, 1});
    //Layers layerss ;
    //layerss = topology.layers ;
    Layers layer = topology.layers;
    int data_size = 4 ;
    int epochs =  100;
    
    
    SimpleNeuralNetWork model(topology) ;
    /*
    Neurons c(1,1);
    Neurons d(1,1);
    Matrix e(1,1, false) ;
    
    std::cout << c.neurons.element[0][0] << "\n" ;
    c.neurons = c.neurons.sigmoid()*2 + c.neurons ;
    std::cout << c.neurons.element[0][0] << "\n" ;
    
    std::cout << c.neurons.element[0][0] << "\n" ;
    c.neurons = model.sigmoid(c.neurons)*3 + c.neurons;
    d.neurons = d.neurons*3 + 1;
    e = d.neurons +1 ;
    std::cout << c.neurons.element[0][0] << "\n" ;
    std::cout << d.neurons.element[0][0] << "\n" ;
    std::cout << e.element[0][0] << "\n" ;
    d.neurons = c.neurons+1 ;
    std::cout << c.neurons.element[0][0] << "\n" ;
    std::cout << d.neurons.element[0][0] << "\n" ;
    std::cout << e.element[0][0] << "\n" ;
   
    exit(0);
     */
    std::vector<std::vector<Matrix>> data ;

    
    //model.train(inputs[0], outputs[0]);
    for(int j = 0; j < epochs ; j++){

        for(int i =0; i < data_size; i++){
             model.train(inputs[i], outputs[i]);

    }

    }
    
    std::vector<Matrix> _weights = model.weights;
    std::vector<Matrix> _bias = model.bias;
    
    

    for(int i =0; i < data_size; i++){
        model.prediction(inputs[i], outputs[i]);

    }

    
    std::cout << "\n" ;
    for(int n = 0; n < _weights.size(); n++){
        std::cout << "W :\n";
        for(int i = 0; i < _weights[n].rows(); i++){
            for(int j = 0; j < _weights[n].columns(); j++){

                std::cout << _weights[n].at(i,j) << " ";
              
            }
            std::cout << "\n";

        }

        std::cout << "\n";

    }

    std::cout << "\n" ;
    for(int n = 0; n < _bias.size(); n++){
        std::cout << "b :\n";
        for(int i = 0; i < _bias[n].rows(); i++){
            for(int j = 0; j < _bias[n].columns(); j++){

                std::cout << _bias[n].at(i,j) << " ";
              
            }
            std::cout << "\n";

        }
        std::cout << "\n";


    }
    

     
    return 0 ;
} 

