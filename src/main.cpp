#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include "SimpleNeuralNetWork.hpp"
#include "Topology.hpp"


typedef std::vector<std::vector<Matrix>> Data ;



int main(){
    
    std::srand(time(0));
  
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
      
    for(long unsigned int i = 0; i < in.size(); i++) {

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

    

    }


    
  

    
    Topology topology({2,3, 1});
    
    Layers layer = topology.layers;
    int data_size = 4 ;
    int epochs =  100;
    
    
    SimpleNeuralNetWork model(topology) ;
    
    std::vector<std::vector<Matrix>> data ;

    
    
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
    for(long unsigned int n = 0; n < _weights.size(); n++){
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
    for(long unsigned int n = 0; n < _bias.size(); n++){
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

