#include "SimpleNeuralNetWork.hpp"
#include <cmath>
SimpleNeuralNetWork::SimpleNeuralNetWork()
{

}

SimpleNeuralNetWork::~SimpleNeuralNetWork()
{

   
}

SimpleNeuralNetWork::SimpleNeuralNetWork(Topology& topology) :_topology(&topology){
    
    
   for(long unsigned int i = 0; i<(*_topology).NeuronsPerLayer.size() -1; i++){
        
        weights.push_back(Matrix(topology.NeuronsPerLayer[i], topology.NeuronsPerLayer[i+1],true)) ;
        bias.push_back(Matrix( 1,topology.NeuronsPerLayer[i+1],true)) ;
        errors.push_back(Matrix(1,topology.NeuronsPerLayer[i+1], false));


   }

    layers = (*_topology).layers ;

}


int SimpleNeuralNetWork::predit(const double &y_hat) {
    double output = 0.5;
    if (y_hat > output ) return 1;


    return 0 ;
}


void SimpleNeuralNetWork::FeedForward(){

    int nbLayer = layers.size() ;
    layers[0].neurons = X ;
   
    for(int l = 0; l < nbLayer -1 ; l++) {
              
        layers[l+1].neurons = layers[l].neurons*(weights[l]) + bias[l];
        layers[l+1].neurons = layers[l+1].neurons.sigmoid() ;
      
    }

  
}

    




void SimpleNeuralNetWork::backPropagation(){
    calcErrors(Y);
    updateWeights();
}



void SimpleNeuralNetWork::updateWeights() { 


int nbLayers  = weights.size()  ;


    for(int l = 0; l < nbLayers; l++){
        Matrix A = layers[l].neurons.transpose() ;
        Matrix G = errors[l] ;
        Matrix learningrate_times_A_transposed_times_dC_over_dW = A * G *learningRate ;
        Matrix learningrate_times_dC_over_dW =  G*learningRate ;
   
        weights[l]  = weights[l] - learningrate_times_A_transposed_times_dC_over_dW;
        bias[l]     = bias[l]    - learningrate_times_dC_over_dW ;
        
    }


}


void SimpleNeuralNetWork::calcErrors(Matrix& Y)
{

int nbLayers  = layers.size() -1 ;
Matrix Y_Hat(layers.back().neurons.rows(), layers.back().neurons.columns(), false); 
Y_Hat = layers.back().neurons ;

errors.back() =   Y_Hat - Y;

for (uint32_t L = nbLayers-1; L > 0; L--){
   
    Matrix Weight_Transposed =  (weights[L]).transpose() ;
    Matrix deriv = layers[L].neurons.sigmoid_deriv(); 
    errors[L-1] =  ((errors[L])*Weight_Transposed).Hadamard(deriv); 
 
    }
    
 
    
  
}

void SimpleNeuralNetWork::train(std::vector<Matrix>& input_data, std::vector<Matrix>& output_data)
{
    

    for (long unsigned int nbInput = 0; nbInput < input_data.size(); nbInput++) {

        X  = input_data[nbInput];
        Y = output_data[nbInput] ;
    
        FeedForward();
        backPropagation();
        
    }
    
}



void SimpleNeuralNetWork::prediction(std::vector<Matrix>& input_data, std::vector<Matrix>& output_data){
    int nbLayer = layers.size() ;
    layers[0].neurons = X ;


    for (int nbInput = 0; nbInput < input_data.size(); nbInput++) {
       
        X  = input_data[nbInput];
        
        std::cout <<"input"<< "("<<nbInput<<")"<<" :"<<  X.at(0,0)<<" " << X.at(0,1)<< std::endl ;
 
        Y = output_data[nbInput] ;
        std::cout <<"output"<< "("<<nbInput<<")"<<" :"<< Y.at(0,0)<<std::endl ;

        for(int l = 0; l < nbLayer -1 ; l++) {
        
            layers[l+1].neurons = layers[l].neurons*(weights[l]) + bias[l];
           
            layers[l+1].neurons = layers[l+1].neurons.sigmoid() ;

            
            

       }

       layers.back().neurons = layers.back().neurons.ApplyFunction(predit);
    
        std::cout << "predicion :" << layers.back().neurons.at(0,0) << "\n";
         std::cout << std::endl ;


    }
   
    

 
   
  

  
}


 Matrix SimpleNeuralNetWork::cost(Matrix& y, Matrix& y_hat)  {
    
    /*
    assert( y_hat.rows()       ==  y.rows()&& "different of numbers of columns cost");
    assert( y_hat.columns()    ==  y.columns()   && "different of numbers of rows cost");
    */
   
    
    
    Matrix output(y.rows(),y.columns(),false);
 
 
    
   /*
    for( unsigned i = 0; i < y.rows(); i++ ) { 
        for( unsigned j = 0; j < y.rows(); j++ ) {

            output.matrix[i][j] = y_hat.matrix[i][j] - y.matrix[i][j] ;
        }
    }
  */

    return output;
}




double SimpleNeuralNetWork::sigmoid (double& m1) {
    
    /*  Returns the value of the sigmoid function f(x) = 1/(1 + e^-x).
        Input: m1, a vector.
        Output: 1/(1 + e^-x) for every element of the input matrix m1.
    */
   
   
    double output ;
 
    output = (double)(1.0 / (1.0 + (double)std::exp(-m1)));
    return output;
}

double SimpleNeuralNetWork::sigmoid_deriv (double& m1){

    double output ;
    output = m1*(1-m1) ;

    return output;

}
double SimpleNeuralNetWork::ReLu(double & x){

   
 
    double output = 0.0;
    if (x <= output) return output;
    
    return x;  
    

}

double SimpleNeuralNetWork::ReLu_deriv (double& x) {

    double output = 0.0;
    if (x > output ) return 1.0;

    return 0.0;  
        
}


