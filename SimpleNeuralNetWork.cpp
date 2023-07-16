#include "SimpleNeuralNetWork.hpp"
#include <cmath>
SimpleNeuralNetWork::SimpleNeuralNetWork()
{

}

SimpleNeuralNetWork::~SimpleNeuralNetWork()
{

  /*
   No need of deleting _topology because topology is deleted in topology.cpp
  */

   

}

SimpleNeuralNetWork::SimpleNeuralNetWork(Topology& topology) :_topology(&topology){
    
    //_topology.NeuronsPerLayer.size() -1 :_topology(topology)
   for(int i = 0; i<(*_topology).NeuronsPerLayer.size() -1; i++){
        
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
        
    
        if (l == 5) {
            std::cout << "W0 :" << "\n" ;  

            std::cout << weights[l].at(0,0) <<" "<< weights[l].at(0,1) <<" "<< weights[l].at(0,2) << "\n";
            std::cout << weights[l].at(1,0) <<" "<< weights[l].at(1,1) <<" "<< weights[l].at(1,2)<< "\n";;
            

        }

        if (l == 0) {
            std::cout << "W1 :" << "\n" ;  

            std::cout << weights[l].at(0,0) << "\n";
            std::cout << weights[l].at(1,0) << "\n";
            //std::cout << weights[l].at(2,0) << "\n";
            

        }
        
        
        layers[l+1].neurons = layers[l].neurons*(weights[l]) + bias[l];
        //layers[l+1].neurons = layers[l+1].neurons + bias[l];
        layers[l+1].neurons = layers[l+1].neurons.sigmoid() ;
      
       
   
    }

  std::cout << "a_2 :"<<layers.back().neurons.at(0,0) << "\n";
  //std::function<int(const int&)>
  //layers.back().neurons = layers.back().neurons.ApplyFunction(std::function <double(const double & m)>{predit});
  //layers.back().neurons = layers.back().neurons.ApplyFunction(predit);
    
   std::cout << "predicion :" << layers.back().neurons.at(0,0) << "\n";
  

  
}

    




void SimpleNeuralNetWork::backPropagation(){
    calcErrors(Y);

    updateWeights();
}



void SimpleNeuralNetWork::updateWeights() { 


int nbLayers  = weights.size()  ;


    for(uint32_t l = 0; l < nbLayers; l++){
        Matrix A = layers[l].neurons.transpose() ;
        Matrix G = errors[l] ;
        Matrix learningrate_times_A_transposed_times_dC_over_dW = A * G *learningRate ;
        Matrix learningrate_times_dC_over_dW =  G*learningRate ;
   
        weights[l]  = weights[l] - learningrate_times_A_transposed_times_dC_over_dW;
        bias[l]     = bias[l]    - learningrate_times_dC_over_dW ;
        if (l == 0) {
            std::cout << "W1_up :" << "\n" ;  

            std::cout << weights[l].at(0,0) <<" "<< learningrate_times_A_transposed_times_dC_over_dW.at(0,0) <<" " << errors[l].at(0,0) << "\n";
            std::cout << weights[l].at(1,0) << " "<< learningrate_times_A_transposed_times_dC_over_dW.at(1,0)<<" " << errors[l].at(0,0)<< "\n";
            //std::cout << weights[l].at(2,0) << "\n";
            

        }
        //layers[l].neurons = layers[l].neurons.transpose() ;
    }


}


void SimpleNeuralNetWork::calcErrors(Matrix& Y)
{


    //std::cout << "\n\nBACKPROP\n" <<std::endl;

    std::cout << "y :" << Y.element[0][0] << std::endl;




int nbLayers  = layers.size() -1 ;
Matrix Y_Hat(layers.back().neurons.rows(), layers.back().neurons.columns(), false); 
Y_Hat = layers.back().neurons ;

errors.back() =   Y_Hat - Y;
//std::cout << nbLayers << std::endl;
for (uint32_t L = nbLayers-1; L > 0; L--){
    //(weights[L]).shape() ;
    Matrix Weight_Transposed =  (weights[L]).transpose() ;
    Matrix deriv = layers[L].neurons.sigmoid_deriv(); 
    errors[L-1] =  ((errors[L])*Weight_Transposed).Hadamard(deriv); 
 
    }
    
 
    
  
}

void SimpleNeuralNetWork::train(std::vector<Matrix>& input_data, std::vector<Matrix>& output_data)
{
    

    for (int nbInput = 0; nbInput < input_data.size(); nbInput++) {
        std::cout << std::endl ;
        X  = input_data[nbInput];
        //inputs.shape();
        std::cout <<"input"<< "("<<nbInput<<")"<<" :"<<  X.at(0,0)<<" " << X.at(0,1)<< std::endl ;
 
        Y = output_data[nbInput] ;
        std::cout <<"output"<< "("<<nbInput<<")"<<" :"<< Y.at(0,0)<<std::endl ;
        //outputs.shape();
        
        
        FeedForward();
        //exit(0);
        backPropagation();
        
    }
    
}



void SimpleNeuralNetWork::prediction(std::vector<Matrix>& input_data, std::vector<Matrix>& output_data){
    int nbLayer = layers.size() ;
    layers[0].neurons = X ;


    for (int nbInput = 0; nbInput < input_data.size(); nbInput++) {
        std::cout << std::endl ;
        X  = input_data[nbInput];
        //inputs.shape();
        std::cout <<"input"<< "("<<nbInput<<")"<<" :"<<  X.at(0,0)<<" " << X.at(0,1)<< std::endl ;
 
        Y = output_data[nbInput] ;
        std::cout <<"output"<< "("<<nbInput<<")"<<" :"<< Y.at(0,0)<<std::endl ;

        for(int l = 0; l < nbLayer -1 ; l++) {
        
            layers[l+1].neurons = layers[l].neurons*(weights[l]) + bias[l];
            //layers[l+1].neurons = layers[l+1].neurons + bias[l];
            layers[l+1].neurons = layers[l+1].neurons.sigmoid() ;

            if (l == 5) {
                std::cout << "W0 :" << "\n" ;  

                std::cout << weights[l].at(0,0) <<" "<< weights[l].at(0,1) <<" "<< weights[l].at(0,2) << "\n";
                std::cout << weights[l].at(1,0) <<" "<< weights[l].at(1,1) <<" "<< weights[l].at(1,2)<< "\n";;
                std::cout << std::endl;

            }
            
        
        
   
       }

       layers.back().neurons = layers.back().neurons.ApplyFunction(predit);
    
        std::cout << "predicion :" << layers.back().neurons.at(0,0) << "\n";
        


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

    return m1;

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


