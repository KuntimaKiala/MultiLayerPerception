# Mutilayer Perception 
A C++ Multilayer Perception no libraries  <br />
I just wanted to create it from strach for c++ studying purposes<br />
XOR


## Usage :
g++ main.cpp SimpleNeuralNetWork.cpp Topology.cpp Neurons.cpp Matrix.cpp -o run  <br />
./run <br />
This command if for training


# NB : 
Since everything is created from scratch even the Matrix library the performance might not be the best


# TODO :
1- After training put the weigths in a file (egs : .txt) <br />
2- Make a predit based on the weigth from training (separate training from prediction) <br />
3- Create a class with all the activation functions <br />
4- User should be able of choose the the activation function from the command line
