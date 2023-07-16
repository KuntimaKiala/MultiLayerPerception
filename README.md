# Mutilayer Perception 
A C++ Multilayer Perception no libraries  <br />
I just wanted to do backpropagation from scratch using c++ for studying purposes<br />
XOR


## Usage :
g++ main.cpp SimpleNeuralNetWork.cpp Topology.cpp Neurons.cpp Matrix.cpp -o run <br />
./run <br />
This command is for training in the end there is a prediction


# NB : 
Since everything is created from scratch even the Matrix library the performance might not be the best


# TODO :
1- After training put the weigths in a file (egs : .txt) <br />
2- Make a predit based on the weigth from training (separate training from prediction) <br />
3- Create a class with all the activation functions <br />
4- User should be able of choose the the activation function from the command line
