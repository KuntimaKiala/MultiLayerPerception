#include "Topology.hpp"

Topology::Topology()
{

}



Topology::Topology(std::vector<int> NPL) : NeuronsPerLayer(NPL){

int numLayer = NeuronsPerLayer.size();

for(int l = 0; l <numLayer; l++) {
    layers.push_back(Neurons(1,NeuronsPerLayer[l]));
}


}

Topology::Topology(int numLayer, std::vector<int> NPL) : NeuronsPerLayer(NPL){

assert(numLayer == NeuronsPerLayer.size());

for(int l = 0; l <numLayer; l++) {
    layers.push_back(Neurons(1,NeuronsPerLayer[l]));
 
    
}

}



Topology::~Topology(){

 
}