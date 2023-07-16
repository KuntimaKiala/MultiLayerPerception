#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#pragma once

#include <vector>
#include "Neurons.hpp"
#include <cassert>

typedef std::vector<Neurons> Layers ;
class Topology
{
public:
    Topology();
    Topology(int numLayer, std::vector<int> NPL);
    Topology(std::vector<int> NPL);
    ~Topology();
    Layers layers;
    std::vector<int> NeuronsPerLayer ;
    

private:

};

#endif