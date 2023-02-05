#include<iostream>
#include "neat.h"

using namespace neat;

int main(){
    cGene A(0,3), B(1, 3), C(2, 3), D(3, 4);

    std::vector<int> ins = {0,1,2}, outs = {4};

    Genotype g; g << A << B << C << D;
    g.setSensors(ins).setOutputs(outs);

    System s(1); s.perturbWeight = 0.2;

    s.mutate(g);

    return 0;
}
