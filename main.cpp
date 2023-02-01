#include<iostream>
#include "neat.h"

using namespace neat;

int main(){

    cGene A(0, 3), B(1, 3), C(2, 3), D(1, 4), E(4, 3), F(0, 4), G(3, 4);

    A.weight = 0.7;
    B.weight = -0.5; B.enable = false;
    C.weight = 0.5;
    D.weight = 0.2;
    E.weight = 0.4;
    F.weight = 0.6;
    G.weight = 0.6;

    Genotype g;

    g << A << B << C << D << E << F << G;

    std::vector<int> in = {0, 1, 2}; std::vector<int> out = {3};

    g.setSensors(in); g.setOutputs(out);

    Phenotype p(g);

    std::vector<double> input = {1, 0, -1};

    std::vector<double> output;

    for(int i=0; i<4; ++i){
        output = p.eval(input);
        std:: cout << output[0] << '\n';
    }

    g.mutateAddConnection();

    Phenotype q(g);

    std::cout << '\n' << '\n';

    for(int i=0; i<4; ++i){
        output = q.eval(input);
        std:: cout << output[0] << '\n';
    }

    return 0;
}
