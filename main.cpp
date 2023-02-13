#include<iostream>
#include "neat.h"

using namespace neat;

int main(){
    ////////////////// PARENT 1 ////////////////////
    cGene a1(1,4), a2(2,4), a3(3,4), a4(2,5), a5(5,4), a8(1,5);

    a1.setInov(1);
    a2.setInov(2).setEnable(false);
    a3.setInov(3);
    a4.setInov(4);
    a5.setInov(5);
    a8.setInov(8);

    Genotype parent1;

    parent1 << a1 << a2 << a3 << a4 << a5 << a8;

    std::vector<int> ins = {1,2,3}, outs = {4};

    parent1.setSensors(ins).setOutputs(outs);

    ////////////////// PARENT 2 //////////////////////

    cGene b5(5,4), b6(5,6), b7(6,4), b9(3,5), b10(1,6);

    b5.setInov(5).setEnable(false);
    b6.setInov(6);
    b7.setInov(7);
    b9.setInov(9);
    b10.setInov(10);

    Genotype parent2;

    parent2 << a1 << a2 << a3 << a4 << b5 << b6 << b7 << b9 << b10;

    parent2.setSensors(ins).setOutputs(outs);

    /////////////////// CROSSOVER ////////////////

    System s(1);

    s.mutate(parent1);
    s.mutate(parent2);

    Genotype child = s.crossover(parent1, parent2, 1, 1);

    std::cout << s.distance(parent1, parent2) << '\n';

    parent1.mutateAddConnection().setInov(11);

    std::cout << s.distance(parent1, parent2) << '\n';

    return 0;
}
