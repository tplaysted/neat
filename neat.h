#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

namespace neat{

const char SENSOR = 0x01;
const char HIDDEN = 0x02;
const char OUTPUT = 0x03;

class cGene;
class Genotype;
class Node;
class Phenotype;
class System;
class nGene;

double sigmoid(double x){
    return 1 / (1 + exp(-4.9 * x));
}

class cGene{
public:

    int inov=-1;
    int inNode=-1;
    int outNode=-1;

    double weight=1;

    bool enable=true;

    cGene()= default;

    cGene(int in, int out){
        this->inNode = in;
        this->outNode = out;
    }

};

class nGene{
public:

    int index=-1;
    char type=HIDDEN;

    nGene()= default;
};

class Genotype{
    public: 

    std::vector<cGene> cGenes; // contiguous array
    std::vector<nGene> nGenes; // indexed array (node k goes in index k)

    std::vector<int> sensors;
    std::vector<int> outputs;

    Genotype &operator << (cGene &g){ // overload << operator for easy attaching of new connect genes
        if(this->nGenes[g.outNode].type == SENSOR){
            throw std::invalid_argument("Connection cannot output to sensor node");
        }

        this->cGenes.push_back(g);
        int in = g.inNode;
        int out = g.outNode;

        if(in >= this->nGenes.size() || out >= this->nGenes.size()){ // resizing node array if needed
            this->nGenes.resize(std::max(in, out) + 1);
        }

        if(this->nGenes[in].index == -1){ // create node if unassigned
            this->nGenes[in].index = in;
        }

        if(this->nGenes[out].index == -1){ // create node if unassigned
            this->nGenes[out].index = out;
        }

        return *this;
    }

    void setSensors(std::vector<int> &s){ // for setting sensor nodes
        for(auto & i : s){
            if(i >= this->nGenes.size()){
                throw std::invalid_argument("Sensor index does not correspond to any node");
            }
            this->nGenes[i].type = SENSOR;
        }
        this->sensors = s;
    }

    void setOutputs(std::vector<int> &s){ // for setting output nodes
        for(auto & i : s){
            if(i >= this->nGenes.size()){
                throw std::invalid_argument("Output index does not correspond to any node");
            }
            this->nGenes[i].type = OUTPUT;
        }
        this->outputs = s;
    }

    cGene& mutateAddConnection(){
        std::vector<int> ins;
        for(auto & node : this->nGenes){ // build list of potential input nodes to choose from
            if(node.index != -1){
                ins.push_back(node.index);
            }
        }

        std::random_device rd;
        std::mt19937 g(rd()); // seed random generator

        std::shuffle(ins.begin(), ins.end(), g); // reorder
        int in = ins[0]; // input node

        std::vector<int> outs;
        for(auto & node : this->nGenes){ // build list of potential input nodes to choose from
            if(node.index != -1 && node.type != SENSOR && node.index != in){ // we don't allow self connections
                outs.push_back(node.index);
            }
        }

        std::shuffle(outs.begin(), outs.end(), g); // reorder
        int out = outs[0]; // output node

        cGene newConnection(in, out); // new connection has weight 1 for now

        *(this) << newConnection;

        return this->cGenes.back(); // return a reference to the newly created gene
    }

    void mutateAddNode(int newNodeVal){
        for(auto & node : this->nGenes){
            if(node.index == newNodeVal){
                throw std::invalid_argument("Node index already exists in genome");
            }
        }

        std::vector<int> temp; temp.resize(this->cGenes.size());

        for(int i=0; i<temp.size(); ++i){ // vector of indexes corresponding to connect genes
            temp[i] = i;
        }

        std::random_device rd;
        std::mt19937 g(rd()); // seed random generator

        std::shuffle(temp.begin(), temp.end(), g); // reorder
        int index = temp[0]; // index of gene where node is to be inserted

        cGene &oldGene = this->cGenes[index]; // get a reference of the gene to be replaced

        cGene newGeneA(oldGene.inNode, newNodeVal); newGeneA.weight = 1;
        cGene newGeneB(newNodeVal, oldGene.outNode); newGeneB.weight = oldGene.weight;

        oldGene.enable = false;

        *this << newGeneA << newGeneB;
    }
};

class Node{
    public:

    double value = 0;

    char type = HIDDEN;

    std::vector<Node*> inputs;
    std::vector<Node*> outputs;

    std::vector<double> weights; // each weight corresponds to an input connection

    bool visited = false;

    double computePass(){ // recursive tree search to get node value
        this->visited = true;

        double out = 0;

        for(int i=0; i < this->inputs.size(); ++i){
            if(this->inputs[i]->type == SENSOR || this->inputs[i]->visited){ // stop conditions
                out += this->weights[i] * this->inputs[i]->value;
            } else{
                out += this->weights[i] * this->inputs[i]->computePass(); // recursive step
            }
        }

        out = sigmoid(out);

        this->value = out;
        return out;
    }

    void resetPass(){ // recursive tree search to reset in preparation for next compute pass
        this->visited = false;

        for(auto & input : this->inputs){
            if(input->type != SENSOR && input->visited){ // search conditions
                input->resetPass();
            }
        }
    }
};

class Phenotype{
    public:

    std::vector<int> sensors; // indices keep track of sensors and outputs
    std::vector<int> outputs;

    std::vector<Node*> nodes; // contiguous array of node pointers

    explicit Phenotype(Genotype& g){
        this->sensors = g.sensors;
        this->outputs = g.outputs;

        // first we must find an upper bound for the number of nodes in the network 

        int numNodes = 0;

        for(auto & genome : g.cGenes){
            if(genome.inNode > numNodes){
                numNodes = genome.inNode;
            }

            if(genome.outNode > numNodes){
                numNodes = genome.outNode;
            }
        }

        numNodes++;

        nodes.resize(numNodes); // resize node pointer container

        for(int & sensor : g.sensors){ // initialise sensor nodes
            if(sensor > numNodes){
                throw std::invalid_argument("Sensor index exceeds number of nodes");
            }

            this->nodes[sensor] = new Node;
            this->nodes[sensor]->type = SENSOR;
        }

        for(int & output : g.outputs){ // initialise output nodes
            if(output > numNodes){
                throw std::invalid_argument("Output index exceeds number of nodes");
            }

            this->nodes[output] = new Node;
            this->nodes[output]->type = OUTPUT;
        }

        // now we iterate over all genes, creating connections and initialising hidden nodes as we go

        for(auto & gene : g.cGenes){
            int in = gene.inNode;
            int out = gene.outNode;

            if(this->nodes[in] == nullptr){ // create new nodes if needed
                this->nodes[in] = new Node; 
            }

            if(this->nodes[out] == nullptr){
                this->nodes[out] = new Node; 
            }

            if(!(gene.enable)){ // connection must be enabled
                continue;
            }

            if(this->nodes[out]->type == SENSOR){ // no inputs to sensor nodes
                continue;
            }

            nodes[out]->inputs.push_back(nodes[in]); // assign connections and weights
            nodes[out]->weights.push_back(gene.weight);

            nodes[in]->outputs.push_back(nodes[out]);
        }
    }

    ~Phenotype(){
        for(auto & node : nodes){
            delete node;
        }
    }

    std::vector<double> eval(std::vector<double> &input){
        std::vector<double> out; // return vector

        if(input.size() != this->sensors.size()){ // argument checking
            throw std::invalid_argument("Input vector length does not match number of sensor nodes");
        }

        for(int i=0; i < input.size(); ++i){
            this->nodes[sensors[i]]->value = input[i]; // set sensor values
        }

        for(int i=0; i < this->outputs.size(); ++i){
            this->nodes[outputs[i]]->computePass(); // compute pass over all the output nodes
            out.push_back(this->nodes[outputs[i]]->value); // fill return vector 
        }

        for(int i=0; i < this->outputs.size(); ++i){
            this->nodes[outputs[i]]->resetPass(); // reset pass over all the output nodes
        }

        return out;
    }
};

class System{
public:
    std::vector<Genotype> population;
    std::vector<cGene> geneList;

    explicit System(int popLimit){
        this->population.resize(popLimit);
    }



};


}