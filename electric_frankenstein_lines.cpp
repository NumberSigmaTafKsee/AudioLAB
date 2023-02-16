struct LineLayer : public Layer
{
    NeuralNetwork * net;
    // calculate the output from the weights
    LineLayer(NeuralNetwork * neuralnet,size_t s,ActivationType a,bool useAutoDiff = true) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet) {

    }
    Matrix Loss(Matrix & m, Matrix & target) {        
        Matrix& w = net->connections[net->connections.size()-1]->weights;                
        //clamp(w);
        floatType y = w(0,0)*w(1,0)+w(2,0);
        Matrix x;
        net->GetOutput(x);
        x(0,0) = y;
        loss = x-target;
        return loss;        
    }
    floatType Formula(Matrix & w) {
        floatType y = w(0,0)*w(1,0)+w(2,0);
        return y;
    }
};

// Eigen::Matrix<Eigen::Matrix<floatType,Eigen::Dynamic,Eigen::Dynamic>,Eigen::Dynamic,Eigen::Dynamic> m(3,3);

void Line(ActivationType atype) {
    floatType y = 0.25;
    std::vector<floatType> examples = {y};
    std::vector<floatType> training = {y};
    
    Matrix e = matrix_new(1,1,examples);
    Matrix t = matrix_new(1,1,training);
    
    
    NeuralNetwork net;
    Layer * input = new Layer(INPUT,1,LINEAR);
    Layer * hidden= new Layer(HIDDEN,3,atype);
    LineLayer * output= new LineLayer(&net,1,LINEAR);
    net.addLayer(input);
    net.addLayer(hidden);
    net.addLayer(output);
    net.connect();

    ParameterSet p(e,t,1000,1);
    p.batch_size=1;
    p.learning_rate = 0.1;
    p.momentum_factor = 0.9;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = GD_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,STOCHASTIC);
    net.ForwardPass(e);
    
    auto m = net.connections[1]->weights;
    std::cout << m(0,0)*m(1,0)+m(2,0) << std::endl;
}


int main(int argc, char * argv[]) {     
    p.setstyle("lines");
    Line(TANH);   
    sleep(15);
}
