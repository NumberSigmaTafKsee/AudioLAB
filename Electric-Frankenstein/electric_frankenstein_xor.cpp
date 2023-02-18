void XOR(ActivationType atype, floatType lt, floatType mf)
{
    std::vector<floatType> examples = {0,0,0,1,1,0,1,1};
    std::vector<floatType> training = {0,1,1,0};
    std::vector<floatType> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<floatType> training_bp = {-1,1,1,-1};

    Matrix e = matrix_new(4,2,examples);
    Matrix t = matrix_new(4,1,training);
        
    std::vector<int64_t> hidden = {16};
    std::vector<ActivationType> activations = {atype};
    Network net(2,hidden,activations,1,LINEAR);
    ParameterSet p(e,t,1000,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = GD_OPTIMIZER;
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.train(p,STOCHASTIC);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix &output = net.GetOutput();
    std::cout << output << std::endl;
}

void XOR2(ActivationType atype, floatType lt, floatType mf)
{
    std::vector<floatType> examples = {0,0,0,1,1,0,1,1};
    std::vector<floatType> training = {0,1,1,0};
    std::vector<floatType> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<floatType> training_bp = {-1,1,1,-1};

    Matrix e = matrix_new(4,2,examples);
    Matrix t = matrix_new(4,1,training);
        
    //std::vector<int64_t> hidden = {5};
    //std::vector<ActivationType> activations = {atype};
    NeuralNetwork net;
    Layer * input = new Layer(INPUT,2,LINEAR);
    Layer * hidden= new Layer(HIDDEN,10,atype);
    Layer * output= new Layer(OUTPUT,1,LINEAR);
    net.addLayer(input);
    net.addLayer(hidden);
    net.addLayer(output);
    net.connect();
    ParameterSet p(e,t,1000,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 0;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = GD_OPTIMIZER;
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.train(p,STOCHASTIC);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix &outs = net.GetOutput();
    std::cout << outs << std::endl;
}    

int main(int argc, char * argv[]) {     
   p.setstyle("lines");
   XOR2(TANH,0.1,0.9);   
   sleep(15);
}
