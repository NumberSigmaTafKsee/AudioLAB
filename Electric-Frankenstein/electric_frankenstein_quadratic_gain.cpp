struct QuadraticGainLayer : public Layer
{
    NeuralNetwork * net;
    Matrix output;

    QuadraticGainLayer(NeuralNetwork * neuralnet, size_t s,ActivationType a,bool useAutoDiff = true) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet){
        output.resize(1,s);
    }
    floatType equation(floatType in, Matrix& w)
    {
        return w(0,0)* in*in + w(1,0)*in + w(2,0);        
    }
    Matrix Loss(Matrix & m, Matrix & target) {        
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;              
        floatType in = input(0,0);                        
        floatType G = equation(in,w);
        auto x = net->GetOutput();
        x(0,0) = G*in;        
        loss = x-target;
        return loss;        
    }
    Matrix& getOutput() {
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;        
        floatType in = input(0,0);                
        floatType y = 0;
        floatType G = equation(in,w);
        output(0,0) = G;
        return output;
    }        
    void Formula(Matrix & e, Matrix & w) {        
        std::vector<floatType> out(e.rows());
        for(size_t i = 0; i < e.rows(); i++) {        
            floatType in = e(i,0);
            floatType G = equation(in,w);
            out[i] = G;
        }    
        p.plot_x(out.data(),out.size(),"neural");
    }    
    void PrintEquation(Matrix & w) {        
        std::cout << "a=" << w(0,0) << std::endl;
        std::cout << "b=" << w(1,0) << std::endl;
        std::cout << "c=" << w(2,0) << std::endl;
    }
};

void QuadraticGain(ActivationType atype) {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    
    for(size_t i = 0; i < 256; i++) {
        examples.push_back((double)i/256.0);
        training.push_back(2*(double)i/256.0);
    }
    p.plot_x(examples.data(),examples.size(),"examples");
    p.plot_x(training.data(),training.size(),"training");

    Matrix e = matrix_new(256,1,examples);
    Matrix t = matrix_new(256,1,training);

    NeuralNetwork net;
    Layer * input = new Layer(INPUT,1,LINEAR);
    Layer * hidden= new Layer(HIDDEN,3,atype);
    QuadraticGainLayer * output= new QuadraticGainLayer(&net,1,LINEAR);
    net.addLayer(input);
    net.addLayer(hidden);
    net.addLayer(output);
    net.connect();
    net.loss_widget = 1e-3;
    
    ParameterSet p(e,t,1000,256);
    p.batch_size=4;
    p.learning_rate = 0.01;
    p.momentum_factor = 0.1;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = ADAM_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,STOCHASTIC);
    net.ForwardPass(e);
    
    auto m = net.connections[1]->weights;
    output->Formula(e,m);
    output->PrintEquation(m);
}

int main(int argc, char * argv[]) {     
    p.setstyle("lines");
    QuadraticGain(TANH);      
    sleep(15);
}
