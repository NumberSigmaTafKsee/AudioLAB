#include "electric-frankenstein.hpp"
#include "Plot.hpp"
#include <complex>
#include "ml_helpers.hpp"

using namespace FX::Distortion;
Plot<floatType> plt;

Lua::LuaJIT luajit("chebystein.lua")

void plot(Matrix & m, const std::string & caption = "")
{    
    std::vector<floatType> r(m.rows());
    for(size_t i = 0; i < m.rows(); i++) r[i] = m(i,0);
    plt.plot_x(r.data(),r.size(),caption.c_str());
}



void ChebyWet()
{
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;    
    for(int i = 0; i < 128; i++) {
        auto sample = std::sin(2*M_PI*(double)i/128.0);
        examples.push_back(sample);
        training.push_back(udo1(sample,1));
    }
    plt.plot_x(examples.data(),examples.size(),"examples");
    plt.plot_x(training.data(),training.size(),"training");
    
    Matrix e = matrix_new(examples.size(),1,examples);
    Matrix t = matrix_new(training.size(),1,training);
        
    NeuralNetwork net;
    ActivationType atype = TANH;
    constexpr size_t CHEBYS=32;
    Layer * input = new Layer(INPUT,1,LINEAR);        
    Layer * hidden1= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden2= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden3= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * output= new Layer(OUTPUT,1,atype);
    net.addLayer(input);        
    net.addLayer(hidden1);    
    net.addLayer(hidden2);    
    net.addLayer(hidden3);    
    net.addLayer(output);
    net.connect();

    net.loss_widget = 1e-6;

    ParameterSet p(e,t,1000,8);
    p.batch_size=4;
    p.learning_rate = 0.1;
    p.momentum_factor = 0.9;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = RMSPROP_OPTIMIZER;    
    net.LastConnection()->bias.setRandom();
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,NONSTOCHASTIC);
    net.ForwardPass(e);
    Matrix x;
    net.GetOutput(x);
    plot(x);
}


int main(int argc, char * argv[]) {     	
    plt.setstyle("lines");
    ChebyWet();
    sleep(15);
}
