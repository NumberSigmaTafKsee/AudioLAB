#include "electric-frankenstein.hpp"
#include "Plot.hpp"
#include "ml_helpers.hpp"

#include <algorithm>
#include <functional>

#include "Analog/MLAnalogSVF.hpp"
#include "Analog/VABlitOscillators.hpp"

#include "LuaJIT.hpp"

Lua::LuaJIT luajit("chebystein.lua");
auto batch_size = Lua::LuaVariable(luajit,"batch_size");
auto learning_rate = Lua::LuaVariable(luajit,"learning_rate");
auto momentum_factor = Lua::LuaVariable(luajit,"momentum_factor");
auto regularization_strength = Lua::LuaVariable(luajit,"regularization_strength");
auto search_time = Lua::LuaVariable(luajit,"search_time");	
auto verbose = Lua::LuaVariable(luajit,"verbose");
auto shuffle = Lua::LuaVariable(luajit,"shuffle");
auto optimizer = Lua::LuaVariable(luajit,"optimizer");
auto stochastic = Lua::LuaVariable(luajit,"stochastic");
auto iters = Lua::LuaVariable(luajit,"iters");
auto cutoff = Lua::LuaVariable(luajit,"cutoff");
auto q = Lua::LuaVariable(luajit,"q");


#define MATRIXRC(w) std::cout << w.rows() << "," << w.cols() << std::endl
#define PRINT(X) std::cout << (X) << std::endl
#define LOOP(x,start,end) for(size_t x = start; x < end; x++)
#define PCQ(c,q) std::cout << c << "," << q << std::endl
std::vector<floatType> examples;
std::vector<floatType> training;
    

typedef floatType DspFloatType;
Plot<floatType> plt;

void clamp(Matrix & m, floatType min=0, floatType max=1.0) {
    for(size_t i = 0; i < m.rows(); i++)
    for(size_t j = 0; j < m.cols(); j++)
    {
        if(m(i,j) < min) m(i,j) = min;
        if(m(i,j) > max) m(i,j) = max;    
    }
}


void plot(Matrix & m, const std::string & caption = "")
{    
    std::vector<floatType> r(m.cols());
    for(size_t i = 0; i < m.cols(); i++) r[i] = m(0,i);
    plt.plot_x(r.data(),r.size(),caption.c_str());
}


void SVF(ActivationType atype) {
   
    Analog::Oscillators::Blit::BlitSaw osc;
    Analog::Oscillators::Blit::BlitSquare osc2;
	Analog::Filters::AnalogSVF::AnalogSVF filter(44100,cutoff.asDouble(),q.asDouble());
	
    for(int i = 0; i < 256; i++) {
        double s = osc2.Tick();
        double x = filter.Tick(s);
        examples.push_back(s);
        training.push_back(x);
    }
    
    //RemoveDCNorm(training);
    plt.plot_x(examples.data(),examples.size(),"examples");
    plt.plot_x(training.data(),training.size(),"training");
	
    Matrix e = matrix_new(1,256,examples);
    Matrix t = matrix_new(1,256,training);

    NeuralNetwork net;
    net.load("svf.net");
    
    /*    
    Layer * input 	 = new Layer(INPUT,256,LINEAR);            
    Layer * hidden1	 = new Layer(HIDDEN,2,LINEAR);            
    Layer * hidden2	 = new Layer(HIDDEN,2,RELU);            
    Layer * synth	 = new Layer(OUTPUT,256,TANH);    
    
    
    net.addLayer(input);            
    net.addLayer(hidden1);    
    net.addLayer(hidden2); 
    net.addLayer(synth);
    net.connect();
    net.loss_widget = 1e-6;
    
    ParameterSet p(e,t,iters.asDouble(),1);
    p.batch_size=batch_size.asDouble();
    p.learning_rate = learning_rate.asDouble();
    p.momentum_factor = momentum_factor.asDouble();
    p.regularization_strength = regularization_strength.asDouble();
    p.search_time=search_time.asDouble();
    p.verbose = verbose.asBool();
    p.shuffle = shuffle.asDouble();
    p.optimizer = optimizer.asDouble();    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,stochastic.asDouble());
    */
      
    auto m = net.LastConnection()->weights;
    net.ForwardPass(e);
    std::cout << net.layers[2]->input << std::endl;
    Matrix out(1,256);
    net.GetOutput(out);    
    plt.plot_x(out.data(),256,"neural");
    net.save("svf.net");
       
}

int main(int argc, char * argv[]) {     
    plt.setstyle("lines");   
    SVF(AUTO);   
    sleep(15);
}


