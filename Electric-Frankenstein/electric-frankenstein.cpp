#include "electric-frankenstein.hpp"
#include "Plot.hpp"
#include <complex>

#include "SoundObject.hpp"
#include "FX/Amplifiers.hpp"
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

using namespace FX::Distortion;

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
    std::vector<floatType> r(m.rows());
    for(size_t i = 0; i < m.rows(); i++) r[i] = m(i,0);
    plt.plot_x(r.data(),r.size(),caption.c_str());
}


floatType chebyshev(floatType x, floatType A[], int order)
{	
    // To = 1
    // T1 = x
    // Tn = 2.x.Tn-1 - Tn-1
    // out = sum(Ai*Ti(x)) , i C {1,..,order}
    floatType Tn_2 = 1.0f;
    floatType Tn_1 = x;
    floatType Tn;
    floatType out  = A[0]*x;
    
    for(int n=2;n<=order;n++)
    {
        Tn    =   (2.0*x*Tn_1 - Tn_2);
        out   +=  A[n-1]*Tn;
        Tn_2  =   Tn_1;
        Tn_1  =   Tn;
    }
    return out;
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
    constexpr size_t CHEBYS=64;
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
	
    ParameterSet p(e,t,iters.asDouble(),8);
    p.batch_size=batch_size.asDouble();
    p.learning_rate = learning_rate.asDouble();
    p.momentum_factor = momentum_factor.asDouble();
    p.regularization_strength = regularization_strength.asDouble();
    p.search_time=search_time.asDouble();
    p.verbose = verbose.asBool();
    p.shuffle = shuffle.asDouble();
    p.optimizer = RMSPROP_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,stochastic.asDouble());
    net.ForwardPass(e);
    Matrix x;
    net.GetOutput(x);
    plot(x);
}




int main(int argc, char * argv[]) {     
	
    plt.setstyle("lines");
    //Cheby(RELU);
    ChebyWet();
    sleep(15);
}
