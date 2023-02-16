#include "electric-frankenstein.hpp"
#include "Plot.hpp"
#include "ml_helpers.hpp"

#include <algorithm>
#include <functional>

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

const int cheby_order = 8;

double quadratic(double s, double a, double b, double c) {
    return a*s*s + b*s + c;
}
double chebyshev(double s, Matrix & w) {
	floatType A[cheby_order];
	for(size_t i = 0; i < cheby_order; i++) A[i] = w(i,0);
	return chebyshev(s,A,cheby_order);
}
double equation(double s, Matrix & w) {
	floatType foo=1.0;    
    for(size_t i = 0; i < w.rows(); i++) {		
		w(i,0) *= foo;
		foo   = sqrt(foo);		
	}
	return -chebyshev(s,w);
}
struct QuadraticLayer : public Layer
{
    NeuralNetwork * net;
    // calculate the output from the weights
    QuadraticLayer(NeuralNetwork * neuralnet,size_t s,ActivationType a,bool useAutoDiff = true) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet) {

    }
    
    void Activate(Matrix &m)
    {		
		getOutput(m);	
    }
    void Grad(Matrix &m)
    {		
		using namespace boost::math::differentiation;
        constexpr unsigned Order = 2;                  // Highest order derivative to be calculated.
        Matrix t = m;
        getOutput(t);
		auto const x = make_fvar<double, Order>((t(0,0)));  // Find derivatives at x=2.                            
		m(0,0) = x.derivative(1);    
		if(std::isnan(m(0,0))) m(0,0) = 0;
		if(std::isinf(m(0,0))) m(0,0) = x < 0? -1:1;            
    } 
    void getOutput(Matrix & m)
    {
		m(0,0) = 0;
		for(size_t i = 0; i < 8; i++) {
			Matrix& w = net->connections[net->connections.size()-i-1]->weights;
			Matrix& input = net->layers[0]->input;   
			floatType in = input(0,0);        
			floatType y = equation(in,w) - equation(0,w);
			m(0,0) += y;
		}
	}
	
    void Formula(Matrix & e, Matrix & w) {        
        std::vector<floatType> out(e.rows());
        for(size_t i = 0; i < e.rows(); i++) {        
            floatType in = e(i,0);
            floatType y = equation(in,w) - equation(0,w);
            out[i] = y;
        }    
        RemoveDCNorm(out);
        plt.plot_x(out.data(),out.size(),"neural");
    }        
};

void Quadratic(ActivationType atype) {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    
    for(size_t i = 0; i < 256; i++) {
        double s = std::cos(M_PI*2.0*(double)i/256.0);
        double x = FX::Distortion::Fold(2*s);
        examples.push_back(s);
        training.push_back(x);
    }
    //RemoveDCNorm(training);
    plt.plot_x(examples.data(),examples.size(),"examples");
    plt.plot_x(training.data(),training.size(),"training");

    Matrix e = matrix_new(256,1,examples);
    Matrix t = matrix_new(256,1,training);

    NeuralNetwork net;
    const int COEFFS=cheby_order;
    Layer * input = new Layer(INPUT,1,LINEAR);
    Layer * h1= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h2= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h3= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);    
    Layer * h4= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h5= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h6= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h7= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h8= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    /*
    Layer * h9= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    Layer * h10= new Layer(HIDDEN,COEFFS,FULLWAVE_RELU);
    */
    QuadraticLayer * output= new QuadraticLayer(&net,1,AUTO);
    net.addLayer(input);    
    net.addLayer(h1);    
    net.addLayer(h2);
    net.addLayer(h3);       
    net.addLayer(h4);   
    net.addLayer(h5);     
    net.addLayer(h6);    
    net.addLayer(h7);
    net.addLayer(h8);   
    /*    
    net.addLayer(h9);
    net.addLayer(h10);   
    */    
    net.addLayer(output);
    net.connect();
    net.loss_widget = 1e-5;
    
    ParameterSet p(e,t,10000,256);
    p.batch_size=8;
    p.learning_rate = 0.01;
    p.momentum_factor = 0.9;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = ADAM_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,STOCHASTIC);
    
    
    auto m = net.LastConnection()->weights;
    output->Formula(e,m);

}

int main(int argc, char * argv[]) {     
    plt.setstyle("lines");   
    Quadratic(AUTO);   
    sleep(15);
}

