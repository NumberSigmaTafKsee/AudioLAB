#include "electric-frankenstein.hpp"
#include "Plot.hpp"
#include "ml_helpers.hpp"

#include <algorithm>
#include <functional>

#define PRINT(X) std::cout << (X) << std::endl
#define LOOP(start,end) for(size_t i = start; i < end; i++)
#define MATRIXRC(w) std::cout << w.rows() << "," << w.cols() << std::endl

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

const int cheby_order = 4;

double quadratic(double s, double a, double b, double c) {
    return a*s*s + b*s + c;
}
double chebyshev(double s, Matrix & w) {
	floatType A[cheby_order];
	for(size_t i = 0; i < cheby_order; i++) A[i] = w(0,i);
	return chebyshev(s,A,cheby_order);
}
double equation(double s, Matrix & w) {		
	return -chebyshev(s,w);
}
struct QuadraticLayer : public Layer
{
    NeuralNetwork * net;
    // calculate the output from the weights
    QuadraticLayer(NeuralNetwork * neuralnet,size_t s,ActivationType a,bool useAutoDiff = true) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet) {

    }    
    
    Matrix Loss(Matrix & output, Matrix & target)
    {
		Matrix w = net->layers[net->layers.size()-2]->input;				
		Matrix input = net->layers[0]->input;
		Matrix m(1,64);
		for(size_t s = 0; s < 64; s++)
		{						
			floatType in = input(0,s);        
			floatType y = equation(in,w) - equation(0,w);
			m(0,s) = y;
		}		
		std::vector<floatType> out(64);
		for(size_t i = 0; i < 64; i++) out[i] = m(0,i);
        RemoveDCNorm(out);
        for(size_t i = 0; i < 64; i++) m(0,i) = out[i];
		return m - target;
	} 
	   	
    void Formula(Matrix & e) {        
        std::vector<floatType> out(64);
        net->ForwardPass(e);        
        output.resize(1,64);
        
        Matrix w = net->layers[net->layers.size()-2]->input;		
		Matrix input = net->layers[0]->input;
		Matrix m(1,64);
		std::cout << w << std::endl;
		for(size_t s = 0; s < 64; s++)
		{						
			floatType in = input(0,s);        
			floatType y = equation(in,w) - equation(0,w);
			m(0,s) = y;
		}		
        for(size_t i = 0; i < 64; i++) out[i] = output(0,i);
        RemoveDCNorm(out);
        plt.plot_x(out.data(),out.size(),"neural");
    }        
    void Activate(Matrix &m)
    {			
		getOutput(m);			
    }
    void Grad(Matrix &m)
    {		
		using namespace boost::math::differentiation;
        constexpr unsigned Order = 2;           
                
        Matrix t = m.eval();
        Matrix& w = net->connections[net->connections.size()-1]->weights;							
		Matrix  x(cheby_order,1);
		std::cout << t << std::endl;
		for(size_t i = 0; i < w.rows(); i++)
		{
			floatType c = 0;
			for(size_t j = 0; j < w.cols(); j++) {
				c += w(i,j);
			}
			c/=w.cols();
			x(i,0) = c;
		}
		//std::cout << x << std::endl;		
		for(size_t i = 0; i < t.rows(); i++)
		{						
			for(size_t j = 0; j < t.cols(); j++) {
				floatType in = t(i,j);
				floatType y = equation(in,x) - equation(0,x);
				t(i,j) = y;
			}
		}
		
        for(size_t i = 0; i < t.rows(); i++)
        {
			for(size_t j = 0; j < t.cols(); j++) {
				auto const x = make_fvar<double, Order>(t(i,j));  
				m(i,j) = x.derivative(1);    
				if(std::isnan(m(i,j))) m(i,j) = 0;
				if(std::isinf(m(i,j))) m(i,j) = x < 0? -1:1;            
			}
		}
    } 
    void getOutput(Matrix & m)
    {		
		Matrix w = net->layers[net->layers.size()-2]->input;				
		Matrix input = net->layers[0]->input;				
		for(size_t s = 0; s < 64; s++)
		{						
			floatType in = input(0,s);        
			floatType y = equation(in,w) - equation(0,w);
			m(0,s) = y;
		}		
		std::vector<floatType> out(64);
		for(size_t i = 0; i < 64; i++) out[i] = m(0,i);
        RemoveDCNorm(out);
        for(size_t i = 0; i < 64; i++) m(0,i) = out[i];       		
	}
	
    
};

void Quadratic(ActivationType atype) {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    
    for(int i = 0; i < 64; i++) {
        double s = std::sin(2*M_PI*i/64.0);
        double x = FX::Distortion::Fold(2*s);
        examples.push_back(s);
        training.push_back(x);
    }
    //RemoveDCNorm(training);
    plt.plot_x(examples.data(),examples.size(),"examples");
    plt.plot_x(training.data(),training.size(),"training");

    Matrix e = matrix_new(1,64,examples);
    Matrix t = matrix_new(1,64,training);

    NeuralNetwork net;
    const int COEFFS=cheby_order;
    Layer * input = new Layer(INPUT,64,LINEAR);
    Layer * h1= new Layer(HIDDEN,64,TANH);    
    Layer * h2= new Layer(HIDDEN,64,TANH);
    Layer * h3= new Layer(HIDDEN,64,TANH);    
    Layer * h4= new Layer(HIDDEN,COEFFS,TANH);
    Layer * h5= new Layer(HIDDEN,COEFFS,TANH);
    Layer * h6= new Layer(HIDDEN,COEFFS,TANH);
    Layer * h7= new Layer(HIDDEN,COEFFS,TANH);
    Layer * h8= new Layer(HIDDEN,COEFFS,TANH);
    Layer * h9= new Layer(HIDDEN,COEFFS,TANH);
    Layer * h10= new Layer(HIDDEN,COEFFS,TANH);	
    QuadraticLayer * output= new QuadraticLayer(&net,64,AUTO);
    net.addLayer(input);    
    net.addLayer(h1);          
    net.addLayer(h2);
    net.addLayer(h3);       
    /*
    net.addLayer(h4);           
    net.addLayer(h5);     
    net.addLayer(h6);    
    net.addLayer(h7);
    net.addLayer(h8);              
    net.addLayer(h9);                  
    */
    net.addLayer(h10);       
    net.addLayer(output);
    net.connect();
    net.loss_widget = 1e-5;
    
    ParameterSet p(e,t,1000,64);
    p.batch_size=1;
    p.learning_rate = 0.001;
    p.momentum_factor = 0.9;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = RMSPROP_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,NONSTOCHASTIC);
    
    
    auto m = net.connections[net.connections.size()-1]->weights;
    output->Formula(e);

}

int main(int argc, char * argv[]) {     
    plt.setstyle("lines");   
    Quadratic(AUTO);   
    sleep(15);
}

