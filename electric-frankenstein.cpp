
#include "electric-frankenstein.hpp"
typedef floatType DspFloatType;
#include "include/Undenormal.hpp"
#include "audio_iir_butterworth.hpp"
#include "Plot.hpp"

Plot<floatType> p;
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
    p.plot_x(r.data(),r.size(),caption.c_str());
}

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
        auto x = net->GetOutput();
        x(0,0) = y;
        loss = x-target;
        return loss;        
    }
    floatType Formula(Matrix & w) {
        floatType y = w(0,0)*w(1,0)+w(2,0);
        return y;
    }
};

struct QuadraticLayer : public Layer
{
    NeuralNetwork * net;
    // calculate the output from the weights
    QuadraticLayer(NeuralNetwork * neuralnet,size_t s,ActivationType a,bool useAutoDiff = true) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet) {

    }
    Matrix Loss(Matrix & m, Matrix & target) {        
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;        
        floatType in = input(0,0);        
        floatType y = w(0,0)*in*in+w(1,0)*in+w(2,0);
        auto x = net->GetOutput();
        x(0,0) = y;
        loss = x-target;
        return loss;        
    }
    void Formula(Matrix & e, Matrix & w) {        
        std::vector<floatType> out(e.rows());
        for(size_t i = 0; i < e.rows(); i++) {        
            floatType in = e(i,0);
            floatType y = w(0,0)*in*in+w(1,0)*in+w(2,0);
            out[i] = y;
        }    
        p.plot_x(out.data(),out.size(),"neural");
    }    
    void PrintEquation(Matrix & w) {        
        std::cout << "a=" << w(0,0) << std::endl;
        std::cout << "b=" << w(1,0) << std::endl;
        std::cout << "c=" << w(2,0) << std::endl;
    }
};

floatType chebyshev(floatType x, floatType A[], int order)
{
    // To = 1
    // T1 = x
    // Tn = 2.x.Tn-1 - Tn-2
    // out = sum(Ai*Ti(x)) , i C {1,..,order}
    floatType Tn_2 = 1.0f;
    floatType Tn_1 = x*A[0];
    floatType Tn;
    floatType out = Tn_1;

    for(int n=2;n<=order;n++)
    {
        Tn    =   2.0*x*Tn_1 - Tn_2;
        out   +=  Tn;
        Tn_2  =   Tn_1;
        Tn_1  =   A[n-1];
    }
    return out;
}

struct ChebyshevLayer : public Layer
{
    NeuralNetwork * net;
    Matrix output;
    int order;
    // calculate the output from the weights
    ChebyshevLayer(NeuralNetwork * neuralnet,int o, size_t s,ActivationType a,bool useAutoDiff = true) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet),order(o) {
        output.resize(1,1);
    }
    inline floatType cheby(int n, floatType x)
    {
        if(n == 0) return 1;
        if(n == 1) return x;
        return 2*x*cheby(n-1,x) - cheby(n-2,x);
    }
    virtual void Activate(Matrix& tmp) {            
        Matrix& w = net->connections[net->connections.size()-1]->weights;        
        Matrix& input = net->layers[0]->input;              
        floatType A[order];                 
        for(size_t i = 0; i < order; i++) {            
            A[i] = w(i,0);            
        }        
        for(size_t i = 0; i < tmp.rows(); i++)
            for(size_t j = 0; j < tmp.cols(); j++)
            {
                tmp(i,j) = std::clamp(chebyshev(input(0,0),A,order) - chebyshev(0,A,order),-1.0,1.0);
            }        
    }
    virtual void Grad(Matrix & tmp) {
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;              
        floatType A[order];
        Matrix m = tmp;        
        for(size_t i = 0; i < order; i++) {
            A[i] = w(i,0);            
        }        
        for(size_t i = 0; i < tmp.rows(); i++)
        for(size_t j = 0; j < tmp.cols(); j++)
        {
            m(i,j) = std::clamp(chebyshev(input(0,0),A,order) - chebyshev(0,A,order),-1.0,1.0);
        }        
        using namespace boost::math::differentiation;
        constexpr unsigned Order = 2;                  // Highest order derivative to be calculated.
                
        for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            auto const x = make_fvar<double, Order>((m(i,j)));  // Find derivatives at                                           
            m(i,j) = x.derivative(1);    
            if(std::isnan(m(i,j))) m(i,j) = 0;
            if(std::isinf(m(i,j))) m(i,j) = x < 0? -1:1;            
        }                
    }    
    
    Matrix Loss(Matrix & m, Matrix & target) {        
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;              
        floatType in = input(0,0);                
        floatType y = 0;
        floatType A[order];
        for(size_t i = 0; i < order; i++) A[i] = w(i,0);
        y = chebyshev(in,A,order);
        auto x = net->GetOutput();
        x(0,0) = y;        
        loss = x-target;
        return loss;        
    }
    Matrix& getOutput() {
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;        
        floatType in = input(0,0);                
        floatType y = 0;
        floatType A[order];
        for(size_t i = 0; i < order; i++) A[i] = w(i,0);
        y = chebyshev(in,A,order);        
        output(0,0) = std::clamp(y,-1.0,1.0);                
        return output;
    }
    
    void Formula(Matrix & e, Matrix & w) {        
        std::vector<floatType> v(e.rows());
        floatType A[order];        
        for(size_t j = 0; j < order; j++) 
        {
            A[j] = w(j,0);            
        }        
        for(size_t i = 0; i < e.rows(); i++) {
            floatType in = e(i,0);              
            v[i] = chebyshev(in,A,order) - chebyshev(0,A,order);            
            v[i] = std::clamp(v[i],-1.0,1.0);
        }        
        p.plot_x(v.data(),v.size(),"chebyshev");        
    }    
    
    void PrintEquation(Matrix & w) {        
        for(size_t i = 1; i < w.rows(); i++)
        {
            std::cout << "c" << i << "=" << w(i,0) << std::endl;
        }
    }
};

struct Biquad
{
    floatType z[3];
    floatType p[3];
    floatType x[2];
    floatType y[2];

    Biquad() {
        x[0] = x[1] = 0;
        y[0] = y[1] = 0;    
    }
    void setCoeffs(floatType Z[3], floatType P[3]) {
        memcpy(z,Z,sizeof(z));
        memcpy(p,P,sizeof(p));
    }
    void clear() {
        memset(x,0,sizeof(x));
        memset(y,0,sizeof(y));
    }
    floatType Tick(floatType I)
    {
        floatType r = I*z[0] + x[0]*z[1] + x[1] * z[2] - y[0]*p[0] - y[1]*p[1];
        x[1] = x[0];
        x[0] = I;
        y[1] = y[0];
        y[0] = r;
        return r;
    }
};

struct ButterworthLowpassFilter
{
    int order;
    floatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthLowpassFilter(int order, floatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5);
    }
    void setFilter(floatType f, floatType Q) {
        auto x = IIRFilters::Butterworth::butterlp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate);
        filter.setCoefficients(c);
    }
    void setCutoff(floatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(floatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    floatType Tick(floatType I, floatType A=1, floatType X=1, floatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    void ProcessBlock(size_t n, floatType * in, floatType * out) {
        for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
    }
    std::vector<floatType> impulse_response(size_t n)
    {
        std::vector<floatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
};

struct BiquadLayer : public Layer
{
    Matrix output;
    NeuralNetwork * net;
    Biquad bq;
    int order;
    // calculate the output from the weights
    BiquadLayer(NeuralNetwork * neuralnet,size_t s,ActivationType a,bool useAutoDiff = false) : Layer(OUTPUT,s,a,useAutoDiff),net(neuralnet),order(s) {
        output.resize(1,s);
    }
    Matrix Loss(Matrix & m, Matrix & target) {        
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        Matrix& input = net->layers[0]->input;                        
        auto x = net->GetOutput();
        bq.clear();
        bq.z[0] = w(0,0);
        bq.z[1] = w(1,0);
        bq.z[2] = w(2,0);
        bq.p[0] = w(3,0);
        bq.p[1] = w(4,0);
        bq.p[2] = w(5,0);
        x(0,0) = bq.Tick(1.0);
        for(size_t i = 1; i < order; i++)
            x(0,i) = bq.Tick(0);        
        loss = x-target;
        return loss;        
    }
    void Formula(Matrix & w) {                
        std::vector<floatType> ir(order);        
        bq.z[0] = w(0,0);
        bq.z[1] = w(1,0);
        bq.z[2] = w(2,0);
        bq.p[0] = w(3,0);
        bq.p[1] = w(4,0);
        bq.p[2] = w(5,0);
        bq.clear();
        ir[0] = bq.Tick(1.0);
        for(size_t i = 1; i < order; i++)
            ir[i] = bq.Tick(0);
        p.plot_x(ir.data(),ir.size(),"biquad");
    }    
    void PrintEquation(Matrix & w) {        
        std::cout << "z0=" << w(0,0) << std::endl;
        std::cout << "z1=" << w(1,0) << std::endl;
        std::cout << "z2=" << w(2,0) << std::endl;
        std::cout << "p0=" << w(3,0) << std::endl;
        std::cout << "p1=" << w(4,0) << std::endl;
        std::cout << "p2=" << w(5,0) << std::endl;
    }
    Matrix& getOutput() {
        Matrix& w = net->connections[net->connections.size()-1]->weights;
        bq.clear();
        bq.z[0] = w(0,0);
        bq.z[1] = w(1,0);
        bq.z[2] = w(2,0);
        bq.p[0] = w(3,0);
        bq.p[1] = w(4,0);
        bq.p[2] = w(5,0);
        output(0,0) = bq.Tick(1.0);
        for(size_t i = 1; i < order; i++)
            output(0,i) = bq.Tick(0);        
        return output;
    }
};


struct SoundProcessorLayer
{
    // SoundProcessor 
    // parameters(I,A,X,Y)
};

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

void Quadratic(ActivationType atype) {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    
    for(size_t i = 0; i < 256; i++) {
        examples.push_back((double)i/256.0);
        training.push_back((double)i/256.0);
    }
    p.plot_x(examples.data(),examples.size(),"examples");
    p.plot_x(training.data(),training.size(),"training");

    Matrix e = matrix_new(256,1,examples);
    Matrix t = matrix_new(256,1,training);

    NeuralNetwork net;
    Layer * input = new Layer(INPUT,1,LINEAR);
    Layer * hidden= new Layer(HIDDEN,3,atype);
    QuadraticLayer * output= new QuadraticLayer(&net,1,LINEAR);
    net.addLayer(input);
    net.addLayer(hidden);
    net.addLayer(output);
    net.connect();
    net.loss_widget = 1e-3;
    
    ParameterSet p(e,t,10000,256);
    p.batch_size=4;
    p.learning_rate = 0.0001;
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
    output->Formula(e,m);
    output->PrintEquation(m);
}

inline DspFloatType udo1(DspFloatType x, DspFloatType g = 1.0)
{
    DspFloatType ax = fabs(x*g);
    if(ax == 0) return 0;
    return -std::clamp(((x/ax)*(1-exp(g*(x*x)/ax))),-1.0,1.0);
}

void ChebyWet()
{
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    floatType A[8] = {0.5,0.25,0.1,0.01,0.001,-0.001,-0.0005,0};
    for(int i = 0; i < 128; i++) {
        examples.push_back(std::sin(2*M_PI*(double)i/128.0));
        training.push_back(chebyshev(std::sin(2*M_PI*(double)i/128.0),A,8));
    }
    p.plot_x(examples.data(),examples.size(),"examples");
    p.plot_x(training.data(),training.size(),"training");
    
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

    ParameterSet p(e,t,10000,8);
    p.batch_size=8;
    p.learning_rate = 0.01;
    p.momentum_factor = 0.9;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = RMSPROP_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,STOCHASTIC);
    net.ForwardPass(e);
    
    auto x = net.GetOutput();
    plot(x);
}

// this is hard
void Cheby(ActivationType atype) {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;

    for(int i = 0; i < 128; i++) {
        examples.push_back(std::sin(2*M_PI*(double)i/128.0));
        training.push_back(udo1(std::sin(2*M_PI*(double)i/128.0),1));
    }
    p.plot_x(examples.data(),examples.size(),"examples");
    p.plot_x(training.data(),training.size(),"training");
    
    Matrix e = matrix_new(examples.size(),1,examples);
    Matrix t = matrix_new(training.size(),1,training);
        
    NeuralNetwork net;
    constexpr size_t CHEBYS=8;
    Layer * input = new Layer(INPUT,1,LINEAR);        
    Layer * hidden1= new Layer(HIDDEN,CHEBYS,atype);            
    Layer * hidden2= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden3= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden4= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden5= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden6= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden7= new Layer(HIDDEN,CHEBYS,atype);        
    Layer * hidden8= new Layer(HIDDEN,CHEBYS,atype);            
    ChebyshevLayer * output = new ChebyshevLayer(&net,CHEBYS,1,LINEAR);
    Layer * out = new Layer(OUTPUT,1,atype);
    //Layer * output= new Layer(OUTPUT,1,atype);        
    net.addLayer(input);        
    
    /*
    net.addLayer(hidden1);                        
    net.addLayer(hidden2);                    
    net.addLayer(hidden3);            
    net.addLayer(hidden4); 
    net.addLayer(hidden5);    
    net.addLayer(hidden6);    
    net.addLayer(hidden7);            
    */  
    net.addLayer(hidden8);                        
    net.addLayer(output);    
    net.connect();

    net.loss_widget = 1e-6;

    ParameterSet p(e,t,1000,8);
    p.batch_size= 128;
    p.learning_rate = 0.001;
    p.momentum_factor = 0.001;
    p.regularization_strength = 0;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = GD_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,STOCHASTIC);
    net.ForwardPass(e);
    
    auto m = net.connections[net.connections.size()-1]->weights;    
    output->Formula(e,m);
    output->PrintEquation(m);
}

void Biquad() {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    ButterworthLowpassFilter lp(2,44100.0);
    int ir_size=36;
    training = examples = lp.impulse_response(ir_size);
    p.plot_x(examples.data(),examples.size(),"butterworth");
    Matrix e = matrix_new(1,ir_size,examples);
    Matrix t = matrix_new(1,ir_size,training);
        
    NeuralNetwork net;
    Layer * input = new Layer(INPUT,ir_size,LINEAR);    
    Layer * hidden1= new Layer(HIDDEN,ir_size,RELU);
    Layer * hidden2= new Layer(HIDDEN,ir_size,RELU);
    Layer * hidden3= new Layer(HIDDEN,ir_size,RELU);
    Layer * hidden4= new Layer(HIDDEN,6,RELU);
    BiquadLayer * output= new BiquadLayer(&net,ir_size,RELU);
    net.addLayer(input);    
    net.addLayer(hidden1);
    net.addLayer(hidden2);
    net.addLayer(hidden3);    
    net.addLayer(hidden4);
    net.addLayer(output);    
    net.connect();
    net.loss_widget = 1e-6;

    ParameterSet p(e,t,10000,1);
    p.batch_size=1;
    p.learning_rate = 0.0001;
    p.momentum_factor = 0.9;
    p.regularization_strength = 0;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = false;
    p.optimizer = GD_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,NONSTOCHASTIC);
    net.ForwardPass(e);
    
    auto m = net.connections[net.connections.size()-1]->weights;    
    output->Formula(m);
    output->PrintEquation(m);
}


int main(int argc, char * argv[]) {     
    p.setstyle("lines");
   //XOR2(TANH,0.1,0.9);
   //Line(TANH);
   //Quadratic(TANH);   
   Cheby(FULLWAVE_SIGMOID);
   //ChebyWet();
   //Biquad();
   sleep(5000);
}
