#pragma once

//#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE
#include "Eigen/Core"
//#include "unsupported/Eigen/AutoDiff"
#include <cfloat>
#include <random>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <functional>
#include <queue>
#include <fstream>

#include <boost/math/differentiation/autodiff.hpp>

#define PROVE(x) std::cout << x << std::endl;

typedef float floatType;
using Matrix = Eigen::Matrix<floatType,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;


enum LayerType {
    INPUT=1,
    OUTPUT=2,
    HIDDEN=3,
};

enum {
    GD_OPTIMIZER,
    ADAM_OPTIMIZER,
    RMSPROP_OPTIMIZER,
};

enum LossFunctionType {
    CROSS_ENTROPY_LOSS=1,
    MEAN_SQUARED_ERROR=2,
};

enum ActivationType {
    LINEAR=0,
    SIGMOID=1,
    RELU=2,
    TANH=3,
    SOFTMAX=4,
    ATAN=5,
    BALLS=6,
    ALGEBRA=7,
    SINWAVE=8,
    FULLWAVE_SIGMOID=9,
    FULLWAVE_RELU=10,
    AUTO=11,    
};

typedef void (*activation)(Matrix & m);
typedef void (*activation_grad)(Matrix & m);
using fvar = boost::math::differentiation::detail::fvar<floatType,2>;


inline Matrix matrix_new(size_t rows, size_t cols, std::vector<floatType> & data) {
    Matrix m(rows,cols);        
    for(size_t i = 0; i < rows; i++)
        for(size_t j = 0; j < cols; j++)
            m(i,j) = data[i*cols + j];
    return m;
}
inline Matrix matrix_create(size_t rows, size_t cols) {
    Matrix  m(rows,cols);
    m.setZero();
    return m;
}
inline Matrix createMatrixZeros(size_t rows, size_t cols) {
    return matrix_create(rows,cols);
}

// all the cool kids call it the hadamard product.
inline Matrix hadamard(Matrix & a, Matrix &b)
{
    return Matrix(a.cwiseProduct(b));
}

/// Neural Network functions


inline void linear(Matrix& input) {
    
}

inline void linear_grad(Matrix& input) {
    input.fill(1);
}   


inline void sigmoid(Matrix & m)
{        
    Eigen::Array<floatType,Eigen::Dynamic,Eigen::Dynamic> t = -m.array();    
    m = (1.0 / (1.0 + t.exp()));    
}

inline void fullwave_sigmoid(Matrix & m)
{        
    Eigen::Array<floatType,Eigen::Dynamic,Eigen::Dynamic> t = -m.array();    
    m = 2.0*(1.0 / (1.0 + t.exp()))-1.0;    
}


inline void sigmoid_grad(Matrix & m)
{    
    Matrix t(-m);
    t = (1.0 / (1.0 + t.array().exp())); 
    m = (t.array() * ( 1.0 - t.array() ));    
}

inline void tanh(Matrix & m)
{
    m = m.array().tanh();
}

inline void tanh_grad(Matrix & m)
{    
    Matrix t = m.array().tanh();
    m = (1.0 - t.array()*t.array());
}

inline void atan(Matrix & m)
{
    m = atan(m.array());
}

inline void atan_grad(Matrix & m)
{
    Matrix t = m.array().atan();
    t *= t;
    m = 1.0 / (1.0+t.array());
}

inline void balls(Matrix & m)
{
    m = 2.0 * atan(tanh(m.array()/2.0));
}

inline void algebra(Matrix & m)
{
    m = m.array() / (sqrt(1.0 + pow(m.array(),2.0)));
}

inline void relu(Matrix & m)
{    
    m = m.cwiseMax(0).eval();   
    m = m.cwiseMin(1.0).eval();    
}

inline void fullwave_relu(Matrix & m)
{    
    m = m.cwiseMax(0).eval();            
    m = 2 * m.array()-1;
}

inline void relu_grad(Matrix & m)
{    
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            floatType x = m(i,j);
            if(x > FLT_MIN) m(i,j) = 1;
            else m(i,j) = 0;
        }    
}

inline void softmax(Matrix & m)
{                
    int i;
    for (i = 0; i < m.rows(); i++){
        floatType summed = 0;
        int j;
        for (j = 0; j < m.cols(); j++){
            summed += std::exp(m(i, j));
        }
        for (j = 0; j < m.cols(); j++){
            m(i, j) =  std::exp(m(i, j)) / summed;
        }
    }
}

inline void sinwave(Matrix & m)
{
    m = (cos(2*M_PI*m.array()));    
}

// random already exists somewhere.
inline floatType randr(floatType min = 0.0f, floatType max = 1.0f) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return distribution(generator);
}


struct RandomNumbers
{
    typedef std::chrono::high_resolution_clock myclock;
    unsigned seed;
    std::default_random_engine generator;

    RandomNumbers() {
        myclock::time_point beginning = myclock::now();
        myclock::duration d = myclock::now() - beginning;
        seed = d.count();    
        generator = std::default_random_engine(seed);
    }

    void set_seed(unsigned s) { seed = s; }
    void reseed() {
        myclock::time_point beginning = myclock::now();
        myclock::duration d = myclock::now() - beginning;
        seed = d.count();    
        generator = std::default_random_engine(seed);
    }
    floatType random(floatType min=0.0f, floatType max=1.0f) {
        std::uniform_real_distribution<double> distribution(min,max);
        return distribution(generator);
    }
};

struct BoxMuller {
    floatType z0,z1;
    bool  do_generate;
    RandomNumbers random;

    BoxMuller() {
        z0=z1=0.0;
        do_generate = false;
    }

    floatType generate() {
        floatType epsilon = FLT_MIN;
        floatType two_pi  = 2 * M_PI;
        do_generate = !do_generate;
        if(!do_generate) return z1;
        floatType u1 = random.random();
        floatType u2 = random.random();
        while(u1 <= epsilon) {
            u1 = randr();            
        }
        while(u2 <= epsilon) {
            u1 = randr();            
        }
        z0 = std::sqrt(-2.0f * std::log(u1)) * std::cos(two_pi * u2);
        z1 = std::sqrt(-2.0f * std::log(u1)) * std::sin(two_pi * u2);
        return z0;
    }
};


struct Connection;
struct Network;

struct Activator
{
    activation       activate_f;
    activation_grad  activate_grad_f;        
    bool             useAutoDiff=true;

    Activator() = default;
    Activator(activation af, activation_grad grad, bool ad) : activate_f(af),activate_grad_f(grad),useAutoDiff(ad) {

    }
    void autodiff(Matrix & m)
    {        
        using namespace boost::math::differentiation;
        constexpr unsigned Order = 2;                  // Highest order derivative to be calculated.
        Matrix t = m;
        activate_f(t);
        
        for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            auto const x = make_fvar<double, Order>((t(i,j)));  // Find derivatives at x=2.                            
            auto const y = x;//deriv(x);
            m(i,j) = y.derivative(1);    
            if(std::isnan(m(i,j))) m(i,j) = 0;
            if(std::isinf(m(i,j))) m(i,j) = x < 0? -1:1;            
        }                    
    } 

    void Activate(Matrix& tmp) {                
        activate_f(tmp);                        
    }
    void Grad(Matrix & tmp) {
        if(!useAutoDiff) activate_grad_f(tmp);        
        else autodiff(tmp);
    }
};



struct Layer {
    LayerType        type;        
    Matrix           input,loss,output,temp;        
    ActivationType   atype;
    Activator        activator;
    size_t 			 size;
	// a neuron is a layer with only 1 input (s=1)
	Layer() = default;
	
    Layer(LayerType t,size_t s,ActivationType a,bool useAutoDiff = true) {
        type = t;       
        size = s; 
        activation       activate_f;
        activation_grad  activate_grad_f;                        
        atype = a;
        switch(a) {
			case AUTO:activate_f = linear;            
                         activate_grad_f = linear_grad;                         
                         useAutoDiff = true;
                         break;
            case LINEAR: activate_f = linear;            
                         activate_grad_f = linear_grad;                         
                         break;
            case SIGMOID: activate_f = sigmoid;
                          activate_grad_f = sigmoid_grad;                          
                          break;                          
            case RELU:  activate_f = relu;
                        activate_grad_f = relu_grad;                        
                        break;
            case FULLWAVE_SIGMOID: activate_f = fullwave_sigmoid;
                          activate_grad_f = tanh_grad;                          
                          break;                          
            case FULLWAVE_RELU:  
                        activate_f = fullwave_relu;
                        activate_grad_f = tanh_grad;                        
                        break;
            case TANH:  activate_f = tanh;
                        activate_grad_f = tanh_grad;                        
                        break;
            case ATAN:  activate_f = atan;
                        activate_grad_f = atan_grad;
                        break;
            case BALLS:  activate_f = balls;
                        activate_grad_f = atan_grad;
                        break;                        
            case ALGEBRA:  activate_f = algebra;
                        activate_grad_f = atan_grad;
                        break;                        
            case SINWAVE: activate_f = sinwave;
                        activate_grad_f = atan_grad;
                        break;                        
            case SOFTMAX:
                        activate_f = softmax;
                        activate_grad_f = linear_grad;
                        break;
			
				
        }        
        activator.activate_f = activate_f;
        activator.activate_grad_f = activate_grad_f;
        activator.useAutoDiff = useAutoDiff;
        input.resize(1,s);
        loss.resize(1,s);
    }    
    ~Layer() {

    }
    //virtual size_t size() const { return 1; }
    virtual Layer* operator[](size_t i) { return this; }
    
    virtual void Activate(Matrix& tmp) {                
		input = tmp.eval(); 
		temp  = tmp.eval();
        activator.Activate(tmp);                                              
    }
    virtual void Grad(Matrix & tmp) {		
        activator.Grad(tmp);
    }    
    // loss = output - target
    virtual Matrix Loss(Matrix & output, Matrix & target) {                
        loss = output - target;
        return loss;
    }

    virtual void getInput(Matrix& m) { m = input.eval(); }
    virtual void getOutput(Matrix& m) { m = input.eval(); }
    virtual void getLoss(Matrix & m) { m = loss.eval();  }
    virtual void Update(Matrix & e) {}    
    
    void write(std::ofstream& f)
    {		
		f << type << std::endl;
		f << atype << std::endl;
		f << size << std::endl;
		f << activator.useAutoDiff << std::endl;
	}
};


struct Connection {
    Layer * from,
          * to;


    Matrix weights;
    Matrix bias;

    Connection(std::ifstream & f, Layer * from, Layer *to)
    {
		this->from = from;
		this->to   = to;
		weights = matrix_create(from->size,to->size);
        bias    = matrix_create(1,to->size);
        for(size_t i = 0; i < weights.rows(); i++)
        for(size_t j = 0; j < weights.cols(); j++)
			f >> weights(i,j);
		for(size_t i = 0; i < bias.rows(); i++)
        for(size_t j = 0; j < bias.cols(); j++)
			f >> bias(i,j);
	}
    Connection(Layer * from, Layer * to) {
        this->from = from;
        this->to   = to;
        weights = matrix_create(from->size,to->size);
        bias    = matrix_create(1,to->size);
        bias.fill(1.0f);

        BoxMuller bm;
        srand(time(NULL));
        for(size_t i = 0; i < weights.rows(); i++)
        {            
            for(size_t j = 0; j < weights.cols(); j++)
            {                
                weights(i,j) = bm.generate()/std::sqrt(weights.rows());                
                //weights(i,j) = (((double)std::rand() / (double)RAND_MAX)-0.5)*2.0;
            }
        }
        
    }
    ~Connection() {

    }
    void write(std::ofstream &f) {
		f << weights << std::endl;
		f << bias << std::endl;
	}
};



struct ParameterSet {
    Matrix data;
    Matrix classes;
    LossFunctionType loss_function;
    size_t batch_size;
    floatType learning_rate;
    floatType search_time;
    floatType regularization_strength;
    floatType momentum_factor;
    size_t max_iters;
    bool shuffle;
    bool verbose;
    size_t ticks;
    floatType gamma1=0.9;
    floatType gamma2=0.995;
    int optimizer = GD_OPTIMIZER;
    bool reshuffle=false;

    ParameterSet( Matrix &d, Matrix &c,
                 size_t epochs, size_t bs,
                 LossFunctionType loss=MEAN_SQUARED_ERROR,
                 floatType lr = 0.01, floatType st = 0.0,
                 floatType rs=0.0,floatType m=0.2, bool s=true, bool v=true) {
            max_iters = epochs;
            data = d;
            classes = c;
            loss_function = loss;
            batch_size = bs;
            learning_rate = lr;
            search_time = st;
            regularization_strength = rs;
            momentum_factor = m;
            shuffle = s;
            verbose = v;
            ticks=10;
    }
};

struct Batch {
    Matrix example;
    Matrix training;

    Batch(Matrix & e, Matrix & c) {        
        Matrix x = e;
        Matrix y = c;
        example    = x.eval();
        training   = y.eval();
    }
    Batch(const Batch & b) {
        Matrix x = b.example;
        Matrix y = b.training;
        example = x.eval();
        training= y.eval();
    }
    Batch& operator = (const Batch & b) {
        Matrix x = b.example;
        Matrix y = b.training;
        example = x.eval();
        training= y.eval();
        return *this;
    }
};

Matrix addToEachRow(const Matrix& m, const Matrix & v)
{
    Matrix r(m);
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols();j++)
            r(i,j) += v(0,j);
    return r;
}

enum {
        NONSTOCHASTIC,
        STOCHASTIC
    };
    
struct Network {
protected:
    Network() = default;

public:    
		
    size_t num_features;
    size_t num_outputs;
    std::vector<Layer*> layers;
    std::vector<Connection*> connections;
    std::vector<std::vector<Batch>> batch_list;
    std::vector<Matrix> errori;
    std::vector<Matrix*> weights;
    std::vector<Matrix*> bias;
    std::vector<Matrix> dWi;
    std::vector<Matrix> dbi;
    std::vector<Matrix> sdw;
    std::vector<Matrix> sdb;
    std::vector<Matrix> vdw;
    std::vector<Matrix> vdb;
    std::vector<Matrix> regi;
    std::vector<Matrix> wTi;
    std::vector<Matrix> errorLastTi;
    std::vector<Matrix> fprimei;
    std::vector<Matrix> inputTi;
    std::vector<Matrix> dWi_avg;
    std::vector<Matrix> dbi_avg;
    std::vector<Matrix> dWi_last;
    std::vector<Matrix> dbi_last;
    
    floatType loss = 1e6;
    floatType loss_widget=1e-6;
    bool delete_data=true;
        
    Network(size_t num_features,
            std::vector<int64_t> & hidden,
            std::vector<ActivationType> & activations,
            size_t num_outputs,
            ActivationType output_activation
            )
    {
        assert(num_features > 0);
        assert(num_outputs > 0);
        this->num_features = num_features;
        this->num_outputs  = num_outputs;
        size_t num_hidden = hidden.size();
        size_t num_layers = 2 + num_hidden;
        layers.resize(num_layers);

        for(size_t i = 0; i < num_layers; i++)
        {
            Layer * ln = NULL;
            if(i == 0)
                ln = new Layer(INPUT,num_features,LINEAR);
            else if(i == num_layers-1)
                ln = new Layer(OUTPUT,num_outputs,output_activation);
            else
                ln = new Layer(HIDDEN,hidden[i-1], activations[i-1]);
            assert(ln != NULL);
            layers[i] = ln;
        }
        size_t num_connections = num_layers-1;
        for(size_t i = 0; i < num_connections; i++)
        {
            assert(layers[i] != NULL);
            assert(layers[i+1]!= NULL);
            Connection * c = new Connection(layers[i],layers[i+1]);
            weights.push_back(&c->weights);
            bias.push_back(&c->bias);
            connections.push_back(c);
        }
    }
    ~Network() {
        if(delete_data)
        {
            for(size_t i = 0; i < layers.size(); i++)
                delete layers[i];
            for(size_t i = 0; i < connections.size(); i++)
                delete connections[i];
        }
    }

    size_t NumLayers() const { return layers.size(); }
    size_t NumConnections() const { return connections.size(); }
    size_t NumInputs() const { return num_features; }
    size_t NumOutputs() const { return num_outputs; }
    
    /*
	Layer* L(size_t i) { return layers[i]; }
	Connection* C()(size_t i) { return connections[i]; }
	Matrix& W(size_t i) { return *weights[i]; }
	Matrix& B(size_t i) { return *bias[i]; }
	*/
    void ForwardPass(Matrix& input) {
        assert(input.cols() == layers[0]->input.cols());
        layers[0]->input = input.eval();
        Matrix tmp,tmp2;        
        for(size_t i = 0; i < connections.size(); i++)
        {         
            tmp  = layers[i]->input*connections[i]->weights;            
            tmp2 = addToEachRow(tmp,connections[i]->bias);
            connections[i]->to->Activate(tmp2);       
        }
    }    
    void BackwardPass(Matrix error, floatType clr = 0.01) {
		size_t layer = layers.size()-1;
		Layer* to = layers[layer];
		Connection* con = connections[layer-1];
		
		errori[layer] = error;
		dWi[layer-1] = con->from->input.transpose() * errori[layer];
		dbi[layer-1] = errori[layer].eval();
												
		for(layer = layers.size()-2; layer > 0; layer--)
		{                                                     
			size_t hidden_layer = layer-1;
			to  = layers[layer];
			con = connections[layer-1];
			wTi[hidden_layer] = connections[layer]->weights.transpose();                            
			errorLastTi[hidden_layer] = errori[layer+1]*wTi[hidden_layer];
			fprimei[hidden_layer] = con->to->input.eval();
			con->to->Grad(fprimei[hidden_layer]);
			errori[layer] = hadamard(errorLastTi[hidden_layer],fprimei[hidden_layer]);
			inputTi[hidden_layer] = con->from->input.transpose();                            
			dWi[hidden_layer] = inputTi[hidden_layer] * errori[layer];
			dbi[hidden_layer] = errori[layer].eval();
		}                                                                
		for(size_t idx=0; idx < connections.size(); idx++) {
			dWi_avg[idx] = dWi[idx] + dWi_avg[idx];
			dbi_avg[idx] = dbi[idx] + dbi_avg[idx];                    
		}                                                                                 
        for(size_t idx = 0; idx < connections.size(); idx++)
        {
            dWi_avg[idx] = dWi_avg[idx] * clr;
            dbi_avg[idx] = dbi_avg[idx] * clr;

            connections[idx]->weights = -dWi_avg[idx] + connections[idx]->weights;
            connections[idx]->bias    = -dbi_avg[idx] + connections[idx]->bias;                    

            dWi_avg[idx].setZero();
            dbi_avg[idx].setZero();
        }
	}	
    floatType CrossEntropyLoss(Matrix& prediction, Matrix& actual, floatType rs) {
        floatType total_err = 0;
        floatType reg_err = 0;
        total_err = (actual * prediction.array().log().matrix()).sum();        
        if(rs > 0)
            for(size_t i = 0; i < connections.size(); i++)
            {
                Matrix & weights = connections[i]->weights;
                reg_err += hadamard(weights,weights).sum();
            }        
        return (-1.0 / actual.rows()*total_err) + rs*0.5*reg_err;
    }
    floatType MeanSquaredError(Matrix& loss, floatType rs) {
        floatType total_err = 0;
        floatType reg_err = 0;                
        total_err = hadamard(loss,loss).sum();        
        if(rs > 0)
            for(size_t i = 0; i < connections.size(); i++)
            {
                Matrix & w = connections[i]->weights;
                reg_err += hadamard(w,w).sum();
            }
            
        return ((0.5 / loss.rows()) * total_err) + (rs*0.5*reg_err);
    }
    virtual void GetInput( Matrix& m) {
        m = layers[0]->input.eval();
    }
    virtual void GetOutput(Matrix& m) {
        LastLayer()->getOutput(m);
    }
    virtual void GetLoss(Matrix& m) {
        LastLayer()->getLoss(m);
    }
    Connection* LastConnection() { return connections[connections.size()-1]; }
    Layer* LastLayer() const { return layers[layers.size()-1]; }
    
    // legacy
    std::vector<int> predict() {
        Layer* output_layer = layers[layers.size()-1];
        std::vector<int> prediction;
        prediction.resize(output_layer->input.rows());
        Matrix & input = output_layer->input;
        for(size_t i = 0; i < input.rows(); i++) {
            int max = 0;
            for(size_t j = 0; j < input.cols(); j++) {
                if(input(i,j) > input(i,max)) max = j;
            }
            prediction[i] = max;
        }
        return prediction;
    }
    floatType accuracy(Matrix & data, Matrix & classes) {
        ForwardPass(data);
        std::vector<int> p = predict();
        floatType num_correct = 0;
        for(size_t i = 0; i < data.rows(); i++) {
            if(classes(i,p[i]) == 1)
                num_correct++;
        }
        return 100*num_correct/classes.rows();
    }
    void shuffle_batches() {
		#pragma omp parallel for
        for(size_t i = 0; i < batch_list.size(); i++) {
              std::random_shuffle(batch_list[i].begin(),batch_list[i].end());
              //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
              //std::shuffle (batch_list[i].begin(),batch_list[i].end(), std::default_random_engine(seed));
        }
    }
    
    void generate_batches(size_t num_batches,
                            size_t batch_size,
                            Matrix & data,
                            Matrix & classes,
                            bool shuffle) {
        size_t rc = 0;
        batch_list.clear();
        for(size_t i = 0; i < num_batches; i++) {
            std::vector<Batch> l;
            size_t cur_batch_size = batch_size;
            if(i == num_batches) {
                if( data.rows() % batch_size != 0) {
                    cur_batch_size = data.rows() % batch_size;
                }
            }
            for(size_t j = 0; j < cur_batch_size; j++) {
                Matrix e = data.row(rc);                
                Matrix c = classes.row(rc);                
                Batch b(e,c);                
                l.push_back(b);
                rc = rc + 1;
                rc = rc % data.rows();
            }
            batch_list.push_back(l);
        }
        //if(shuffle) shuffle_batches();
    }
    void gd_optimize(const ParameterSet & ps, size_t epoch, size_t rows)
    {
        floatType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }
        
        floatType clr = currentLearningRate / rows;
        #pragma omp parallel for
        for(size_t idx = 0; idx < connections.size(); idx++)
        {
            dWi_avg[idx] = dWi_avg[idx] * clr;
            dbi_avg[idx] = dbi_avg[idx] * clr;

            if(ps.regularization_strength > 0) {
                regi[idx] = connections[idx]->weights * ps.regularization_strength;
                dWi_avg[idx] = regi[idx] + dWi_avg[idx];
            }
            if(ps.momentum_factor > 0) {
                dWi_last[idx] = dWi_last[idx] * ps.momentum_factor;
                dbi_last[idx] = dbi_last[idx] * ps.momentum_factor;                    
                dWi_avg[idx] = (dWi_last[idx] + dWi_avg[idx]);
                dbi_avg[idx] = (dbi_last[idx] + dbi_avg[idx]);         
            }           

            connections[idx]->weights = -dWi_avg[idx] + connections[idx]->weights;
            connections[idx]->bias    = -dbi_avg[idx] + connections[idx]->bias;                    

            if(ps.regularization_strength > 0 || ps.momentum_factor > 0) {
                dWi_last[idx] = dWi_avg[idx].eval();
                dbi_last[idx] = dbi_avg[idx].eval();                
            }
            dWi_avg[idx].setZero();
            dbi_avg[idx].setZero();
        }
    }

    void adam_optimize(const ParameterSet & ps, size_t epoch, size_t rows)
    {   
        floatType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }
        
        floatType alpha = currentLearningRate / rows;
        #pragma omp parallel for
        for(size_t i = 0; i < connections.size(); i++)
        {
            
            dWi_avg[i] = (alpha/10.0)*dWi_avg[i];
            dbi_avg[i] = (alpha/10.0)*dbi_avg[i];
                                                
            if(ps.regularization_strength > 0) {
                regi[i]     = connections[i]->weights*ps.regularization_strength;
                dWi_avg[i]  = regi[i] + dWi_avg[i];                        
            }
            if(ps.momentum_factor > 0) {            
                dWi_last[i]=(ps.momentum_factor/10.0)*dWi_last[i];
                dbi_last[i]=(ps.momentum_factor/10.0)*dbi_last[i];
                dWi_avg[i] =dWi_last[i]+dWi_avg[i];
                dbi_avg[i] =dbi_last[i]+dbi_avg[i];
            }
            
            sdw[i] = ps.gamma1 * sdw[i] + (1-ps.gamma1)*dWi_avg[i];
            sdb[i] = ps.gamma1 * sdb[i] + (1-ps.gamma1)*dbi_avg[i];
                        
            vdw[i] = ps.gamma2 * vdw[i] + (1-ps.gamma2) * hadamard(dWi_avg[i],dWi_avg[i]); 
            vdb[i] = ps.gamma2 * vdb[i] + (1-ps.gamma2) * hadamard(dbi_avg[i],dbi_avg[i]); 
            
            Matrix mdw_corr = sdw[i] / (1 - pow(ps.gamma1,epoch+1));
            Matrix mdb_corr = sdb[i] / (1 - pow(ps.gamma1,epoch+1));

            Matrix vdw_corr = vdw[i] / (1 - pow(ps.gamma2,epoch+1));
            Matrix vdb_corr = vdb[i] / (1 - pow(ps.gamma2,epoch+1));

            vdw_corr = vdw_corr.cwiseMax(1e-08);
            vdb_corr = vdb_corr.cwiseMax(1e-08);
            
            Matrix m1 = (alpha / (sqrt(vdw_corr.array())));
            Matrix m2 = (alpha / (sqrt(vdb_corr.array())));
            
            connections[i]->weights = connections[i]->weights - hadamard(m1,mdw_corr);
            connections[i]->bias    = connections[i]->bias    - hadamard(m2,mdb_corr);    
            
            if(ps.momentum_factor > 0 || ps.regularization_strength > 0)
            {
                dWi_last[i] = dWi_avg[i].eval();
                dbi_last[i] = dbi_avg[i].eval();                                
            }      
            
            dWi_avg[i].setZero();
            dbi_avg[i].setZero();
        }
    }       
    void rmsprop_optimize(const ParameterSet & ps, size_t epoch, size_t rows)
    {                    

        floatType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }
        
        floatType alpha = currentLearningRate / rows;
        #pragma omp parallel for
        for(size_t i = 0; i < connections.size(); i++)
        {
            
            dWi_avg[i] = alpha*dWi_avg[i];
            dbi_avg[i] = alpha*dbi_avg[i];
                                                
            if(ps.regularization_strength > 0) {
                regi[i]     = connections[i]->weights*ps.regularization_strength;
                dWi_avg[i]  = regi[i] + dWi_avg[i];                        
            }
            if(ps.momentum_factor > 0) {            
                dWi_last[i]=ps.momentum_factor*dWi_last[i];
                dbi_last[i]=ps.momentum_factor*dbi_last[i];            
                dWi_avg[i] =dWi_last[i]+dWi_avg[i];
                dbi_avg[i] =dbi_last[i]+dbi_avg[i];
            }
                

            vdw[i] = ps.gamma1 * vdw[i] + (1-ps.gamma1)*hadamard(dWi_avg[i],dWi_avg[i]);
            vdb[i] = ps.gamma1 * vdb[i] + (1-ps.gamma1)*hadamard(dbi_avg[i],dbi_avg[i]);

            Matrix vdw_corr = vdw[i] / (1 - pow(ps.gamma1,epoch+1));
            Matrix vdb_corr = vdb[i] / (1 - pow(ps.gamma1,epoch+1));

            vdw_corr = vdw_corr.cwiseMax(1e-08);
            vdb_corr = vdb_corr.cwiseMax(1e-08);
            
            Matrix m1 = (alpha / (sqrt(vdw_corr.array())));
            Matrix m2 = (alpha / (sqrt(vdb_corr.array())));
            
            connections[i]->weights = connections[i]->weights - hadamard(m1,dWi_avg[i]);
            connections[i]->bias    = connections[i]->bias    - hadamard(m2,dbi_avg[i]);    
            
            if(ps.momentum_factor > 0 || ps.regularization_strength > 0)
            {
                dWi_last[i] = dWi_avg[i].eval();
                dbi_last[i] = dbi_avg[i].eval();                                
            }           
            dWi_avg[i].setZero();
            dbi_avg[i].setZero();
        }
    }                        

    void learn(size_t batch, size_t cur_batch_size)
    {
        for(size_t training = 0; training < cur_batch_size; training++)
        {
            Matrix& example = batch_list[batch][training].example;
            Matrix& target  = batch_list[batch][training].training;                    
            ForwardPass(example);
            
            size_t layer = layers.size()-1;
            Layer* to = layers[layer];
            Connection* con = connections[layer-1];
            
            errori[layer] = to->Loss(to->input,target);   
            dWi[layer-1] = con->from->input.transpose() * errori[layer];
            dbi[layer-1] = errori[layer].eval();
                                                    
            for(layer = layers.size()-2; layer > 0; layer--)
            {                                                     
                size_t hidden_layer = layer-1;
                to  = layers[layer];
                con = connections[layer-1];
                wTi[hidden_layer] = connections[layer]->weights.transpose();                            
                errorLastTi[hidden_layer] = errori[layer+1]*wTi[hidden_layer];
                fprimei[hidden_layer] = con->to->input.eval();
                con->to->Grad(fprimei[hidden_layer]);
                errori[layer] = hadamard(errorLastTi[hidden_layer],fprimei[hidden_layer]);
                inputTi[hidden_layer] = con->from->input.transpose();                            
                dWi[hidden_layer] = inputTi[hidden_layer] * errori[layer];
                dbi[hidden_layer] = errori[layer].eval();
            }                                                                         
            for(size_t idx=0; idx < connections.size(); idx++) {
                dWi_avg[idx] = dWi[idx] + dWi_avg[idx];
                dbi_avg[idx] = dbi[idx] + dbi_avg[idx];                    
            }                                                     
        }
    }
    virtual void verbosity(const ParameterSet & ps, size_t epoch, Matrix & data, Matrix & classes) {
        if(ps.verbose == true) {
            if(epoch % ps.ticks == 0 || epoch <= 1) {
                Matrix tmp;
                ForwardPass(data);
                if(ps.loss_function == CROSS_ENTROPY_LOSS) {
                    GetOutput(tmp);
                    printf("EPOCH: %ld loss is %f\n",epoch, CrossEntropyLoss(tmp,classes,ps.regularization_strength));
                }
                else {
                    GetLoss(tmp);
                    printf("EPOCH: %ld loss is %f\n",epoch, loss=MeanSquaredError(tmp,ps.regularization_strength));                    
                }
            }
        }
    }
    void create_matrix()
    {
        Matrix beforeOutputT = createMatrixZeros(layers[layers.size()-2]->size,1);
        		
        for(size_t i = 0; i < connections.size(); i++) {
			assert(layers[i]);
			assert(connections[i]);
            errori.push_back(createMatrixZeros(1,layers[i]->size));
            dWi.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
            dbi.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            sdw.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
            sdb.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            vdw.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
            vdb.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            regi.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
        }
        
        errori.push_back(createMatrixZeros(1,LastLayer()->size));
        size_t num_hidden = layers.size()-2;
        
        for(size_t k = 0; k < num_hidden; k++)
        {
            wTi.push_back(createMatrixZeros(connections[k+1]->weights.cols(),connections[k+1]->weights.rows()));
            errorLastTi.push_back(createMatrixZeros(1,wTi[k].cols()));
            fprimei.push_back(createMatrixZeros(1,connections[k]->to->size));
            inputTi.push_back(createMatrixZeros(connections[k]->from->size,1));
        }
        for(size_t i = 0; i < connections.size(); i++) {
            dWi_avg.push_back(createMatrixZeros(connections[i]->weights.rows(),connections[i]->weights.cols()));
            dbi_avg.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            dWi_last.push_back(createMatrixZeros(connections[i]->weights.rows(),connections[i]->weights.cols()));
            dbi_last.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
        }
    }

    // if you turn shuffle off it is the same as non-stochastic
    void stochastic(ParameterSet & ps) {
        
        Matrix & data = ps.data;
        Matrix & classes = ps.classes;
        size_t num_batches = data.rows() / ps.batch_size;
        if(data.rows() % ps.batch_size != 0) num_batches++;
        size_t epoch = 0;


        create_matrix();
        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);

        while(epoch <= ps.max_iters && loss > loss_widget) {
            if(ps.shuffle) {
                shuffle_batches();
            }             
            
            for(size_t batch = 0; batch < num_batches; batch++) 
            {
                size_t cur_batch_size = ps.batch_size;
                
                if(batch == num_batches) {
                    if(data.rows() % ps.batch_size != 0) {
                        cur_batch_size = data.rows() % ps.batch_size;
                    }
                }
                learn(batch,cur_batch_size);                
                if(ps.reshuffle) {
                    shuffle_batches();
                }                        
            }
            if(ps.optimizer == ADAM_OPTIMIZER) adam_optimize(ps,epoch,data.rows());
            else if(ps.optimizer == RMSPROP_OPTIMIZER) rmsprop_optimize(ps,epoch,data.rows());
            else gd_optimize(ps,epoch,data.rows());
            verbosity(ps,epoch,data,classes);
            epoch++;
        }
    }

    // the only difference is shuffle is always off
    void nonstochastic(ParameterSet & ps) {
        
        Matrix & data = ps.data;
        Matrix & classes = ps.classes;
        size_t num_batches = data.rows();        
        size_t epoch = 0;

        create_matrix();
        ps.batch_size=1;
        ps.shuffle=false;        
        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);
        while(epoch <= ps.max_iters && loss > loss_widget) {                    
            for(size_t batch=0; batch < num_batches; batch++) learn(batch,ps.batch_size);                            
            if(ps.optimizer == ADAM_OPTIMIZER) adam_optimize(ps,epoch,data.rows());
            else if(ps.optimizer == RMSPROP_OPTIMIZER) rmsprop_optimize(ps,epoch,data.rows());
            else gd_optimize(ps,epoch,data.rows());
            verbosity(ps,epoch,data,classes);
            epoch++;
        }
    }

    void train(ParameterSet & ps, int type) {
        if(type == NONSTOCHASTIC) nonstochastic(ps);
        else stochastic(ps);
    }
    
    void save(const char * filename) {
		std::ofstream f;
		f.open(filename,std::ofstream::trunc);
		f << num_features << std::endl;
		f << num_outputs  << std::endl;
		f << layers.size() << std::endl;
		f << connections.size() << std::endl;
		for(size_t i = 0; i < layers.size(); i++)
			layers[i]->write(f);
		for(size_t i = 0; i < connections.size(); i++)
			connections[i]->write(f);
		f.close();
	}
	void load(const char * filename) {
		std::ifstream f;
		size_t layer_size,connection_size;
		f.open(filename);
		f >> num_features;
		f >> num_outputs;
		f >> layer_size;		
		f >> connection_size;
		layers.resize(layer_size);
		for(size_t i = 0; i < layers.size(); i++)
		{
			int     type;            
			int		atype;			
			size_t 	size;
			bool 	autodiff;
			f >> type;
			f >> atype;
			f >> size;			
			f >> autodiff;
			layers[i] = new Layer((LayerType)type,size,(ActivationType)atype,autodiff);
		}
		connections.resize(connection_size);
		for(size_t i = 0; i < connections.size(); i++)
		{
			connections[i] = new Connection(f,layers[i],layers[i+1]);
			weights.push_back(&connections[i]->weights);
            bias.push_back(&connections[i]->bias);            
		}
		f.close();
	}
	
};


struct NeuralNetwork : public Network
{       
    NeuralNetwork() : Network()
    {
        delete_data = false;
    }
    ~NeuralNetwork() {

    }
    void addLayer(Layer * pL)
    {
        layers.push_back(pL);
    }
    void connect() {
        size_t num_layers = layers.size();
        size_t num_connections = num_layers;
        num_features = layers[0]->size;
        num_outputs  = layers[layers.size()-1]->size;
        for(size_t i = 0; i < num_connections-1; i++)
        {
            assert(layers[i] != NULL);
            assert(layers[i+1]!= NULL);
            Connection * c = new Connection(layers[i],layers[i+1]);
            connections.push_back(c);
        }
    }            
};
