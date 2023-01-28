#include "Octopus.hpp"
#include <cfloat>
#include <random>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

#define PROVE(x) std::cout << x << std::endl;
Octave::Octopus interp;
typedef float floatType;

enum LayerType {
    INPUT=1,
    OUTPUT=2,
    HIDDEN=3,
};

enum LossFunctionType {
    CROSS_ENTROPY_LOSS=1,
    MEAN_SQUARED_ERROR=2,
};


typedef void (*activation)(MatrixXf & m);
typedef void (*activation_grad)(MatrixXf & m);


inline MatrixXf matrix_new(size_t rows, size_t cols, std::vector<floatType> & data) {
    MatrixXf m(rows,cols);        
    for(size_t i = 0; i < rows; i++)
        for(size_t j = 0; j < cols; j++)
            m(i,j) = data[i*cols + j];
    return m;
}
inline MatrixXf matrix_create(size_t rows, size_t cols) {
    MatrixXf  m(rows,cols);
    m.fill(0);
    return m;
}
inline MatrixXf createMatrixZeros(size_t rows, size_t cols) {
    return matrix_create(rows,cols);
}

// all the cool kids call it the hadamard product.
MatrixXf hadamard(MatrixXf & a, MatrixXf &b)
{    
    ValueList v;
    v(0) = a;
    v(1) = b;
    v = interp("hadamard",v,1);
    return v(0).matrix_value();
}

/// Neural Network functions


void linear(MatrixXf& input) {
    
}

void linear_grad(MatrixXf& input) {
    input.fill(1);
}   


void sigmoid(MatrixXf & m)
{            
    ValueList v;
    v(0) = m;    
    v = interp("sigmoid",v,1);
    m = v(0).matrix_value();
}

void fullwave_sigmoid(MatrixXf & m)
{                
    sigmoid(m);
    m = 2.0*m-1.0;
}

void sigmoid_grad(MatrixXf & m)
{   
    ValueList v;
    v(0) = m;    
    v = interp("sigmoid_grad",v,1);
    m = v(0).matrix_value();
}


void tanh(MatrixXf & m)
{
    ValueList v;
    v(0) = m;
    v = Octave::tanh(v);
    m = v(0).matrix_value();
}


void tanh_grad(MatrixXf & m)
{        
    MatrixXf t = m;
    tanh(t);
    m = (1 - hadamard(t,t));    
}

void atan(MatrixXf & m)
{
    ValueList v;
    v(0) = m;
    v = Octave::atan(v);
    m = v(0).matrix_value();
}

void atan_grad(MatrixXf & m)
{
    MatrixXf t = m;
    atan(t);
    ValueList v;
    v(0) = m.rows();
    v(1) = m.cols();
    m = 1.0/(1.0+hadamard(t,t));
}

void balls(MatrixXf & m)
{    
    ValueList v;
    v(0) = m;
    v = interp("gunderballs.m",v,1);
    m = v(0).matrix_value();
}
void algebra(MatrixXf & m)
{
    ValueList v;
    v(0) = m;
    v = interp("algebra",v,1);
    m = v(0).matrix_value();
}



void relu(MatrixXf & m)
{    
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            if(m(i,j) < 0) m(i,j) = 0;
            if(m(i,j) > 1) m(i,j) = 1;
        }    
}

void fullwave_relu(MatrixXf & m)
{    
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            if(m(i,j) < -1) m(i,j) = -1;
            if(m(i,j) > 1) m(i,j) = 1;
        }    
}

void relu_grad(MatrixXf & m)
{    
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            floatType x = m(i,j);
            if(x > FLT_MIN) m(i,j) = 1;
            else m(i,j) = 0;
        }    
}


void softmax(MatrixXf & m)
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

void sinwave(MatrixXf & m)
{
    ValueList v;
    v(0) = 2*M_PI*m;
    v = Octave::sin(v);
    m = v(0).matrix_value();
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

// random already exists somewhere.
floatType randr(floatType min = 0.0f, floatType max = 1.0f) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return distribution(generator);
}


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
};

struct Connection;





struct Layer {
    LayerType        type;
    size_t           size;
    ActivationType   atype;
    MatrixXf        input;
    activation       activate_f;
    activation_grad  activate_grad_f;
    bool is_analyzer=false;
    bool useAutoDiff=false;
    
    
    std::function<MatrixXf (Layer *)> analyzer_func;

    
    void autodiff(MatrixXf & m)
    {        
        using namespace boost::math::differentiation;
        constexpr unsigned Order = 2;                  // Highest order derivative to be calculated.
        for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            // y = f(x)
            auto const y = make_fvar<double, Order>(m(i,j));  // Find derivatives at x=2.                
            m(i,j) = y.derivative(1);    
        }                        
    } 

    Layer(LayerType t, size_t s, ActivationType a) {
        type = t;
        size = s;
        atype = a;
        switch(a) {
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
            case FULLWAVE_RELU:  activate_f = fullwave_relu;
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
        input = MatrixXf(1,size);
    }
    ~Layer() {

    }
    bool isAnalyzer() const { return is_analyzer; }

    void setAnalyzer(std::function< MatrixXf (Layer *)>  func)
    {
        is_analyzer = true;
        analyzer_func = func;
    }
    MatrixXf runAnalyzer() {
        return analyzer_func(this);
    }
    void Activate(MatrixXf& tmp) {        
        input = tmp;
        activate_f(input);                
        //tanh(input);
    }
    void Grad(MatrixXf & tmp) {
        if(!useAutoDiff) activate_grad_f(tmp);        
        else autodiff(tmp);

        // this always works for all sigmoid functions
        // it works with fullwave sigmoid and relu
        //tanh_grad(tmp);    
    }

};


struct Connection {
    Layer * from,
          * to;


    MatrixXf weights;
    MatrixXf bias;

    
    Connection(Layer * from, Layer * to) {
        this->from = from;
        this->to   = to;
        weights = matrix_create(from->size,to->size);
        bias    = matrix_create(1,to->size);
        bias.fill(1.0f);

        BoxMuller bm;
        for(size_t i = 0; i < weights.rows(); i++)
            for(size_t j = 0; j < weights.cols(); j++)
                weights(i,j) = bm.generate()/std::sqrt(weights.rows());
        
    }
    ~Connection() {

    }

};

std::string stringify(MatrixXf & m)
{
    std::stringstream ss;
    ss << m;
    return ss.str();
}

enum {
    GD_OPTIMIZER,
    ADAM_OPTIMIZER,
    RMSPROP_OPTIMIZER,
};

struct ParameterSet {
    MatrixXf data;
    MatrixXf classes;
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

    ParameterSet( MatrixXf &d, MatrixXf &c,
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
    MatrixXf example;
    MatrixXf training;

    Batch(MatrixXf & e, MatrixXf & c) {        
        MatrixXf x = e;
        MatrixXf y = c;
        example    = x;
        training   = y;
    }
    Batch(const Batch & b) {
        MatrixXf x = b.example;
        MatrixXf y = b.training;
        example = x;
        training= y;
    }
    Batch& operator = (const Batch & b) {
        MatrixXf x = b.example;
        MatrixXf y = b.training;
        example = x;
        training= y;
        return *this;
    }
};

MatrixXf addToEachRow(const MatrixXf& m, const MatrixXf & v)
{
    MatrixXf r(m);
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
    std::vector<MatrixXf> errori;
    std::vector<MatrixXf*> weights;
    std::vector<MatrixXf*> bias;
    std::vector<MatrixXf> dWi;
    std::vector<MatrixXf> dbi;
    std::vector<MatrixXf> sdw;
    std::vector<MatrixXf> sdb;
    std::vector<MatrixXf> vdw;
    std::vector<MatrixXf> vdb;
    std::vector<MatrixXf> regi;
    std::vector<MatrixXf> wTi;
    std::vector<MatrixXf> errorLastTi;
    std::vector<MatrixXf> fprimei;
    std::vector<MatrixXf> inputTi;
    std::vector<MatrixXf> dWi_avg;
    std::vector<MatrixXf> dbi_avg;
    std::vector<MatrixXf> dWi_last;
    std::vector<MatrixXf> dbi_last;
    
    floatType loss = 1e6;
    floatType loss_widget=1e-6;

    

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
                ln = new Layer(INPUT,num_features, LINEAR);
            else if(i == num_layers-1)
                ln = new Layer(OUTPUT,num_outputs,output_activation);
            else
                ln = new Layer(HIDDEN, hidden[i-1], activations[i-1]);
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
        for(size_t i = 0; i < layers.size(); i++)
            delete layers[i];
        for(size_t i = 0; i < connections.size(); i++)
            delete connections[i];
    }

    size_t NumLayers() const { return layers.size(); }
    size_t NumConnections() const { return connections.size(); }
    size_t NumInputs() const { return num_features; }
    size_t NumOutputs() const { return num_outputs; }
    size_t LastLayer() const { return layers.size()-1; }

    void ForwardPass(MatrixXf& input) {
        assert(input.cols() == layers[0]->input.cols());
        layers[0]->input = input;
        MatrixXf tmp,tmp2;        
        for(size_t i = 0; i < connections.size(); i++)
        {         
            if(layers[i]->isAnalyzer()) {
                tmp = layers[i]->runAnalyzer();                
            }
            else {
                tmp  = layers[i]->input*connections[i]->weights;            
            }
            //tmp2 = tmp.addToEachRow(connections[i]->bias);            
            tmp2 = addToEachRow(tmp,connections[i]->bias);
            connections[i]->to->Activate(tmp2);       
        }
    }
    floatType CrossEntropyLoss(MatrixXf& prediction, MatrixXf& actual, floatType rs) {
        floatType total_err = 0;
        floatType reg_err = 0;
        total_err = Octave::sum(actual * Octave::log(prediction))(0,0);
        if(rs > 0)
            for(size_t i = 0; i < connections.size(); i++)
            {
                MatrixXf & weights = connections[i]->weights;
                reg_err += Octave::sum(hadamard(weights,weights))(0,0);
            }        
        return (-1.0 / actual.rows()*total_err) + rs*0.5*reg_err;
    }
    floatType MeanSquaredError(MatrixXf& prediction, MatrixXf & actual, floatType rs) {
        floatType total_err = 0;
        floatType reg_err = 0;
        MatrixXf tmp = actual - prediction;
        
        total_err = Octave::sum(hadamard(tmp,tmp))(0,0);
        
        if(rs > 0)
            for(size_t i = 0; i < connections.size(); i++)
            {
                MatrixXf & w = connections[i]->weights;
                reg_err += Octave::sum(hadamard(w,w))(0,0);
            }
            
        return ((0.5 / actual.rows()) * total_err) + (rs*0.5*reg_err);
    }
    MatrixXf& GetInput() {
        return layers[0]->input;
    }
    MatrixXf& GetOutput() {
        return layers[LastLayer()]->input;
    }
    // legacy
    std::vector<int> predict() {
        Layer* output_layer = layers[layers.size()-1];
        std::vector<int> prediction;
        prediction.resize(output_layer->input.rows());
        MatrixXf & input = output_layer->input;
        for(size_t i = 0; i < input.rows(); i++) {
            int max = 0;
            for(size_t j = 0; j < input.cols(); j++) {
                if(input(i,j) > input(i,max)) max = j;
            }
            prediction[i] = max;
        }
        return prediction;
    }
    floatType accuracy(MatrixXf & data, MatrixXf & classes) {
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
        for(size_t i = 0; i < batch_list.size(); i++) {
              std::random_shuffle(batch_list[i].begin(),batch_list[i].end());
              //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
              //std::shuffle (batch_list[i].begin(),batch_list[i].end(), std::default_random_engine(seed));
        }
    }
    
    void generate_batches(size_t num_batches,
                            size_t batch_size,
                            MatrixXf & data,
                            MatrixXf & classes,
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
                MatrixXf e = data.row(rc);                
                MatrixXf c = classes.row(rc);                
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
                dWi_last[idx] = dWi_avg[idx];
                dbi_last[idx] = dbi_avg[idx];                
            }
            dWi_avg[idx].fill(0);
            dbi_avg[idx].fill(0);
        }
    }
    /*
    void adam_optimize(const ParameterSet & ps, size_t epoch, size_t rows)
    {   
        floatType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }
        
        floatType alpha = currentLearningRate / rows;
        
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
            
            MatrixXf mdw_corr = sdw[i] / (1 - pow(ps.gamma1,epoch+1));
            MatrixXf mdb_corr = sdb[i] / (1 - pow(ps.gamma1,epoch+1));

            MatrixXf vdw_corr = vdw[i] / (1 - pow(ps.gamma2,epoch+1));
            MatrixXf vdb_corr = vdb[i] / (1 - pow(ps.gamma2,epoch+1));

            
            ValueList v;
            v(0) = vdw_corr;
            v = interp.eval("cwisemax",v,1e-08);
            vdw_corr = v(0).matrix_value();
            v(1) = vdb_corr;
            v = interp.eval("cwisemax",v,1e-08);
            vdb_corr = v(0).matrix_value();
            
            MatrixXf m1 = (alpha / (Octave::sqrt(vdw_corr)));
            MatrixXf m2 = (alpha / (Octave::sqrt(vdb_corr)));
            
            connections[i]->weights = connections[i]->weights - hadamard(m1,mdw_corr);
            connections[i]->bias    = connections[i]->bias    - hadamard(m2,mdb_corr);    
            
            if(ps.momentum_factor > 0 || ps.regularization_strength > 0)
            {
                dWi_last[i] = dWi_avg[i];
                dbi_last[i] = dbi_avg[i];                                
            }      
            
            dWi_avg[i].fill(0);
            dbi_avg[i].fill(0);
        }
    }       
    void rmsprop_optimize(const ParameterSet & ps, size_t epoch, size_t rows)
    {                    

        floatType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }
        
        floatType alpha = currentLearningRate / rows;
        
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

            MatrixXf vdw_corr = vdw[i] / (1 - pow(ps.gamma1,epoch+1));
            MatrixXf vdb_corr = vdb[i] / (1 - pow(ps.gamma1,epoch+1));

            ValueList v;
            v(0) = vdw_corr;
            v = interp.eval("cwisemax",v,1e-08);
            vdw_corr = v(0).matrix_value();
            v(1) = vdb_corr;
            v = interp.eval("cwisemax",v,1e-08);
            vdb_corr = v(0).matrix_value();
            
            
            MatrixXf m1 = (alpha / (Octave::sqrt(vdw_corr)));
            MatrixXf m2 = (alpha / (Octave::sqrt(vdb_corr)));
            
            connections[i]->weights = connections[i]->weights - hadamard(m1,dWi_avg[i]);
            connections[i]->bias    = connections[i]->bias    - hadamard(m2,dbi_avg[i]);    
            
            if(ps.momentum_factor > 0 || ps.regularization_strength > 0)
            {
                dWi_last[i] = dWi_avg[i];
                dbi_last[i] = dbi_avg[i];                                
            }           
            dWi_avg[i].fill(0);
            dbi_avg[i].fill(0);
        }
    }                        
    */
    void learn(size_t batch, size_t cur_batch_size)
    {
        for(size_t training = 0; training < cur_batch_size; training++)
        {
            MatrixXf& example = batch_list[batch][training].example;
            MatrixXf& target  = batch_list[batch][training].training;                    
            ForwardPass(example);
            
            size_t layer = layers.size()-1;
            Layer* to = layers[layer];
            Connection* con = connections[layer-1];
            
            errori[layer] = to->input - target;                                         
            dWi[layer-1]  = con->from->input.transpose() * errori[layer];
            dbi[layer-1]  = errori[layer];
                                                    
            for(layer = layers.size()-2; layer > 0; layer--)
            {                                                     
                size_t hidden_layer = layer-1;
                to  = layers[layer];
                con = connections[layer-1];
                wTi[hidden_layer] = connections[layer]->weights.transpose();                            
                errorLastTi[hidden_layer] = errori[layer+1]*wTi[hidden_layer];
                fprimei[hidden_layer] = con->to->input;
                
                con->to->Grad(fprimei[hidden_layer]);
                
                errori[layer] = hadamard(errorLastTi[hidden_layer],fprimei[hidden_layer]);
                inputTi[hidden_layer] = con->from->input.transpose();                            
                dWi[hidden_layer] = inputTi[hidden_layer] * errori[layer];
                dbi[hidden_layer] = errori[layer];
                
            }                                                                
            for(size_t idx=0; idx < connections.size(); idx++) {
                dWi_avg[idx] = dWi[idx] + dWi_avg[idx];
                dbi_avg[idx] = dbi[idx] + dbi_avg[idx];                    
            }         
            
        }
    }
    void verbosity(const ParameterSet & ps, size_t epoch, MatrixXf & data, MatrixXf & classes) {
        if(ps.verbose == true) {
            if(epoch % ps.ticks == 0 || epoch <= 1) {
                ForwardPass(data);
                if(ps.loss_function == CROSS_ENTROPY_LOSS) {
                    printf("EPOCH: %ld loss is %f\n",epoch, CrossEntropyLoss(GetOutput(),classes,ps.regularization_strength));
                }
                else {
                    printf("EPOCH: %ld loss is %f\n",epoch, loss=MeanSquaredError(GetOutput(),classes,ps.regularization_strength));
                }
            }
        }
    }
    void create_matrix()
    {
        MatrixXf beforeOutputT = createMatrixZeros(layers[layers.size()-2]->size,1);
        

        for(size_t i = 0; i < connections.size(); i++) {
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
        errori.push_back(createMatrixZeros(1,layers[LastLayer()]->size));
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
        
        MatrixXf & data = ps.data;
        MatrixXf & classes = ps.classes;
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
            //if(ps.optimizer == ADAM_OPTIMIZER) adam_optimize(ps,epoch,data.rows());
            //else if(ps.optimizer == RMSPROP_OPTIMIZER) rmsprop_optimize(ps,epoch,data.rows());
            //else gd_optimize(ps,epoch,data.rows());
            gd_optimize(ps,epoch,data.rows());
            verbosity(ps,epoch,data,classes);
            epoch++;
        }


    }

    // the only difference is shuffle is always off
    void nonstochastic(ParameterSet & ps) {
        
        MatrixXf & data = ps.data;
        MatrixXf & classes = ps.classes;
        size_t num_batches = data.rows();        
        size_t epoch = 0;

        create_matrix();
        ps.batch_size=1;
        ps.shuffle=false;        
        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);
        while(epoch <= ps.max_iters && loss > loss_widget) {                    
            for(size_t batch=0; batch < num_batches; batch++) learn(batch,ps.batch_size);                            
            //if(ps.optimizer == ADAM_OPTIMIZER) adam_optimize(ps,epoch,data.rows());
            //else if(ps.optimizer == RMSPROP_OPTIMIZER) rmsprop_optimize(ps,epoch,data.rows());
            //else gd_optimize(ps,epoch,data.rows());
            gd_optimize(ps,epoch,data.rows());
            verbosity(ps,epoch,data,classes);
            epoch++;
        }
    }

    void train(ParameterSet & ps, int type) {
        if(type == NONSTOCHASTIC) nonstochastic(ps);
        else stochastic(ps);
    }
};

struct NeuralNetwork : public Network
{
    
    NeuralNetwork() = default;

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

void XOR2(ActivationType atype, floatType lt, floatType mf)
{
    std::vector<floatType> examples = {0,0,0,1,1,0,1,1};
    std::vector<floatType> training = {0,1,1,0};
    std::vector<floatType> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<floatType> training_bp = {-1,1,1,-1};

    MatrixXf e = matrix_new(4,2,examples);
    MatrixXf t = matrix_new(4,1,training);
        
    //std::vector<int64_t> hidden = {5};
    //std::vector<ActivationType> activations = {atype};
    NeuralNetwork net;
    Layer * input = new Layer(INPUT,2,LINEAR);
    Layer * hidden= new Layer(HIDDEN,5,atype);
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
    MatrixXf &outs = net.GetOutput();
    std::cout << outs << std::endl;
}    

void Macho(ActivationType atype, floatType lt, floatType mf)
{
    std::vector<floatType> examples;
    std::vector<floatType> training;
    size_t n=1024;
    int64_t sizey = 128;
    for(size_t i = 0; i < n; i++)
    {
        examples.push_back(2*((double)rand() / (double)RAND_MAX)-1);
        training.push_back(2*((double)rand() / (double)RAND_MAX)-1);
    }

    MatrixXf e = matrix_new(1,n,examples);
    MatrixXf t = matrix_new(1,n,training);
        
    std::vector<int64_t> hidden = {sizey};
    std::vector<ActivationType> activations = {atype};
    Network net(n,hidden,activations,n,atype);
    ParameterSet p(e,t,1000,1);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 0;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = false;
    p.optimizer = GD_OPTIMIZER;
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.train(p,NONSTOCHASTIC);
    
    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    MatrixXf &output = net.GetOutput();
    //std::cout << output << std::endl;
    int s =0;
    for(size_t i = 0; i < output.cols(); i++)
    {        
        if( fabs(output(0,i) - training[i]) < 1e-3) s++;
    }
    std::cout << s << std::endl;  
}

void XOR(ActivationType atype, floatType lt, floatType mf)
{
    std::vector<floatType> examples = {0,0,0,1,1,0,1,1};
    std::vector<floatType> training = {0,1,1,0};
    std::vector<floatType> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<floatType> training_bp = {-1,1,1,-1};

    MatrixXf e = matrix_new(4,2,examples);
    MatrixXf t = matrix_new(4,1,training);
        
    std::vector<int64_t> hidden = {4};
    std::vector<ActivationType> activations = {atype};
    Network net(2,hidden,activations,1,LINEAR);
    ParameterSet p(e,t,10000,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = true;
    p.optimizer = RMSPROP_OPTIMIZER;
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.train(p,STOCHASTIC);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    MatrixXf &output = net.GetOutput();
    std::cout << output << std::endl;
}

// Eigen::MatrixXf<Eigen::MatrixXf<float,Eigen::Dynamic,Eigen::Dynamic>,Eigen::Dynamic,Eigen::Dynamic> m(3,3);

int main(int argc, char * argv[]) {     
    Macho(TANH,0.001,0.95);
    //XOR(SIGMOID,0.1,0.9);
    //XOR(RELU,0.1,0.9); 
    //XOR(FULLWAVE_SIGMOID,0.01,0.9);
    //XOR(FULLWAVE_RELU,0.1,0.9);             
    //XOR(TANH,0.1,0.9);   
    //XOR(ATAN,0.1,0.9);         
    //XOR(BALLS,0.1,0.9);         
    //XOR(ALGEBRA,0.1,0.9);                                    
    //XOR(SINWAVE,0.001,0.9);                                
}
