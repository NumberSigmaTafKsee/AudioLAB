#include "viper.hpp"

#include <cfloat>
#include <random>
#include <algorithm>

using namespace Viper;


typedef float ValueType;

enum LayerType {
    INPUT=1,
    OUTPUT=2,
    HIDDEN=3,
};


enum LossFunctionType {
    CROSS_ENTROPY_LOSS=1,
    MEAN_SQUARED_ERROR=2,
};

typedef void (*activation)(Matrix<ValueType> & m);
typedef void (*activation_grad)(Matrix<ValueType> & m);

Cublas _cublas;
Cublas *cublas=&_cublas;

inline Matrix<ValueType> matrix_new(size_t rows, size_t cols, std::vector<ValueType> & data) {
    Matrix<ValueType> m(rows,cols,data);
    return m;
}
inline Matrix<ValueType> matrix_create(size_t rows, size_t cols) {
    Matrix<ValueType>  m(rows,cols);
    m.zero();
    return m;
}
inline Matrix<ValueType> createMatrixZeros(size_t rows, size_t cols) {
    return matrix_create(rows,cols);
}
inline void linear_act(Matrix<ValueType>& input) {
    
}
inline void linear_grad_act(Matrix<ValueType>& input) {
    input.fill(1.0f);
}
inline void softmax_act(Matrix<ValueType>& input) {
    input.softmax();    
}
inline void tanh_act(Matrix<ValueType>& input) {
    input.tanh();    
}
inline void tanh_grad_act(Matrix<ValueType>& input) {
    input.tanh_deriv();    
}
inline void sigmoid_act(Matrix<ValueType>& input) {            
    input.sigmoid();    
}
inline void sigmoid_grad_act(Matrix<ValueType>& input) {
    input.sigmoid_deriv();
}
inline void relu_act(Matrix<ValueType>& input) {    
    input.relu();    
}
inline void relu_grad_act(Matrix<ValueType>& input) {    
    input.relu_deriv();    
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
    ValueType random(ValueType min=0.0f, ValueType max=1.0f) {
        std::uniform_real_distribution<double> distribution(min,max);
        return distribution(generator);
    }
};


// random already exists somewhere.
ValueType randr(ValueType min = 0.0f, ValueType max = 1.0f) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);                
    return distribution(generator);
}


struct BoxMuller {
    ValueType z0,z1;
    bool  do_generate;
    
    BoxMuller() {
        z0=z1=0.0;
        do_generate = false;
    }

    ValueType generate() {
        ValueType epsilon = FLT_MIN;
        ValueType two_pi  = 2 * M_PI;
        do_generate = !do_generate;
        if(!do_generate) return z1;
        ValueType u1 = randr();
        ValueType u2 = randr();
        while(u1 <= epsilon) {
            u1 = randr();
            u2 = randr();
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
};

struct Connection;
struct Network;

struct Layer {
    LayerType        type;
    size_t           size;
    ActivationType   atype;
    Matrix<ValueType>        input;
    activation       activate_f;
    activation_grad  activate_grad_f;
    Network *  neural;

    Layer(LayerType t, size_t s, ActivationType a) {
        type = t;
        size = s;
        atype = a;
        switch(a) {
            case LINEAR: activate_f = linear_act;
                         activate_grad_f = linear_grad_act;
                         break;
            case SIGMOID: activate_f = sigmoid_act;
                          activate_grad_f = sigmoid_grad_act;
                          break;                          
            case RELU:  activate_f = relu_act;
                        activate_grad_f = relu_grad_act;
                        break;
            case TANH:  activate_f = tanh_act;
                        activate_grad_f = tanh_grad_act;
                        break;
            case SOFTMAX:
                        activate_f = softmax_act;
                        activate_grad_f = linear_grad_act;
                        break;
        }
        input = Matrix<ValueType>(1,size);
        input.zero();
    }
    ~Layer() {

    }    
    void Activate(Matrix<ValueType>& tmp) {     
        input = tmp;
        activate_f(input);                        
    }
    void Grad(Matrix<ValueType> & tmp) {
        activate_grad_f(tmp);        
    }

    void neural_operation(Connection * c,size_t h, size_t l);
};



struct Connection {
    Layer * from,
          * to;

    Matrix<ValueType> &weights;
    Matrix<ValueType> &bias;
    Matrix<ValueType> input;

    Connection(Layer * from, Layer * to, Matrix<ValueType> &w, Matrix<ValueType> & b) 
    : weights(w),bias(b)
    {
        this->from = from;
        this->to   = to;
        //weights = matrix_create(from->size,to->size);
        //bias    = matrix_create(1,to->size);
        bias.fill(1.0f);

        BoxMuller bm;
        for(size_t i = 0; i < weights.rows(); i++)
            for(size_t j = 0; j < weights.cols(); j++)
                weights.set(i,j,bm.generate()/std::sqrt(weights.rows()));
        
        weights.upload_device();
    }
    ~Connection() {

    }
    void print() {
        weights.print();
        bias.print();
    }
};


struct ParameterSet {
    Matrix<ValueType> data;
    Matrix<ValueType> classes;
    LossFunctionType loss_function;
    size_t batch_size;
    ValueType learning_rate;
    ValueType search_time;
    ValueType regularization_strength;
    ValueType momentum_factor;
    size_t max_iters;
    bool shuffle;
    bool verbose;
    bool turbo;
    size_t ticks;
    ValueType decay=1.0;
    ValueType gamma1=0.9,gamma2=0.995;

    ParameterSet( Matrix<ValueType> &d, Matrix<ValueType> &c, 
                 size_t epochs, size_t bs,       
                 LossFunctionType loss=MEAN_SQUARED_ERROR,
                 ValueType lr = 0.01, ValueType st = 0.0,
                 ValueType rs=0.0,ValueType m=0.2, bool s=true, bool v=true, size_t _ticks=250,bool turbo_mode=true) {
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
            turbo = turbo_mode;
            ticks = _ticks;
    }
};

struct Batch {
    Matrix<ValueType> example;
    Matrix<ValueType> training;
    Matrix<ValueType> reverse_example;
    Matrix<ValueType> reverse_training;

    Batch(Matrix<ValueType> & e, Matrix<ValueType> & c) {
        example    = e.eval();
        training   = c.eval();
        reverse_example = example.t();
        reverse_training = training.t();        
    }    
    Batch(const Batch & b) {
        example.resize(b.example.M,b.example.N);
        example.copy(b.example);
        training.resize(b.training.M,b.training.N);
        training.copy(b.training);
        reverse_example = example.t();
        reverse_training = training.t();        
    }
    Batch& operator = (const Batch & b) {
        example.resize(b.example.M,b.example.N);
        example.copy(b.example);
        training.resize(b.training.M,b.training.N);
        training.copy(b.training);
        reverse_example  = example.t();
        reverse_training = training.t();        
        return *this;
    }
    void print() {
        std::cout << "--------------------------\n";
        example.print();
        training.print();
    }
};
                    
void neural_operation(Matrix<ValueType> & errorLastTi, Matrix<ValueType>& errori, Matrix<ValueType> &weights, Matrix<ValueType> inputfrom, Matrix<ValueType> & fprimei,
                      activation_grad grad, Matrix<ValueType> inputto, Matrix<ValueType>& dWi, Matrix<ValueType> & dbi)
{                      
    //errorLastTi = errori * weights.t();            
    mul_matrix(errorLastTi,errori,weights.t());
    fprimei.copy(inputfrom);
    grad(fprimei);    
    //dbi = hadamard(errorLastTi,fprimei);    
    hadamard_fast(dbi,errorLastTi,fprimei);
    //dWi = inputto * errori;                                                        
    mul_matrix(dWi,inputto,dbi);
    //dbi.copy(errori);        
}

void neural_optimize(Matrix<ValueType>& dWi_last, Matrix<ValueType>& dbi_last, ValueType momentum_factor,
                     Matrix<ValueType>& dWi_avg, Matrix<ValueType>& dbi_avg, ValueType clr,
                     Matrix<ValueType>& weights, Matrix<ValueType>& bias,
                     Matrix<ValueType>& regi, ValueType regularization_strength)
{
    //dWi_last = dWi_last * momentum_factor;                                                                                            
    //dbi_last = dbi_last * momentum_factor;                                                            
    mul_const(dWi_last,momentum_factor,dWi_last);
    mul_const(dbi_last,momentum_factor,dbi_last);
    
    //dWi_avg = dWi_avg * clr;            
    //dbi_avg = dbi_avg * clr;              
    mul_const(dWi_avg,clr,dWi_avg);
    mul_const(dbi_avg,clr,dbi_avg);
        
    //regi    = weights * regularization_strength;                                                                                                                                                
    //dWi_avg = regi + dWi_avg;                       
    mul_const(weights,regularization_strength,regi);
    add_matrix(dWi_avg,regi,dWi_avg);                        

    //dWi_avg = (dWi_last + dWi_avg);
    //dbi_avg = (dbi_last + dbi_avg);
    add_matrix(dWi_avg,dWi_last,dWi_avg);
    add_matrix(dbi_avg,dbi_last,dbi_avg);

    //weights = (-dWi_avg) + weights;            
    //bias    = (-dbi_avg) + bias;                                                    
    //add_matrix(weights,-dWi_avg,weights);
    //add_matrix(bias,-dbi_avg,bias);
    sub_matrix(weights,weights,dWi_avg);
    sub_matrix(bias,bias,dbi_avg);
    
    dWi_last.copy(dWi_avg);
    dbi_last.copy(dbi_avg);
}


void adam_optimize( ValueType gamma1, ValueType gamma2, 
                    Matrix<ValueType> & dWi_avg, Matrix<ValueType> & dbi_avg,                    
                    Matrix<ValueType> & dWi_last, Matrix<ValueType> & dbi_last,                    
                    ParameterSet& ps, Matrix<ValueType> & regi,
                    Matrix<ValueType> & sdw, Matrix<ValueType> & sdb, Matrix<ValueType> & vdw,
                    Matrix<ValueType> & vdb, int epoch, ValueType alpha,
                    Matrix<ValueType> & weights, Matrix<ValueType> & bias)
{   
    if(ps.momentum_factor != 0) {            
        mul_const(dWi_last,ps.momentum_factor,dWi_last);
        mul_const(dbi_last,ps.momentum_factor,dbi_last);
    }
        
    mul_const(dWi_avg,alpha,dWi_avg);
    mul_const(dbi_avg,alpha,dbi_avg);

    if(ps.regularization_strength != 0) {
        mul_const(weights,ps.regularization_strength,regi);
        add_matrix(dWi_avg,regi,dWi_avg);                        
    }
 
    if(ps.momentum_factor != 0 or ps.regularization_strength != 0) {
        add_matrix(dWi_avg,dWi_last,dWi_avg);
        add_matrix(dbi_avg,dbi_last,dbi_avg);
    }

    
    sdw = gamma1 * sdw + (1-gamma1)*dWi_avg;
    sdb = gamma1 * sdb + (1-gamma1)*dbi_avg;
                
    vdw = gamma2 * vdw + (1-gamma2) * hadamard(dWi_avg,dWi_avg);  //*/pow(dWi_avg,(ValueType)2.0);    
    vdb = gamma2 * vdb + (1-gamma2) * hadamard(dbi_avg,dbi_avg);  //*/pow(dbi_avg,(ValueType)2.0);
    
    Matrix<ValueType> mdw_corr = sdw / (1 - pow(gamma1,epoch+1));
    Matrix<ValueType> mdb_corr = sdb / (1 - pow(gamma1,epoch+1));
    Matrix<ValueType> vdw_corr = vdw / (1 - pow(gamma2,epoch+1));
    Matrix<ValueType> vdb_corr = vdb / (1 - pow(gamma2,epoch+1));
        
    weights = weights - hadamard((((ValueType)1.0/alpha) * (sqrt(vdw_corr + (ValueType)1e-08))),mdw_corr);
    bias    = bias    - hadamard((((ValueType)1.0/alpha) * (sqrt(vdb_corr + (ValueType)1e-08))),mdb_corr);    
    
    if(ps.momentum_factor != 0 or ps.regularization_strength != 0)
    {
        dWi_last.copy(dWi_avg);
        dbi_last.copy(dbi_avg);
    }
}                        

void rmsprop_optimize( ValueType gamma,
                    Matrix<ValueType> & dWi_avg, Matrix<ValueType> & dbi_avg,                    
                    Matrix<ValueType> & dWi_last, Matrix<ValueType> & dbi_last,                    
                    ParameterSet& ps, Matrix<ValueType> & regi,
                    Matrix<ValueType> & vdw,Matrix<ValueType> & vdb, int epoch, ValueType alpha,
                    Matrix<ValueType> & weights, Matrix<ValueType> & bias)
{                    
    if(ps.momentum_factor != 0) {            
        mul_const(dWi_last,ps.momentum_factor,dWi_last);
        mul_const(dbi_last,ps.momentum_factor,dbi_last);
    }
        
    mul_const(dWi_avg,alpha,dWi_avg);
    mul_const(dbi_avg,alpha,dbi_avg);

    if(ps.regularization_strength != 0) {
        mul_const(weights,ps.regularization_strength,regi);
        add_matrix(dWi_avg,regi,dWi_avg);                        
    }
 
    if(ps.momentum_factor != 0 or ps.regularization_strength != 0) {
        add_matrix(dWi_avg,dWi_last,dWi_avg);
        add_matrix(dbi_avg,dbi_last,dbi_avg);
    }

    vdw = gamma * vdw + (1-gamma)*hadamard(dWi_avg,dWi_avg);
    vdb = gamma * vdb + (1-gamma)*hadamard(dbi_avg,dbi_avg);

    weights = weights - hadamard((alpha/sqrt(vdw + 1e-08)),dWi_avg);
    bias    = bias    - hadamard((alpha/sqrt(vdb + 1e-08)),dbi_avg);

    if(ps.momentum_factor != 0 or ps.regularization_strength != 0)
    {
        dWi_last.copy(dWi_avg);
        dbi_last.copy(dbi_avg);
    }
}                        

            
// need to get gpu reduction
ValueType sumsq(Matrix<ValueType> &tmp) {
    ValueType total = 0;
    tmp.download_host();
    for(size_t i = 0; i < tmp.size(); i++)
        total += tmp[i]*tmp[i];
    return total;
}

// process by layers (list)
// process by connections (graph)


struct Network {
    size_t num_features;
    size_t num_outputs;
    std::vector<Layer*> layers;
    std::vector<Connection*> connections;
    std::vector<std::vector<Batch>> batch_list;

    // all the matrices are in the neural as it is much faster like this.
    std::vector<Matrix<ValueType>> errori;
    std::vector<Matrix<ValueType>> weights;
    std::vector<Matrix<ValueType>> bias;
    std::vector<Matrix<ValueType>> dWi;
    std::vector<Matrix<ValueType>> dbi;
    std::vector<Matrix<ValueType>> sdw;
    std::vector<Matrix<ValueType>> sdb;
    std::vector<Matrix<ValueType>> vdw;
    std::vector<Matrix<ValueType>> vdb;
    std::vector<Matrix<ValueType>> regi;
    std::vector<Matrix<ValueType>> wTi;
    std::vector<Matrix<ValueType>> errorLastTi;
    std::vector<Matrix<ValueType>> fprimei;
    std::vector<Matrix<ValueType>> inputTi;
    std::vector<Matrix<ValueType>> dWi_avg;
    std::vector<Matrix<ValueType>> dbi_avg;
    std::vector<Matrix<ValueType>> dWi_last;
    std::vector<Matrix<ValueType>> dbi_last;

    float loss = 1e6;
    float loss_widget=1e-6;
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
            ln->neural = this;
            layers[i] = ln;
        }
        size_t num_connections = num_layers-1;
        weights.resize(num_connections);
        bias.resize(num_connections);
        for(size_t i = 0; i < num_connections; i++)
        {
            assert(layers[i] != NULL);
            assert(layers[i+1]!= NULL);
            weights[i] =  Matrix<ValueType>(layers[i]->size,layers[i+1]->size);
            bias[i] = Matrix<ValueType>(1,layers[i+1]->size);            
            Connection * c = new Connection(layers[i],layers[i+1],weights[i],bias[i]);
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

    void ForwardPass(Matrix<ValueType>& input) {
        assert(input.cols() == layers[0]->input.cols());
        layers[0]->input = input.eval();
        Matrix<ValueType> tmp;       
        #pragma unroll
        for(size_t i = 0; i < connections.size(); i++)
        {                      
            //tmp = layers[i]->input * connections[i]->weights;       
            mul_matrix(tmp,layers[i]->input,connections[i]->weights);       
            tmp.addToEachRow(connections[i]->bias);                                                
            connections[i]->to->Activate(tmp);            
        }        
    }
    ValueType CrossEntropyLoss(Matrix<ValueType>& prediction, Matrix<ValueType>& actual, ValueType rs) {
        ValueType total_err = 0;
        ValueType reg_err = 0;        
        total_err = (actual * log2(prediction)).sum();
        
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix<ValueType> & weights = connections[i]->weights;
            reg_err += (hadamard(weights,weights)).sum();
        }
        
        return (-1.0f / actual.rows()*total_err) + rs*0.5f*reg_err;
    }
    ValueType MeanSquaredError(Matrix<ValueType>& prediction, Matrix<ValueType> & actual, ValueType rs) {
        ValueType total_err = 0;
        ValueType reg_err = 0;
        Matrix<ValueType> tmp = actual - prediction;
        //total_err = hadamard(tmp,tmp).sum();
        hadamard_fast(tmp,tmp,tmp);
        total_err = tmp.sum();
        
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix<ValueType> & w = connections[i]->weights;
            reg_err += (hadamard(w,w)).sum();            
        }
        
        return ((0.5f / actual.rows()) * total_err) + (rs*0.5f*reg_err);
    }
    Matrix<ValueType>& GetInput() {
        return layers[0]->input;
    }
    Matrix<ValueType>& GetOutput() {
        return layers[LastLayer()]->input;
    }
    // legacy
    std::vector<int> predict() {
        Layer* output_layer = layers[layers.size()-1];
        std::vector<int> prediction;
        prediction.resize(output_layer->input.rows());
        Matrix<ValueType> & input = output_layer->input;        
        for(size_t i = 0; i < input.rows(); i++) {
            int max = 0;
            for(size_t j = 0; j < input.cols(); j++) {
                if(input.get(i,j) > input.get(i,max)) max = j;
            }
            prediction[i] = max;
        }
        return prediction;
    }
    ValueType accuracy(Matrix<ValueType> & data, Matrix<ValueType> & classes) {
        ForwardPass(data);
        std::vector<int> p = predict();
        ValueType num_correct = 0;
        for(size_t i = 0; i < data.rows(); i++) {
            if(classes.get(i,p[i]) == 1)
                num_correct++;
        }
        return 100*num_correct/classes.rows();
    }
    void shuffle_batches() {
        for(size_t i = 0; i < batch_list.size(); i++)
            std::random_shuffle(batch_list[i].begin(),batch_list[i].end());        
    }
    
    void generate_batches(size_t num_batches,
                            size_t batch_size,
                            Matrix<ValueType> & data,
                            Matrix<ValueType> & classes,
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
                Matrix<ValueType> e = data.get_row(rc);                
                Matrix<ValueType> c = classes.get_row(rc);                                                
                Batch b(e,c);                
                l.push_back(b);
                rc = rc + 1;
                rc = rc % data.rows();
            }
            batch_list.push_back(l);
        }        
        if(shuffle) shuffle_batches();
    }
    std::vector<Batch> generate_list( Matrix<ValueType> & data,
                        Matrix<ValueType> & classes,
                        bool shuffle) {
        size_t rc = 0;
        std::vector<Batch> batches;
        int num_batches = data.rows();
        
        for(size_t i = 0; i < num_batches; i++) {                        
            Matrix<ValueType> e = data.get_row(rc);                
            Matrix<ValueType> c = classes.get_row(rc);                                                
            Batch b(e,c);                
            batches.push_back(b);
            rc = rc + 1;
            rc = rc % data.rows();
        }         
         if(shuffle) {
             std::random_shuffle(batches.begin(),batches.end());
         }                        
         return batches;
    }
    void train_pattern(Batch & batch) {
        Matrix<ValueType>& example = batch.example;
        Matrix<ValueType>& target  = batch.training;                                        
        Matrix<ValueType> beforeOutputT;             
        ForwardPass(example);                                             

        int layer = layers.size()-1;        
        Layer* to = layers[layer];
        Connection* con = connections[layer-1];
        size_t hidden_layer = layer-1;                        
                
        sub_matrix(errori[layer],to->input,target);
        dbi[hidden_layer].copy(errori[layer]);        
        mul_matrix(dWi[hidden_layer],con->from->input.t(),errori[layer]);
                
        for( layer = layers.size()-2; layer > 0; layer--)   
        {                   
            hidden_layer = layer-1;                        
            to = layers[layer];
            con = connections[layer-1];    
            to->neural_operation(con,hidden_layer,layer);                    
        }               
        
        #pragma unroll
        for(size_t idx=0; idx < connections.size(); idx++) {                                                                       
            //dWi_avg[idx] = dWi[idx] + dWi_avg[idx];                                                            
            add_matrix(dWi_avg[idx],dWi_avg[idx],dWi[idx]);
            //dbi_avg[idx] = dbi[idx] + dbi_avg[idx];      
            add_matrix(dbi_avg[idx],dbi_avg[idx],dbi[idx]);                                       
        }                                                         
    }   
    void reverse_train_pattern(Batch & batch) {
        Matrix<ValueType>& example = batch.reverse_example;
        Matrix<ValueType>& target  = batch.reverse_training;                                        
        Matrix<ValueType> beforeOutputT;             
        ForwardPass(example);                                             

        int layer = layers.size()-1;        
        Layer* to = layers[layer];
        Connection* con = connections[layer-1];
        size_t hidden_layer = layer-1;                        
                
        sub_matrix(errori[layer],to->input,target);
        dbi[hidden_layer].copy(errori[layer]);        
        mul_matrix(dWi[hidden_layer],con->from->input.t(),errori[layer]);
                
        for( layer = layers.size()-2; layer > 0; layer--)   
        {                   
            hidden_layer = layer-1;                        
            to = layers[layer];
            con = connections[layer-1];    
            to->neural_operation(con,hidden_layer,layer);                    
        }               
        
        #pragma unroll
        for(size_t idx=0; idx < connections.size(); idx++) {                                                                       
            //dWi_avg[idx] = dWi[idx] + dWi_avg[idx];                                                            
            add_matrix(dWi_avg[idx],dWi_avg[idx],dWi[idx]);
            //dbi_avg[idx] = dbi[idx] + dbi_avg[idx];      
            add_matrix(dbi_avg[idx],dbi_avg[idx],dbi[idx]);                                       
        }                                                         
    }   
    void adam_update(ParameterSet& ps, size_t epoch) {
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        ValueType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }                
        ValueType clr = currentLearningRate / data.rows();        
        for(size_t idx = 0; idx < connections.size(); idx++)
        {                                       
            adam_optimize(ps.gamma1,ps.gamma2,dWi_avg[idx],dbi_avg[idx],                            
                            dWi_last[idx],dbi_last[idx],ps,regi[idx],                            
                            sdw[idx],sdb[idx],vdw[idx],vdb[idx],
                            epoch,clr,weights[idx],bias[idx]);                  
        }                               
    }

    void rmsprop_update(ParameterSet& ps, size_t epoch) {
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        ValueType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }                
        ValueType clr = currentLearningRate / data.rows();        
        for(size_t idx = 0; idx < connections.size(); idx++)
        {                                       
            rmsprop_optimize(ps.gamma1,dWi_avg[idx],dbi_avg[idx],                            
                            dWi_last[idx],dbi_last[idx],  
                            ps,regi[idx],                          
                            vdw[idx],vdb[idx],
                            epoch,clr,weights[idx],bias[idx]);            
        }                               
    }

                    

    void train_batch(size_t batch, size_t training) {        
        train_pattern(batch_list[batch][training]);
    }    
    
    void update(ParameterSet& ps, size_t epoch) {
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        ValueType currentLearningRate = ps.learning_rate;
        if(ps.search_time != 0) {
            currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
        }                
        ValueType clr = currentLearningRate / data.rows();        
        for(size_t idx = 0; idx < connections.size(); idx++)
        {                                       
            neural_optimize(dWi_last[idx], dbi_last[idx], ps.momentum_factor,
                            dWi_avg[idx],dbi_avg[idx],clr,connections[idx]->weights,
                            connections[idx]->bias,regi[idx],ps.regularization_strength);                                    
        }        
    }
    void report(size_t epoch, ParameterSet & ps) {
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;

        if(ps.verbose == true) {            
            ForwardPass(data);
            if(ps.loss_function == CROSS_ENTROPY_LOSS) {
                printf("EPOCH: %ld loss is %f\n",epoch, loss=CrossEntropyLoss(GetOutput(),classes,ps.regularization_strength));
            }
            else {
                printf("EPOCH: %ld loss is %f\n",epoch, loss=MeanSquaredError(GetOutput(),classes,ps.regularization_strength));
            }        
        }
    }
    void clear_matrix() {
        errori.clear();
        dWi.clear();
        dbi.clear();
        regi.clear();

        for(size_t i = 0; i < connections.size(); i++) {            
            errori.push_back(createMatrixZeros(1,layers[i]->size));

            dWi.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
            sdw.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));                                                        
            vdw.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));                                            

            dbi.push_back(createMatrixZeros(1,connections[i]->bias.cols()));

            sdb.push_back(createMatrixZeros(1,connections[i]->bias.cols()));            
            vdb.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            
            regi.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
        }
        errori.push_back(createMatrixZeros(1,layers[LastLayer()]->size));
        size_t num_hidden = layers.size()-2;        

        wTi.clear();
        errorLastTi.clear();
        fprimei.clear();
        inputTi.clear();
        for(size_t k = 0; k < num_hidden; k++)
        {
            wTi.push_back(createMatrixZeros(connections[k+1]->weights.cols(),connections[k+1]->weights.rows()));
            errorLastTi.push_back(createMatrixZeros(1,wTi[k].cols()));
            fprimei.push_back(createMatrixZeros(1,connections[k]->to->size));
            inputTi.push_back(createMatrixZeros(connections[k]->from->size,1));
        }
        dWi_avg.clear();
        dbi_avg.clear();
        dWi_last.clear();
        dbi_last.clear();

        for(size_t i = 0; i < connections.size(); i++) {
            dWi_avg.push_back(createMatrixZeros(connections[i]->weights.rows(),connections[i]->weights.cols()));            
            dbi_avg.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            dWi_last.push_back(createMatrixZeros(connections[i]->weights.rows(),connections[i]->weights.cols()));            
            dbi_last.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
        }    
    }
    void batch(ParameterSet & ps) {
        
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        
        clear_matrix();        
        size_t num_batches = data.rows() / ps.batch_size;

        if(data.rows() % ps.batch_size != 0) num_batches++;

        size_t epoch = 0;

        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);
        
        while(epoch <= ps.max_iters && loss > loss_widget) {
            if(ps.shuffle) {
                shuffle_batches();
            }                                               
            for(size_t batch = 0; batch < num_batches; batch++) {            
                size_t cur_batch_size = ps.batch_size;                
                
                if(batch == num_batches) {
                    if(data.rows() % ps.batch_size != 0) {
                        cur_batch_size = data.rows() % ps.batch_size;
                    }
                }                   
                #pragma unroll                             
                for(size_t training = 0; training < cur_batch_size; training++)
                {                                               
                    train_batch(batch,training);                                                                                                    
                }                                          
            }                     
            adam_update(ps,epoch);            
            epoch++;                
            if(epoch % ps.ticks == 0 || epoch <=1 ) report(epoch,ps);                                                                                         
        }
    }
    void reverse_batch(ParameterSet & ps) {
        
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        
        clear_matrix();        
        size_t num_batches = data.rows() / ps.batch_size;

        if(data.rows() % ps.batch_size != 0) num_batches++;

        size_t epoch = 0;

        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);
        
        while(epoch <= ps.max_iters && loss > loss_widget) {
            if(ps.shuffle) {
                shuffle_batches();
            }                                               
            for(size_t batch = 0; batch < num_batches; batch++) {            
                size_t cur_batch_size = ps.batch_size;                
                
                if(batch == num_batches) {
                    if(data.rows() % ps.batch_size != 0) {
                        cur_batch_size = data.rows() % ps.batch_size;
                    }
                }                   
                #pragma unroll                             
                for(size_t training = 0; training < cur_batch_size; training++)
                {                                               
                    train_batch(batch,training);                                                                                                    
                }                                          
            }                     
            adam_update(ps,epoch);            
            epoch++;                
            if(epoch % ps.ticks == 0 || epoch <=1 ) report(epoch,ps);                                                                                         
        }
    }
    void batch_all(ParameterSet & ps) {
        
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        
        clear_matrix();        
        size_t num_batches = data.rows() / ps.batch_size;

        if(data.rows() % ps.batch_size != 0) num_batches++;

        size_t epoch = 0;

        generate_batches(num_batches, ps.batch_size, data, classes,false);
        
        while(epoch <= ps.max_iters && loss > loss_widget) {
            //if(ps.shuffle) {
            //    shuffle_batches();
            //}                                               
            for(size_t batch = 0; batch < num_batches; batch++) {            
                size_t cur_batch_size = ps.batch_size;                
                
                if(batch == num_batches) {
                    if(data.rows() % ps.batch_size != 0) {
                        cur_batch_size = data.rows() % ps.batch_size;
                    }
                }                   
                #pragma unroll                             
                for(size_t training = 0; training < cur_batch_size; training++)
                {                                               
                    train_batch(batch,training);                                                                                                    
                }                                          
            }                     
            update(ps,epoch);            
            epoch++;                
            if(epoch % ps.ticks == 0 || epoch <=1 ) report(epoch,ps);                                                                                         
        }
    }
    void train(ParameterSet & ps) {
        
        Matrix<ValueType> & data = ps.data;
        Matrix<ValueType> & classes = ps.classes;
        clear_matrix();
        
        size_t num_batches = data.rows() / ps.batch_size;

        if(data.rows() % ps.batch_size != 0) num_batches++;

        size_t epoch = 0;

        //generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);
        //std::vector<Batch> batches = generate_list(data,classes,ps.shuffle);
        size_t batch = 0;

        while(epoch <= ps.max_iters && loss > loss_widget) {                        
            Batch b(data,classes);
            train_pattern(b);                                       
            adam_update(ps,epoch); 
            ps.momentum_factor *= ps.decay;
            epoch++;            
            if(epoch % ps.ticks == 0 || epoch == 1) report(epoch,ps);                                          
        }
    }
};

void Layer::neural_operation(Connection * con, size_t hidden_layer, size_t layer)
{
    ::neural_operation( neural->errorLastTi[hidden_layer],
                      neural->errori[layer+1],
                      neural->connections[layer]->weights,            
                      con->to->input.eval(),                            
                      neural->fprimei[hidden_layer],
                      con->to->activate_grad_f,
                      con->from->input.t(),
                      neural->dWi[hidden_layer],
                      neural->dbi[hidden_layer]);
}

void XOR(ActivationType atype, ValueType lt, ValueType mf)
{

    std::vector<ValueType> examples = {0,0,0,1,1,0,1,1};
    std::vector<ValueType> training = {0,1,1,0};
    std::vector<ValueType> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<ValueType> training_bp = {-1,1,1,-1};

    Matrix<ValueType> e = matrix_new(4,2,examples);
    Matrix<ValueType> t = matrix_new(4,1,training);
    
    std::vector<int64_t> hidden = {16};
    std::vector<ActivationType> activations = {atype};
    Network net(2,hidden,activations,1,LINEAR);
    int nX = 1000;    
    ParameterSet p(e,t,nX,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 1e-6;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = false;
    p.ticks=10;
    net.loss_widget = 1e-9;
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.batch(p);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix<ValueType> &output = net.GetOutput();
    output.print();
    //clear_memory();
}

void RunXORTest() {    
    //for(size_t i = 0; i < 10; i++)
    XOR(TANH,0.1,0.0001);        
    
}

template<typename T>
std::ostream& operator << (std::ostream & o, const complex<T>& c)
{
    o << "(" << c.real() << "," << c.imag() << ")" << std::endl;
    return o;
}
template<typename T>
std::ostream& operator << (std::ostream & o, const ComplexVector<T> & c)
{
    for(size_t i = 0; i < c.size(); i++)
        o << "(" << c.real() << "," << c.imag() << "),";
    return o;
}

void RunTest(ActivationType atype, ValueType lt, ValueType mf)
{
    std::vector<ValueType> examples;
    std::vector<ValueType> training;
    size_t n = 128;
    int64_t sizey = 64;
    for(size_t i = 0; i < n; i++)
    {
        examples.push_back(2*((ValueType)rand()/(ValueType)RAND_MAX)-1);
        training.push_back(2*((ValueType)rand()/(ValueType)RAND_MAX)-1);
    }
    Matrix<ValueType> e = matrix_new(1,n,examples);
    Matrix<ValueType> t = matrix_new(1,n,training);
    std::vector<int64_t> hidden = {sizey,sizey,sizey};
    std::vector<ActivationType> activations = {atype,atype,atype};
    Network net(n,hidden,activations,n,atype);
    int nX = 1000;    
    ParameterSet p(e,t,nX,1);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    // this never does anything useful and is removed
    p.regularization_strength = 0;     
    p.search_time = 1000;
    p.verbose = true;
    p.shuffle = false;    
    p.ticks=1;
    // if you change this it may blow up or get worse
    p.decay=1;
    net.loss_widget=1e-6;
    std::cout << "Cranium Online" << std::endl;
    net.train(p);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix<ValueType> &output = net.GetOutput();
    output.download_host();
    //output.print();
    int s =0;
    for(size_t i = 0; i < output.cols(); i++)
    {        
        if( fabs(output(0,i) - training[i]) < 1e-3) s++;
    }
    std::cout << s << std::endl;        
}

/*
Matrix<Matrix<ValueType>> a(3,3);
    for(size_t i = 0; i < 3; i++) 
        for(size_t j = 0; j < 3; j++)
        {
            a(i,j) = Matrix<ValueType>(3,3);
            a(i,j).fill(1);
        }
    a.print();
*/    
int main(int argc, char * argv[]) {       
    //srand(time(NULL));
    //RunTest(TANH,0.01,0.0);              
    //RunXORTest();                    
    Matrix<float> a(1,1),b(3,3),c;
    a.fill(1);
    b.fill(2);
    c = a * b;
    c.print();
}