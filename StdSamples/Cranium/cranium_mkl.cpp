//#define EIGEN_USE_MKL_ALL
#include "carlo_mkl.hpp"
#include <cfloat>
#include <random>
#include <functional>
#include <algorithm>
#include <chrono>

using namespace Casino::MKL;

template<typename T> 
std::ostream& operator << (std::ostream & stream, const Matrix<float> & m) {
    for(size_t i = 0; i < m.M; i++)
    {
        for(size_t j = 0; j < m.N; j++)
        {
            stream << m(i,j) << ",";
        }
        stream << std::endl;
    }
    return stream;
}
enum LayerType {
    INPUT=1,
    OUTPUT=2,
    HIDDEN=3,
};


enum LossFunctionType {
    CROSS_ENTROPY_LOSS=1,
    MEAN_SQUARED_ERROR=2,
};

typedef void (*activation)(Matrix<float> & m);
typedef void (*activation_grad)(Matrix<float> & m);


inline Matrix<float> matrix_new(size_t rows, size_t cols, std::vector<float> & data) {
    Matrix<float> m(rows,cols);        
    for(size_t i = 0; i < rows; i++)
        for(size_t j = 0; j < cols; j++)
            m(i,j) = data[i*cols + j];
    return m;
}
inline Matrix<float> matrix_create(size_t rows, size_t cols) {
    Matrix<float>  m(rows,cols);
    m.zero();
    return m;
}
inline Matrix<float> createMatrixZeros(size_t rows, size_t cols) {
    return matrix_create(rows,cols);
}


/// Neural Network functions

template<typename T>
void linear(Matrix<T>& input) {
    
}
template<typename T>
void linear_grad(Matrix<T>& input) {
    input.fill(1);
}   

template<typename T>
void sigmoid(Matrix<T> & m)
{        
    Vector<T> x = m;
    Vector<T> r;
    r = (1.0f / (1.0f + exp(-x)));        
    memcpy(m.data(),r.data(),m.size()*sizeof(float));
}

template<typename T>
void sigmoid_grad(Matrix<T> & m)
{    
    Vector<T> x = m;
    Vector<T> r;
    r = (x * ( 1.0f - x));    
    memcpy(m.data(),r.data(),m.size()*sizeof(float));
}

template<typename T>
void tanh(Matrix<T> & m)
{    
    m = Casino::MKL::tanh(m);
}

template<typename T>
void tanh_grad(Matrix<T> & m)
{    
    Vector<T> x = m;
    Vector<T> r = (1.0f - (x*x));    
    memcpy(m.data(),r.data(),m.size()*sizeof(float));
}

template<typename T>
void relu(Matrix<T> & m)
{    
    for(size_t i = 0; i < m.M; i++)
        for(size_t j = 0; j < m.N; j++)
        {
            if(m(i,j) < 0) m(i,j) = 0;
            if(m(i,j) > 1) m(i,j) = 1;
        }
}

template<typename T>
void relu_grad(Matrix<T> & m)
{    
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            T x = m(i,j);
            if(x > FLT_MIN) m(i,j) = 1;
            else m(i,j) = 0;
        }    
}

template<typename T>
void atan(Matrix<T> & m)
{
    m = Casino::MKL::atan(m);
}

template<typename T>
void atan_grad(Matrix<T> & m)
{
    m = 1.0f / (1.0f+hadamard(m,m));
}

template<typename T>
void erf(Matrix<T> & m)
{
    m = Casino::MKL::erf(m);
}

template<typename T>
void erf_grad(Matrix<T> & m)
{
    m = (2.0f/(float)sqrt(M_PI))*Casino::MKL::exp(-hadamard(m,m));
}

template<typename T>
void clamp(Matrix<T> & m)
{
    for(size_t i = 0; i < m.M; i++)
        for(size_t j = 0; j < m.N; j++)
        {
            if(m(i,j) < T(-1)) m(i,j) = T(-1.0);
            if(m(i,j) > T(1)) m(i,j) = T(1.0);
        }
}

template<typename T>
void softmax(Matrix<T> & m)
{                
    T summed = sum(exp(m.eval()));
    m = exp(m.eval())/summed;    
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
    float random(float min=0.0f, float max=1.0f) {
        std::uniform_real_distribution<double> distribution(min,max);
        return distribution(generator);
    }
};

// random already exists somewhere.
float randr(float min = 0.0f, float max = 1.0f) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return distribution(generator);
}


struct BoxMuller {
    float z0,z1;
    bool  do_generate;
    RandomNumbers random;

    BoxMuller() {
        z0=z1=0.0;
        do_generate = false;
    }

    float generate() {
        float epsilon = FLT_MIN;
        float two_pi  = 2 * M_PI;
        do_generate = !do_generate;
        if(!do_generate) return z1;
        float u1 = random.random();
        float u2 = random.random();
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
    ERF=6,
};

struct Connection;

struct Layer {
    LayerType        type;
    size_t           size;
    ActivationType   atype;
    Matrix<float>    input;
    activation       activate_f;
    activation_grad  activate_grad_f;
    bool             is_analyzer=false;
    std::function<Matrix<float> (Layer *)> analyzer_func;

    Matrix<float> wTi;
    Matrix<float> errorLastTi;
    Matrix<float> fprimei;
    Matrix<float> inputTi;

    Layer(LayerType t, size_t s, ActivationType a) {
        type = t;
        size = s;
        atype = a;
        switch(a) {
            case LINEAR: activate_f = linear<float>;
                         activate_grad_f = linear_grad<float>;
                         break;
            case SIGMOID: activate_f = sigmoid<float>;
                          activate_grad_f = sigmoid_grad<float>;
                          break;                          
            case RELU:  activate_f = relu<float>;
                        activate_grad_f = relu_grad<float>;
                        break;
            case TANH:  activate_f = tanh<float>;
                        activate_grad_f = tanh_grad<float>;
                        break;
            case ATAN:  activate_f = atan<float>;
                        activate_grad_f = atan_grad<float>;
                        break;                        
            case ERF:  activate_f = erf<float>;
                        activate_grad_f = erf_grad<float>;
                        break;                                                
            case SOFTMAX:
                        activate_f = softmax<float>;
                        activate_grad_f = linear_grad<float>;
                        break;
        }
        input = Matrix<float>(1,size);
    }
    ~Layer() {

    }
    // evaluator X => encoder => parameters => evaluator => decoder => X, loss = evaluate(target - Y)
    // critic    X => encoder => parameters => critic => parameters => decoder => Y,loss = evaluate(target - Y)
    // analyzer  X => encoder => parameters => analyzer => output/decoder => Y,loss = target- output
    // function  X => encoder => parameters => function => output/decoder => loss = Y,target- output
    // algorithm X => encoder => parameters => algorithms => output/decoder => loss = Y,target- output

    bool isAnalyzer() const { return is_analyzer; }

    void setAnalyzer(std::function< Matrix<float> (Layer *)>  func)
    {
        is_analyzer = true;
        analyzer_func = func;
    }
    Matrix<float> runAnalyzer() {
        return analyzer_func(this);
    }
    void Activate(Matrix<float>& tmp) {
        input = tmp.eval();
        activate_f(input);        
    }
    void Grad(Matrix<float> & tmp) {
        activate_grad_f(tmp);        
    }

};


struct Connection {
    Layer * from,
          * to;


    Matrix<float> weights;
    Matrix<float> bias;

    Matrix<float> errori;
    Matrix<float> dWi;
    Matrix<float> dbi;
    Matrix<float> regi;
    Matrix<float> dWi_avg;
    Matrix<float> dbi_avg;
    Matrix<float> dWi_last;
    Matrix<float> dbi_last;

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


struct ParameterSet {
    Matrix<float> data;
    Matrix<float> classes;
    LossFunctionType loss_function;
    size_t batch_size;
    float learning_rate;
    float search_time;
    float regularization_strength;
    float momentum_factor;
    size_t max_iters;
    bool shuffle;
    bool verbose;

    ParameterSet( Matrix<float> &d, Matrix<float> &c,
                 size_t epochs, size_t bs,
                 LossFunctionType loss=MEAN_SQUARED_ERROR,
                 float lr = 0.01, float st = 0.0,
                 float rs=0.0,float m=0.2, bool s=true, bool v=true) {
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
    }
};

struct Batch {
    Matrix<float> example;
    Matrix<float> training;

    Batch(Matrix<float> & e, Matrix<float> & c) {
        example    = e;
        training   = c;
    }
    Batch(const Batch & b) {
        example = b.example;
        training= b.training;
    }
    Batch& operator = (const Batch & b) {
        example = b.example;
        training= b.training;
        return *this;
    }
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

    void ForwardPass(Matrix<float>& input) {
        assert(input.cols() == layers[0]->input.cols());
        layers[0]->input = input.eval();
        Matrix<float> tmp,tmp2;        
        for(size_t i = 0; i < connections.size(); i++)
        {                                
            if(layers[i]->isAnalyzer()) {
                tmp = layers[i]->runAnalyzer();                
            }
            else {
                tmp  = layers[i]->input*connections[i]->weights;            
            }
            tmp2 = tmp.addToEachRow(connections[i]->bias);            
            connections[i]->to->Activate(tmp2);       
        }
    }
    float CrossEntropyLoss(Matrix<float>& prediction, Matrix<float>& actual, float rs) {
        float total_err = 0;
        float reg_err = 0;
        Matrix<float> temp = log2(prediction);        
        total_err = sum(hadamard<float>(actual,temp));        
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix<float> & weights = connections[i]->weights;
            reg_err += sum(hadamard<float>(weights,weights));
        }
        return (-1.0f / actual.rows()*total_err) + rs*0.5f*reg_err;
    }
    float MeanSquaredError(Matrix<float>& prediction, Matrix<float> & actual, float rs) {
        float total_err = 0;
        float reg_err = 0;
        Matrix<float> tmp = actual - prediction;
        total_err = sum(hadamard<float>(tmp,tmp));
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix<float> & w = connections[i]->weights;
            reg_err += sum(hadamard<float>(w,w));
        }
        return ((0.5f / actual.rows()) * total_err) + (rs*0.5f*reg_err);
    }
    Matrix<float>& GetInput() {
        return layers[0]->input;
    }
    Matrix<float>& GetOutput() {
        return layers[LastLayer()]->input;
    }
    // legacy
    std::vector<int> predict() {
        Layer* output_layer = layers[layers.size()-1];
        std::vector<int> prediction;
        prediction.resize(output_layer->input.rows());
        Matrix<float> & input = output_layer->input;
        for(size_t i = 0; i < input.rows(); i++) {
            int max = 0;
            for(size_t j = 0; j < input.cols(); j++) {
                if(input(i,j) > input(i,max)) max = j;
            }
            prediction[i] = max;
        }
        return prediction;
    }
    float accuracy(Matrix<float> & data, Matrix<float> & classes) {
        ForwardPass(data);
        std::vector<int> p = predict();
        float num_correct = 0;
        for(size_t i = 0; i < data.rows(); i++) {
            if(classes(i,p[i]) == 1)
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
                            Matrix<float> & data,
                            Matrix<float> & classes,
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
                Matrix<float> e = data.row(rc);                
                Matrix<float> c = classes.row(rc);                
                Batch b(e,c);                
                l.push_back(b);
                rc = rc + 1;
                rc = rc % data.rows();
            }
            batch_list.push_back(l);
        }
        if(shuffle) shuffle_batches();
    }
    void optimize(ParameterSet & ps) {
        std::vector<Matrix<float>> errori;
        std::vector<Matrix<float>> dWi;
        std::vector<Matrix<float>> dbi;
        std::vector<Matrix<float>> regi;

        Matrix<float> beforeOutputT = createMatrixZeros(layers[layers.size()-2]->size,1);
        Matrix<float> & data = ps.data;
        Matrix<float> & classes = ps.classes;

        for(size_t i = 0; i < connections.size(); i++) {
            errori.push_back(createMatrixZeros(1,layers[i]->size));
            dWi.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
            dbi.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            regi.push_back(createMatrixZeros(connections[i]->weights.rows(),
                                            connections[i]->weights.cols()));
        }
        errori.push_back(createMatrixZeros(1,layers[LastLayer()]->size));
        size_t num_hidden = layers.size()-2;
        std::vector<Matrix<float>> wTi;
        std::vector<Matrix<float>> errorLastTi;
        std::vector<Matrix<float>> fprimei;
        std::vector<Matrix<float>> inputTi;
        for(size_t k = 0; k < num_hidden; k++)
        {
            wTi.push_back(createMatrixZeros(connections[k+1]->weights.cols(),connections[k+1]->weights.rows()));
            errorLastTi.push_back(createMatrixZeros(1,wTi[k].cols()));
            fprimei.push_back(createMatrixZeros(1,connections[k]->to->size));
            inputTi.push_back(createMatrixZeros(connections[k]->from->size,1));
        }
        std::vector<Matrix<float>> dWi_avg;
        std::vector<Matrix<float>> dbi_avg;
        std::vector<Matrix<float>> dWi_last;
        std::vector<Matrix<float>> dbi_last;
        for(size_t i = 0; i < connections.size(); i++) {
            dWi_avg.push_back(createMatrixZeros(connections[i]->weights.rows(),connections[i]->weights.cols()));
            dbi_avg.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
            dWi_last.push_back(createMatrixZeros(connections[i]->weights.rows(),connections[i]->weights.cols()));
            dbi_last.push_back(createMatrixZeros(1,connections[i]->bias.cols()));
        }
        size_t num_batches = data.rows() / ps.batch_size;

        if(data.rows() % ps.batch_size != 0) num_batches++;

        size_t epoch = 0;

        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);

        while(epoch <= ps.max_iters) {
            if(ps.shuffle) {
                shuffle_batches();
            }
            epoch++;
            
            for(size_t batch = 0; batch < num_batches; batch++) 
            {
                size_t cur_batch_size = ps.batch_size;
                
                if(batch == num_batches) {
                    if(data.rows() % ps.batch_size != 0) {
                        cur_batch_size = data.rows() % ps.batch_size;
                    }
                }
                for(size_t training = 0; training < cur_batch_size; training++)
                {
                    Matrix<float>& example = batch_list[batch][training].example;
                    Matrix<float>& target  = batch_list[batch][training].training;       
                    
                    ForwardPass(example);
                    
                    for(size_t layer = layers.size()-1; layer > 0; layer--)
                    {
                        Layer* to = layers[layer];
                        Connection* con = connections[layer-1];
                        if(layer == layers.size()-1) {
                            errori[layer] = to->input - target;                             
                            beforeOutputT = con->from->input.t();                            
                            dWi[layer-1] = beforeOutputT * errori[layer];
                            dbi[layer-1] = errori[layer].eval();                            
                        }                        
                        else {                            
                            size_t hidden_layer = layer-1;
                            wTi[hidden_layer] = connections[layer]->weights.t();                            
                            errorLastTi[hidden_layer] = errori[layer+1]*wTi[hidden_layer];
                            fprimei[hidden_layer] = con->to->input.eval();                            
                            con->to->Grad(fprimei[hidden_layer]);                                                                                                                
                            errori[layer] = hadamard(errorLastTi[hidden_layer],fprimei[hidden_layer]);                            
                            inputTi[hidden_layer] = con->from->input.t();                            
                            dWi[layer-1] = inputTi[hidden_layer] * errori[layer];
                            dbi[layer-1] = errori[layer].eval();                                                                                                                                            
                        }                        
                    }                                        
                    for(size_t idx=0; idx < connections.size(); idx++) {
                        dWi_avg[idx] = dWi[idx] + dWi_avg[idx];
                        dbi_avg[idx] = dbi[idx] + dbi_avg[idx];                                               
                    }                     
                }
            }
            float currentLearningRate = ps.learning_rate;
            if(ps.search_time != 0) {
                currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
            }
            
            float clr = currentLearningRate / data.rows();
            
            for(size_t idx = 0; idx < connections.size(); idx++)
            {
                dWi_avg[idx] = dWi_avg[idx] * clr;
                dbi_avg[idx] = dbi_avg[idx] * clr;                
                regi[idx] = connections[idx]->weights * ps.regularization_strength;
                dWi_avg[idx] = regi[idx] + dWi_avg[idx];
                dWi_last[idx] = dWi_last[idx] * ps.momentum_factor;
                dbi_last[idx] = dbi_last[idx] * ps.momentum_factor;
                dWi_avg[idx] = (dWi_last[idx] + dWi_avg[idx]);
                dbi_avg[idx] = (dbi_last[idx] + dbi_avg[idx]);                    
                dWi_avg[idx] = dWi_avg[idx] * -1.0f;
                dbi_avg[idx] = dbi_avg[idx] * -1.0f;                
                connections[idx]->weights = dWi_avg[idx] + connections[idx]->weights;
                connections[idx]->bias    = dbi_avg[idx] + connections[idx]->bias;                                                                    
                dWi_last[idx] = dWi_avg[idx]*-1.0f;
                dbi_last[idx] = dbi_avg[idx]*-1.0f;
                dWi_avg[idx].zero();
                dbi_avg[idx].zero();
                regi[idx].zero();                
            }
            
            
            if(ps.verbose == true) {
                if(epoch % 250 == 0 || epoch <= 1) {
                    ForwardPass(data);
                    if(ps.loss_function == CROSS_ENTROPY_LOSS) {
                        printf("EPOCH: %ld loss is %f\n",epoch, CrossEntropyLoss(GetOutput(),classes,ps.regularization_strength));
                    }
                    else {
                        printf("EPOCH: %ld loss is %f\n",epoch, MeanSquaredError(GetOutput(),classes,ps.regularization_strength));
                    }
                }
            }
        }
    }
};

/*
Autodiff    

Inputs/Outputs
    Digital
    Integer
    String
    Real
    Complex
    Vector
    Matrix

Layers
    FIR
    IIR
    Delay
    Recurrent
    Matrix
    Tree
    Map
    Graph
    Finite State Machine
    Logic
    If-Else
    Analyzer
    DSP
    Filter    
    Evaluator
    Critic
    Program
    Algorithm
    Function
    Equation
    
Genetic Programmer
    Dynamic Hidden
    Online Real-Time Learning
    Continual Learning
    Expand/Grow
    Shrink/Reduce
    Evolving Hidden Matrix
    Evolving Hidden Graph
*/



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

void XOR2(ActivationType atype, float lt, float mf)
{
    std::vector<float> examples = {0,0,0,1,1,0,1,1};
    std::vector<float> training = {0,1,1,0};
    std::vector<float> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<float> training_bp = {-1,1,1,-1};

    Matrix<float> e = matrix_new(4,2,examples);
    Matrix<float> t = matrix_new(4,1,training);
        
    //std::vector<int64_t> hidden = {5};
    //std::vector<ActivationType> activations = {atype};
    NeuralNetwork net;
    Layer * input = new Layer(INPUT,2,LINEAR);
    Layer * hidden= new Layer(HIDDEN,5,atype);
    Layer * output= new Layer(OUTPUT,1,atype);
    net.addLayer(input);
    net.addLayer(hidden);
    net.addLayer(output);
    net.connect();
    ParameterSet p(e,t,1000,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 0;
    p.verbose = true;
    p.shuffle = true;
    //it is for softmax
    //softmax is always nan
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.optimize(p);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix<float> &outs = net.GetOutput();
    outs.print();
}    

void XOR(ActivationType atype, float lt, float mf)
{

    std::vector<float> examples = {0,0,0,1,1,0,1,1};
    std::vector<float> training = {0,1,1,0};
    std::vector<float> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<float> training_bp = {-1,1,1,-1};

    Matrix<float> e = matrix_new(4,2,examples);
    Matrix<float> t = matrix_new(4,1,training);
        
    std::vector<int64_t> hidden = {16};
    std::vector<ActivationType> activations = {atype};
    Network net(2,hidden,activations,1,atype);
    ParameterSet p(e,t,1000,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 0;
    p.verbose = true;
    p.shuffle = true;
    p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.optimize(p);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix<float> &output = net.GetOutput();
    output.print();
}


int main(int argc, char * argv[]) {
    XOR2(SIGMOID,0.1,0.9);
    XOR2(TANH,0.1,0.9);
    XOR2(RELU,0.1,0.9);  
    XOR2(ATAN,0.1,0.9);  
    XOR2(ERF,0.1,0.9);      
   
   /*
   Matrix<float> a(3,3),b(3,3),c;
   a.fill(1);
   b.fill(2);
   c = a+b;
   c.print();
   */
   /*
   Vector<float> a(5),b(5),c;
   a.fill(1);
   b.fill(2);
   c = a + b;
   c.print();
   c = b - a;
   c.print();
   c = a * b;
   c.print();
   c = a / b;
   c.print();
   */
}
