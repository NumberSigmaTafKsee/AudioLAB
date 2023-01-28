//#define EIGEN_USE_MKL_ALL
#include "simple_armadillo.hpp"
#include <cfloat>
#include <random>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <functional>
#include <cassert>

using namespace arma;
using Matrix = Mat<double>;

enum LayerType {
    INPUT=1,
    OUTPUT=2,
    HIDDEN=3,
};


enum LossFunctionType {
    CROSS_ENTROPY_LOSS=1,
    MEAN_SQUARED_ERROR=2,
};

typedef void (*activation)(Matrix & m);
typedef void (*activation_grad)(Matrix & m);


inline Matrix matrix_new(size_t rows, size_t cols, std::vector<double> & data) {
    Matrix m(rows,cols);
    for(size_t i = 0; i < rows; i++)
        for(size_t j = 0; j < cols; j++)
            m(i,j) = data[i*cols + j];    
    return m;
}
inline Matrix matrix_create(size_t rows, size_t cols) {
    Matrix  m(rows,cols);
    m.zeros();
    return m;
}
inline Matrix createMatrixZeros(size_t rows, size_t cols) {
    return matrix_create(rows,cols);
}

// all the cool kids call it the hadamard product.
Matrix hadamard(Matrix & a, Matrix &b)
{
    //return Matrix(a.cwiseProduct(b));
    return a % b;
}

/// Neural Network functions


void linear(Matrix& input) {
    
}

void linear_grad(Matrix& input) {
    input.fill(1);
}   


void sigmoid(Matrix & m)
{        
    //Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> t = m.array();
    Matrix t = -1.0f*m;
    m = (1 / (1 + exp(t)));    
}


void sigmoid_grad(Matrix & m)
{    
    Matrix r(m);    
    m = (r % ( 1 - r ));    
}


void tanh(Matrix & m)
{
    m = arma::tanh(m);
}


void tanh_grad(Matrix & m)
{    
    Matrix t = m;
    m = (1 - (t*t));    
}


void relu(Matrix & m)
{        
    for(size_t i = 0; i < m.n_rows; i++)
        for(size_t j = 0; j < m.n_cols; j++)
        {
            if(m(i,j) < 0) m(i,j) = 0;
            if(m(i,j) > 1) m(i,j) = 1;
        }
}


void relu_grad(Matrix & m)
{    
    for(size_t i = 0; i < m.n_rows; i++)
        for(size_t j = 0; j < m.n_cols; j++)
        {
            double x = m(i,j);
            if(x > FLT_MIN) m(i,j) = 1;
            else m(i,j) = 0;
        }    
}


void softmax(Matrix & m)
{                
    int i;
    for (i = 0; i < m.n_rows; i++){
        double summed = 0;
        int j;
        for (j = 0; j < m.n_cols; j++){
            summed += std::exp(m(i,j));
        }
        for (j = 0; j < m.n_cols; j++){
            m(i,j) =  std::exp(m(i,j)) / summed;
        }
    }
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
    double random(double min=0.0f, double max=1.0f) {
        std::uniform_real_distribution<double> distribution(min,max);
        return distribution(generator);
    }
};

// random already exists somewhere.
double randr(double min = 0.0f, double max = 1.0f) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return distribution(generator);
}


struct BoxMuller {
    double z0,z1;
    bool  do_generate;
    RandomNumbers random;

    BoxMuller() {
        z0=z1=0.0;
        do_generate = false;
    }

    double generate() {
        double epsilon = FLT_MIN;
        double two_pi  = 2 * M_PI;
        do_generate = !do_generate;
        if(!do_generate) return z1;
        double u1 = random.random();
        double u2 = random.random();
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
};

struct Connection;

struct Layer {
    LayerType        type;
    size_t           size;
    ActivationType   atype;
    Matrix        input;
    activation       activate_f;
    activation_grad  activate_grad_f;
    bool is_analyzer=false;
    
    std::function<Matrix (Layer *)> analyzer_func;

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
            case TANH:  activate_f = tanh;
                        activate_grad_f = tanh_grad;
                        break;
            case SOFTMAX:
                        activate_f = softmax;
                        activate_grad_f = linear_grad;
                        break;
        }
        input = Matrix(1,size);
    }
    ~Layer() {

    }
    bool isAnalyzer() const { return is_analyzer; }

    void setAnalyzer(std::function< Matrix (Layer *)>  func)
    {
        is_analyzer = true;
        analyzer_func = func;
    }
    Matrix runAnalyzer() {
        return analyzer_func(this);
    }
    void Activate(Matrix& tmp) {        
        input = tmp.eval();
        activate_f(input);        
    }
    void Grad(Matrix & tmp) {
        activate_grad_f(tmp);        
    }

};


struct Connection {
    Layer * from,
          * to;


    Matrix weights;
    Matrix bias;

    
    Connection(Layer * from, Layer * to) {
        this->from = from;
        this->to   = to;
        weights = matrix_create(from->size,to->size);
        bias    = matrix_create(1,to->size);
        bias.fill(1.0f);

        BoxMuller bm;
        
        for(size_t i = 0; i < weights.n_rows; i++)
            for(size_t j = 0; j < weights.n_cols; j++)
                weights(i,j) = bm.generate()/std::sqrt(weights.n_cols);
        
    }
    ~Connection() {

    }

};

std::string stringify(Matrix & m)
{
    std::stringstream ss;
    ss << m;
    return ss.str();
}

struct ParameterSet {
    Matrix data;
    Matrix classes;
    LossFunctionType loss_function;
    size_t batch_size;
    double learning_rate;
    double search_time;
    double regularization_strength;
    double momentum_factor;
    size_t max_iters;
    bool shuffle;
    bool verbose;

    ParameterSet( Matrix &d, Matrix &c,
                 size_t epochs, size_t bs,
                 LossFunctionType loss=MEAN_SQUARED_ERROR,
                 double lr = 0.01, double st = 0.0,
                 double rs=0.0,double m=0.2, bool s=true, bool v=true) {
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
    for(size_t i = 0; i < m.n_rows; i++)
        for(size_t j = 0; j < m.n_cols;j++)
            r(i,j) += v(0,j);
    return r;
}

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

    void ForwardPass(Matrix& input) {
        //assert(input.n_cols == layers[0]->input.n_cols);
        layers[0]->input = input.eval();
        Matrix tmp,tmp2;        
        for(size_t i = 0; i < connections.size(); i++)
        {         
            if(layers[i]->isAnalyzer()) {
                tmp = layers[i]->runAnalyzer();                
            }
            else {
                tmp  = (layers[i]->input*connections[i]->weights);
            }
            //tmp2 = tmp.addToEachRow(connections[i]->bias);            
            tmp2 = addToEachRow(tmp,connections[i]->bias);
            connections[i]->to->Activate(tmp2);       
        }
    }
    double CrossEntropyLoss(Matrix& prediction, Matrix& actual, double rs) {
        double total_err = 0;
        double reg_err = 0;
        total_err = accu(actual % (log(prediction)));
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix & weights = connections[i]->weights;
            reg_err += accu(weights % weights);
        }
        return (-1.0f / actual.n_rows*total_err) + rs*0.5f*reg_err;
    }
    double MeanSquaredError(Matrix& prediction, Matrix & actual, double rs) {
        double total_err = 0;
        double reg_err = 0;
        Matrix tmp = actual - prediction;
        total_err = accu(tmp % tmp);
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix & w = connections[i]->weights;
            reg_err += accu(w%w);
        }
        return ((0.5f / actual.n_rows) * total_err) + (rs*0.5f*reg_err);
    }
    Matrix& GetInput() {
        return layers[0]->input;
    }
    Matrix& GetOutput() {
        return layers[LastLayer()]->input;
    }
    // legacy
    std::vector<int> predict() {
        Layer* output_layer = layers[layers.size()-1];
        std::vector<int> prediction;
        prediction.resize(output_layer->input.n_cols);
        Matrix & input = output_layer->input;
        for(size_t i = 0; i < input.n_rows; i++) {
            int max = 0;
            for(size_t j = 0; j < input.n_cols; j++) {
                if(input(i,j) > input(max,i)) max = j;
            }
            prediction[i] = max;
        }
        return prediction;
    }
    double accuracy(Matrix & data, Matrix & classes) {
        ForwardPass(data);
        std::vector<int> p = predict();
        double num_correct = 0;
        for(size_t i = 0; i < data.n_cols; i++) {
            if(classes(i,p[i]) == 1)
                num_correct++;
        }
        return 100*num_correct/classes.n_rows;
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
                            Matrix & data,
                            Matrix & classes,
                            bool shuffle) {
        size_t rc = 0;
        batch_list.clear();
        for(size_t i = 0; i < num_batches; i++) {
            std::vector<Batch> l;
            size_t cur_batch_size = batch_size;
            if(i == num_batches) {
                if( data.n_rows % batch_size != 0) {
                    cur_batch_size = data.n_rows % batch_size;
                }
            }
            for(size_t j = 0; j < cur_batch_size; j++) {
                Matrix e = data.row(rc);                
                Matrix c = classes.row(rc);                
                Batch b(e,c);                
                l.push_back(b);
                rc = rc + 1;
                rc = rc % data.n_rows;
            }
            batch_list.push_back(l);
        }
        //if(shuffle) shuffle_batches();
    }
    void optimize(ParameterSet & ps) {
        std::vector<Matrix> errori;
        std::vector<Matrix> dWi;
        std::vector<Matrix> dbi;
        std::vector<Matrix> regi;

        Matrix beforeOutputT = createMatrixZeros(1,layers[layers.size()-2]->size);
        Matrix & data = ps.data;
        Matrix & classes = ps.classes;

        for(size_t i = 0; i < connections.size(); i++) {
            errori.push_back(createMatrixZeros(layers[i]->size,1));
            dWi.push_back(createMatrixZeros(connections[i]->weights.n_rows,
                                            connections[i]->weights.n_cols));
            dbi.push_back(createMatrixZeros(1,connections[i]->bias.n_cols));
            regi.push_back(createMatrixZeros(connections[i]->weights.n_rows,
                                            connections[i]->weights.n_cols));
        }
        errori.push_back(createMatrixZeros(1,layers[LastLayer()]->size));
        size_t num_hidden = layers.size()-2;
        std::vector<Matrix> wTi;
        std::vector<Matrix> errorLastTi;
        std::vector<Matrix> fprimei;
        std::vector<Matrix> inputTi;
        for(size_t k = 0; k < num_hidden; k++)
        {
            wTi.push_back(createMatrixZeros(connections[k+1]->weights.n_cols,connections[k+1]->weights.n_rows));
            errorLastTi.push_back(createMatrixZeros(1,wTi[k].n_cols));
            fprimei.push_back(createMatrixZeros(1,connections[k]->to->size));
            inputTi.push_back(createMatrixZeros(connections[k]->from->size,1));
        }
        std::vector<Matrix> dWi_avg;
        std::vector<Matrix> dbi_avg;
        std::vector<Matrix> dWi_last;
        std::vector<Matrix> dbi_last;
        for(size_t i = 0; i < connections.size(); i++) {
            dWi_avg.push_back(createMatrixZeros(connections[i]->weights.n_rows,connections[i]->weights.n_cols));
            dbi_avg.push_back(createMatrixZeros(1,connections[i]->bias.n_cols));
            dWi_last.push_back(createMatrixZeros(connections[i]->weights.n_rows,connections[i]->weights.n_cols));
            dbi_last.push_back(createMatrixZeros(1,connections[i]->bias.n_cols));
        }
        size_t num_batches = data.n_rows / ps.batch_size;

        if(data.n_rows % ps.batch_size != 0) num_batches++;

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
                    if(data.n_rows % ps.batch_size != 0) {
                        cur_batch_size = data.n_rows % ps.batch_size;
                    }
                }
                for(size_t training = 0; training < cur_batch_size; training++)
                {
                    Matrix& example = batch_list[batch][training].example;
                    Matrix& target  = batch_list[batch][training].training;                    
                    ForwardPass(example);
                    
                    for(size_t layer = layers.size()-1; layer > 0; layer--)
                    {
                        Layer* to = layers[layer];
                        Connection* con = connections[layer-1];
                        
                        if(layer == layers.size()-1) {
                            errori[layer] = to->input - target;                             
                            beforeOutputT = con->from->input.t();                            
                            dWi[layer-1] = (beforeOutputT * errori[layer]);
                            dbi[layer-1] = errori[layer];
                        }                                                                                           
                        
                        else {                                                       
                            size_t hidden_layer = layer-1;
                            wTi[hidden_layer] = connections[layer]->weights.t();                            
                            errorLastTi[hidden_layer] = (errori[layer+1]*wTi[hidden_layer]);
                            
                            fprimei[hidden_layer] = con->to->input.eval();
                            
                            con->to->Grad(fprimei[hidden_layer]);                            
                            errori[layer] = hadamard(errorLastTi[hidden_layer],fprimei[hidden_layer]);
                            
                            inputTi[hidden_layer] = con->from->input.t();                            
                            dWi[hidden_layer] = (inputTi[hidden_layer] * errori[layer]);
                            dbi[hidden_layer] = errori[layer];                                                        
                        }                                                                         
                    }                                        
                    for(size_t idx=0; idx < connections.size(); idx++) {
                        dWi_avg[idx] = dWi[idx] + dWi_avg[idx];
                        dbi_avg[idx] = dbi[idx] + dbi_avg[idx];                    
                    } 
                                                                                                                
                }
                
            }
            double currentLearningRate = ps.learning_rate;
            if(ps.search_time != 0) {
                currentLearningRate = ps.learning_rate / (1.0f + (epoch / ps.search_time));
            }
            
            double clr = currentLearningRate / data.n_rows;
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
                dWi_avg[idx].zeros();
                dbi_avg[idx].zeros();
                regi[idx].zeros();
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

void XOR2(ActivationType atype, double lt, double mf)
{
    std::vector<double> examples = {0,0,0,1,1,0,1,1};
    std::vector<double> training = {0,1,1,0};
    std::vector<double> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<double> training_bp = {-1,1,1,-1};

    Matrix e = matrix_new(4,2,examples);
    Matrix t = matrix_new(4,1,training);
        
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
    p.regularization_strength = 0.0001;
    p.verbose = true;
    p.shuffle = true;
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.optimize(p);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix &outs = net.GetOutput();
    std::cout << outs << std::endl;
}    
void XOR(ActivationType atype, double lt, double mf)
{
    std::vector<double> examples = {0,0,0,1,1,0,1,1};
    std::vector<double> training = {0,1,1,0};
    std::vector<double> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<double> training_bp = {-1,1,1,-1};

    Matrix e = matrix_new(4,2,examples);
    Matrix t = matrix_new(4,1,training);
    
    std::vector<int64_t> hidden = {5};
    std::vector<ActivationType> activations = {atype};
    Network net(2,hidden,activations,1,atype);
    ParameterSet p(e,t,1000,4);
    p.learning_rate = lt;
    p.momentum_factor = mf;
    p.regularization_strength = 0;
    p.verbose = true;
    p.shuffle = true;
    
    //p.loss_function = CROSS_ENTROPY_LOSS;
    std::cout << "Cranium Online" << std::endl;
    net.optimize(p);

    std::cout << "Ready." << std::endl;    
    net.ForwardPass(e);
    Matrix &output = net.GetOutput();
    std::cout << output << std::endl;
    
}

int main(int argc, char * argv[]) {  
   XOR(SIGMOID,0.1,0.9);
   XOR(TANH,0.1,0.9);
   XOR(RELU,0.1,0.9);         
}
