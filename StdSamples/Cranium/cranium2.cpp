// it's actually a little slower using mkl.
//#define EIGEN_USE_MKL_ALL
#include "Eigen.h"
#include <cfloat>
#include <random>
#include <algorithm>
#include <chrono>

using namespace SimpleEigen;

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
    Matrix<float> m(data,rows,cols);        
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

// all the cool kids call it the hadamard product.
template<typename T>
Matrix<T> hadamard(Matrix<T> & a, Matrix<T> &b)
{
    return Matrix<T>(a.cwiseProduct(b.matrix));
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
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.matrix.array();
    t = -t;
    m.matrix = (1 / (1 + t.exp())).matrix();    
}

template<typename T>
void sigmoid_grad(Matrix<T> & m)
{    
    Matrix<T> r(m);
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.matrix.array();
    m.matrix = (t * ( 1 - t )).matrix();    
}

template<typename T>
void tanh(Matrix<T> & m)
{
    m.matrix = m.matrix.array().tanh().matrix();
}

template<typename T>
void tanh_grad(Matrix<T> & m)
{    
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.matrix.array();
    m.matrix = (1 - (t*t)).matrix();    
}

template<typename T>
void relu(Matrix<T> & m)
{    
    m.matrix = m.matrix.cwiseMax(0).eval();    
}

template<typename T>
void relu_grad(Matrix<T> & m)
{    
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
        {
            T x = m.matrix(i,j);
            if(x > FLT_MIN) m.matrix(i,j) = 1;
            else m.matrix(i,j) = 0;
        }    
}

template<typename T>
void softmax(Matrix<T> & m)
{                
    int i;
    for (i = 0; i < m.rows(); i++){
        float summed = 0;
        int j;
        for (j = 0; j < m.cols(); j++){
            summed += std::exp(m(i, j));
        }
        for (j = 0; j < m.cols(); j++){
            m(i, j) =  std::exp(m(i, j)) / summed;
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
};

struct Connection;
struct Layer;

struct Connection {
    Layer * from,
          * to;

    Matrix<float> weights;
    Matrix<float> bias;
    Matrix<float> wTi;
    Matrix<float> dWi;
    Matrix<float> dbi;
    Matrix<float> regi;
    Matrix<float> dWi_avg;
    Matrix<float> dbi_avg;
    Matrix<float> dWi_last;
    Matrix<float> dbi_last;

    Connection(Layer * from, Layer * to);
    ~Connection() = default;

    void integrate() {                
        dWi_avg = dWi_avg + dWi;                 
        dbi_avg = dbi_avg + dbi;                                                                                                      
    }
    void optimize(float clr, float rs, float mf) {                        
        dWi_avg = dWi_avg * clr;        
        dbi_avg = dbi_avg * clr;        
        regi = weights * rs;
        dWi_avg = regi + dWi_avg;        
        dWi_last = dWi_last * mf;        
        dbi_last = dbi_last * mf;                
        dWi_avg = (dWi_last + dWi_avg)*-1.0f;        
        dbi_avg = (dbi_last + dbi_avg)*-1.0f;
        weights = dWi_avg + weights;                            
        bias    = dbi_avg + bias;                                                             
        dWi_last = dWi_avg * -1.0f;        
        dbi_last = dbi_avg * -1.0f;                                                
        dWi_avg.zero();        
        dbi_avg.zero();        
        regi.zero();     
    }
};

struct Layer {
    LayerType        type;
    size_t           size;
    ActivationType   atype;
    Matrix<float>        input;
    activation       activate_f;
    activation_grad  activate_grad_f;

    Matrix<float> errori;    
    Matrix<float> errorLastTi;
    Matrix<float> fprimei;
    Matrix<float> inputTi;
    Matrix<float> wTi;    
    Matrix<float> beforeOutputT;    

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
            case SOFTMAX:
                        activate_f = softmax<float>;
                        activate_grad_f = linear_grad<float>;
                        break;
        }
        input = Matrix<float>(1,size);
    }
    ~Layer() = default;

    void train_output(Matrix<float>& target, Connection * cconnect, Connection* pconnect) {        
        errori = input - target;
        pconnect->dWi = pconnect->from->input.t()  * errori;
        pconnect->dbi = errori;        
    }
    void train_hidden(Layer * next, Layer * prev, Connection * cconnect, Connection * nconnect, Connection * pconnect) {                   
        prev->errorLastTi = next->errori * cconnect->weights.t();
        prev->fprimei = pconnect->to->input.eval();
        pconnect->to->Grad(prev->fprimei);
        errori = hadamard(prev->errorLastTi,prev->fprimei);
        pconnect->dWi = pconnect->from->input.t() * errori;
        pconnect->dbi = errori;     
    }

    void Activate(Matrix<float>& tmp) {
        input = tmp.eval();
        activate_f(input);        
    }
    void Grad(Matrix<float> & tmp) {
        activate_grad_f(tmp);        
    }

};


Connection::Connection(Layer * from, Layer * to) {
    this->from = from;
    this->to   = to;
    weights = matrix_create(from->size,to->size);
    bias    = matrix_create(1,to->size);
    bias.fill(1.0f);

    BoxMuller bm;
    for(size_t i = 0; i < weights.rows(); i++)
        for(size_t j = 0; j < weights.cols(); j++)
            weights.set(i,j,bm.generate()/std::sqrt(weights.rows()));
    
    dWi  = createMatrixZeros(weights.rows(),weights.cols());
    dbi  = createMatrixZeros(1,bias.cols());
    regi = createMatrixZeros(weights.rows(),weights.cols());
    dWi_avg = createMatrixZeros(weights.rows(),weights.cols());
    dbi_avg = createMatrixZeros(1,bias.cols());
    dWi_last = createMatrixZeros(weights.rows(),weights.cols());
    dbi_last = createMatrixZeros(1,bias.cols());            
}
std::string stringify(Matrix<float> & m)
{
    std::stringstream ss;
    ss << m.matrix;
    return ss.str();
}

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
            tmp  = layers[i]->input*connections[i]->weights;            
            tmp2 = tmp.addToEachRow(connections[i]->bias);            
            connections[i]->to->Activate(tmp2);       
        }
    }
    float CrossEntropyLoss(Matrix<float>& prediction, Matrix<float>& actual, float rs) {
        float total_err = 0;
        float reg_err = 0;
        total_err = (actual * log(prediction)).sum();
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix<float> & weights = connections[i]->weights;
            reg_err += hadamard<float>(weights,weights).sum();
        }
        return (-1.0f / actual.rows()*total_err) + rs*0.5f*reg_err;
    }
    float MeanSquaredError(Matrix<float>& prediction, Matrix<float> & actual, float rs) {
        float total_err = 0;
        float reg_err = 0;
        Matrix<float> tmp = actual - prediction;
        total_err = hadamard<float>(tmp,tmp).sum();
        for(size_t i = 0; i < connections.size(); i++)
        {
            Matrix<float> & w = connections[i]->weights;
            reg_err += hadamard<float>(w,w).sum();
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
                if(input.get(i,j) > input.get(i,max)) max = j;
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
        Matrix<float> & data = ps.data;
        Matrix<float> & classes = ps.classes;
        
        size_t num_batches = data.rows() / ps.batch_size;

        if(data.rows() % ps.batch_size != 0) num_batches++;

        size_t epoch = 0;

        generate_batches(num_batches, ps.batch_size, data, classes, ps.shuffle);

        while(epoch <= ps.max_iters) {
            if(ps.shuffle) {
                shuffle_batches();
            }
            epoch++;
            if(epoch >= ps.max_iters) break;
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
                        Connection* pconnect = connections[layer-1];
                        Connection* cconnect = connections[layer];                        
                        if(layer == layers.size()-1) {                                                                                 
                            to->train_output(target,cconnect,pconnect);
                        }                                                                                             
                        else {                                                                                                                  
                            Layer* next = layers[layer+1];
                            Layer* prev = layers[layer-1];
                            Connection* nconnect = connections[layer+1];                            
                            to->train_hidden(next,prev,cconnect,nconnect,pconnect);
                        }                                                                    
                    }                      
                    for(size_t idx=0; idx < connections.size(); idx++) {                                                                      
                        Connection* connect = connections[idx];
                        connect->integrate();
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
                Connection* connect = connections[idx];                    
                connect->optimize(clr,ps.regularization_strength,ps.momentum_factor);                    
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
    Network net(2,hidden,activations,1,LINEAR);
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
    Matrix<float> &output = net.GetOutput();
    output.print();
}



int main(int argc, char * argv[]) {
   XOR(SIGMOID,0.1,0.9);
   XOR(TANH,0.1,0.9);
   XOR(RELU,0.1,0.9);   
}
