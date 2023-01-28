%module minidnn
%{
#include "MiniDNN.h"

%}

%include "stdint.i"
%include "std_vector.i"


%inline %{

enum ActivationType
{
    IDENTITY,
    RELU,
    MISH,
    SIGMOID,
    SOFTMAX,
    TANH
};

struct Layer 
{
    MiniDNN::Layer * layer;

    Layer() { layer = NULL; }
};

struct ConvolutionalLayer : public Layer 
{

    ConvolutionalLayer(ActivationType type, const int in_width, const int in_height,
                      const int in_channels, const int out_channels,
                      const int window_width, const int window_height) 
    {     
        switch(type)
        {
            case IDENTITY: layer = new MiniDNN::Convolutional<MiniDNN::Identity>(in_width, in_height, in_channels, out_channels, window_width, window_height);
                           break;
            case RELU: layer = new MiniDNN::Convolutional<MiniDNN::ReLU>(in_width, in_height, in_channels, out_channels, window_width, window_height);
                           break;
            case MISH: layer = new MiniDNN::Convolutional<MiniDNN::Mish>(in_width, in_height, in_channels, out_channels, window_width, window_height);
                           break;
            case SIGMOID: layer = new MiniDNN::Convolutional<MiniDNN::Sigmoid>(in_width, in_height, in_channels, out_channels, window_width, window_height);
                           break;
            case SOFTMAX: layer = new MiniDNN::Convolutional<MiniDNN::Softmax>(in_width, in_height, in_channels, out_channels, window_width, window_height);
                           break;
            case TANH: layer = new MiniDNN::Convolutional<MiniDNN::Tanh>(in_width, in_height, in_channels, out_channels, window_width, window_height);
                           break;
            default:
                assert(true == false);
        }
    }             
};

struct FullyConnectedLayer : public Layer
{
    FullyConnectedLayer(ActivationType type, const int in_size, const int out_size)
    {     
        switch(type)
        {
            case IDENTITY: layer = new MiniDNN::FullyConnected<MiniDNN::Identity>(in_size,out_size);
                           break;
            case RELU: layer = new MiniDNN::FullyConnected<MiniDNN::ReLU>(in_size,out_size);
                           break;
            case MISH: layer = new MiniDNN::FullyConnected<MiniDNN::Mish>(in_size,out_size);
                           break;
            case SIGMOID: layer = new MiniDNN::FullyConnected<MiniDNN::Sigmoid>(in_size,out_size);
                           break;
            case SOFTMAX: layer = new MiniDNN::FullyConnected<MiniDNN::Softmax>(in_size,out_size);
                           break;
            case TANH: layer = new MiniDNN::FullyConnected<MiniDNN::Tanh>(in_size,out_size);
                           break;
        }
    }
};

struct MaxpoolingLayer : public Layer 
{
    MaxpoolingLayer(ActivationType type, const int in_width, const int in_height, const int in_channels,
                   const int pooling_width, const int pooling_height)
    {        
        switch(type)
        {
            case IDENTITY: layer = new MiniDNN::MaxPooling<MiniDNN::Identity>(in_width,in_height,in_channels, pooling_width, pooling_height);            
                           break;
            case RELU: layer = new MiniDNN::MaxPooling<MiniDNN::ReLU>(in_width,in_height,in_channels, pooling_width, pooling_height);            
                           break;
            case MISH: layer = new MiniDNN::MaxPooling<MiniDNN::Mish>(in_width,in_height,in_channels, pooling_width, pooling_height);            
                           break;
            case SIGMOID: layer = new MiniDNN::MaxPooling<MiniDNN::Sigmoid>(in_width,in_height,in_channels, pooling_width, pooling_height);            
                           break;
            case SOFTMAX: layer = new MiniDNN::MaxPooling<MiniDNN::Softmax>(in_width,in_height,in_channels, pooling_width, pooling_height);            
                           break;
            case TANH: layer = new MiniDNN::MaxPooling<MiniDNN::Tanh>(in_width,in_height,in_channels, pooling_width, pooling_height);            
                           break;

        }
    }
};

struct Optimizer
{
    MiniDNN::Optimizer *opt;

    Optimizer()
    {
        opt = NULL;
    }
    ~Optimizer() 
    {
        if(opt) delete opt;
    }
    
};

struct AdaGradOptimizer : public Optimizer
{
    AdaGradOptimizer(double lrate=0.001, double eps=1e-6)
    {
        if(opt) delete opt;
        opt = new MiniDNN::AdaGrad(lrate,eps);
    }
};

struct AdamOptimizer : public Optimizer
{
    AdamOptimizer(double lrate=1e-3, double eps=1e-6, double beta1=0.9, double beta2=0.999)
    {
        if(opt) delete opt;
        opt = new MiniDNN::Adam(lrate,eps,beta1,beta2);        
    }
};

struct RMSPropOptimizer : public Optimizer 
{
    RMSPropOptimizer(double lrate=1e-3, double eps=1e-6, double gamma=0.0)
    {
        if(opt) delete opt;
        opt = new MiniDNN::RMSProp(lrate,eps,gamma);        
    }
};

struct SGDOptimizer : public Optimizer 
{
    SGDOptimizer(double lrate=0.001, double eps=1e-6)
    {
        if(opt) delete opt;
        opt = new MiniDNN::SGD(lrate,eps);
    }
};

struct Output
{
    MiniDNN::Output * out;
};

struct BinaryClassEntropyOutput : public Output
{
    BinaryClassEntropyOutput() { out = new MiniDNN::BinaryClassEntropy(); }
};

struct MultiClassEntropyOutput : public Output
{
    MultiClassEntropyOutput() { out = new MiniDNN::MultiClassEntropy(); }
};

struct RegressionMSEOutput : public Output 
{
    RegressionMSEOutput() { out = new MiniDNN::RegressionMSE(); }
};


struct Network
{
    MiniDNN::Network * net;
    MiniDNN::VerboseCallback verbosecb;

    Network(bool verbose = true)
    {
        net = new MiniDNN::Network();
        if(verbose) net->set_callback(verbosecb);
    }

    void add_layer(Layer * layer)
    {
        net->add_layer(layer->layer);
    }

    void set_output(Output & output)
    {
        net->set_output(output.out);
    }

    int num_layers() { return net->num_layers(); }
    
    void init(double mu=0.0, double sigma = 0.01)    
    {
        net->init(mu,sigma);
    }
    void fit(Optimizer * opt, const Eigen::MatrixXd &x, const Eigen::MatrixXd& y, int batch_size, int epoch, int seed=-1)
    {                
        net->fit(*opt->opt,x, y, batch_size, epoch, seed);        
    }
    
    Eigen::MatrixXd predict(const Eigen::MatrixXd& x)
    {
        return net->predict(x);
    }
    
    void export_net(const char * folder, const char * filename)
    {
        net->export_net(folder,filename);
    }
    void read_net(const char * folder, const char * filename)
    {
        net->read_net(folder,filename);
    }
};
%}

