// Gain = Y/X
// Transfer Function = Polynomial/Coefficients
// IIR = coefficients, 1/2 poles/zeros, biquads
// Delay = 1 + g*z^-N = N taps

#include "cranium.hpp"

struct Biquad
{
    DspFloatType z[3];
    DspFloatType p[3];
    DspFloatType x[2];
    DspFloatType y[2];

    Biquad() {
        x[0] = x[1] = 0;
        y[0] = y[1] = 0;    
    }
    void setCoeffs(DspFloatType Z[3], DspFloatType P[3]) {
        memcpy(z,Z,sizeof(z));
        memcpy(p,P,sizeof(p));
    }

    DspFloatType Tick(DspFloatType I)
    {
        DspFloatType r = I*z[0] + x[0]*z[1] + x[1] * z[2] - y[0]*p[0] - y[1]*p[1];
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
    DspFloatType sampleRate,frequency,q;    
    IIRFilters::BiquadFilterCascade filter;

    ButterworthLowpassFilter(int order, DspFloatType sr)
    {
        this->order = order;
        this->sampleRate = sr;        
        setFilter(1000.0,0.5)
    }
    void setFilter(DspFloatType f, DspFloatType Q) {
        auto x = IIRFilters::biquadlp(order,Q);
        frequency = f;
        q = Q;
        auto c = IIRFilters::AnalogBiquadCascade(x,f,sampleRate)
        filter.setCoefficients(c);
    }
    void setCutoff(DspFloatType f) {
        if(f <= 0) return;
        if(f >= sampleRate/2) return;
        setFilter(f,q);
    }
    void setQ(DspFloatType Q) {
        if(Q < 0.5) Q = 0.5;
        if(Q > 999) Q = 999;
        setFilter(frequency,Q);
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
        return filter.Tick(I,A,X,Y);
    }
    std::vector<DspFloatType> impulse_response(size_t n)
    {
        std::vector<DspFloatType> r(n);
        r[0] = Tick(1.0);
        for(size_t i = 1; i < n; i++) r[i] = Tick(0);
        return r;
    }
};


void XOR(ActivationType atype, DspFloatTypeType lt, DspFloatTypeType mf)
{
    std::vector<DspFloatTypeType> examples = {0,0,0,1,1,0,1,1};
    std::vector<DspFloatTypeType> training = {0,1,1,0};
    std::vector<DspFloatTypeType> examples_bp = {-1,-1,-1,1,1,-1,1,1};
    std::vector<DspFloatTypeType> training_bp = {-1,1,1,-1};

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


// Eigen::Matrix<Eigen::Matrix<DspFloatType,Eigen::Dynamic,Eigen::Dynamic>,Eigen::Dynamic,Eigen::Dynamic> m(3,3);

int main(int argc, char * argv[]) {     
    
}
