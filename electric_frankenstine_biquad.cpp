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
    std::vector<floatType> frequencyResponse() {
        const int ir_size = order;
        AudioDSP::FFTPlanRealDouble fftPlan(ir_size);
        std::vector<std::complex<floatType>> fr(ir_size);
        std::vector<floatType> r(ir_size);        
        r[0] = bq.Tick(1.0);
        for(size_t i = 1; i < order; i++)
            r[i] = bq.Tick(0);        
        AudioDSP::fft(fftPlan,r.data(),fr.data());
        fr[0] = fr[ir_size-1] = 0;
        for(size_t i = 1; i < ir_size-1; i++) r[i] = std::abs(fr[i]);        
        return r;
    }
    Matrix Loss(Matrix & m, Matrix & target) {                
        Matrix& input = net->layers[0]->input;                                
        Matrix x;
        getOutput(x);        
        loss = x-target;
        return loss;        
    }
    void PrintEquation(Matrix & w) {        
        Matrix& w1 = net->connections[net->connections.size()-7]->weights;
        Matrix& w2 = net->connections[net->connections.size()-6]->weights;
        Matrix& w3 = net->connections[net->connections.size()-5]->weights;
        Matrix& w4 = net->connections[net->connections.size()-4]->weights;
        Matrix& w5 = net->connections[net->connections.size()-3]->weights;
        Matrix& w6 = net->connections[net->connections.size()-2]->weights;

        std::cout << "z0=" << w1(0,0) << std::endl;
        std::cout << "z1=" << w2(1,0) << std::endl;
        std::cout << "z2=" << w3(2,0) << std::endl;
        std::cout << "p0=" << w4(3,0) << std::endl;
        std::cout << "p1=" << w5(4,0) << std::endl;
        std::cout << "p2=" << w6(5,0) << std::endl;
    }
    Matrix& getOutput() {
        Matrix& w1 = net->connections[net->connections.size()-7]->weights;
        Matrix& w2 = net->connections[net->connections.size()-6]->weights;
        Matrix& w3 = net->connections[net->connections.size()-5]->weights;
        Matrix& w4 = net->connections[net->connections.size()-4]->weights;
        Matrix& w5 = net->connections[net->connections.size()-3]->weights;
        Matrix& w6 = net->connections[net->connections.size()-2]->weights;

        Matrix& w = net->connections[net->connections.size()-1]->weights;
        bq.clear();
        bq.z[0] = w1(0,0);
        bq.z[1] = w2(1,0);
        bq.z[2] = w3(2,0);
        bq.p[0] = w4(3,0);
        bq.p[1] = w5(4,0);
        bq.p[2] = w6(5,0);
        std::vector<floatType> r = frequencyResponse();
        for(size_t i = 0; i < order; i++) output(0,i) = r[i];
        return output;
    }    
};

void Biquad() {
    floatType y = 0.25;
    std::vector<floatType> examples;
    std::vector<floatType> training;
    std::vector<std::complex<floatType>> fr;
    ButterworthLowpassFilter lp(2,44100.0);
    int ir_size=128;
    training = examples = lp.impulse_response(ir_size);
    AudioDSP::FFTPlanRealDouble fftPlan(ir_size);
    fr.resize(ir_size);
    AudioDSP::fft(fftPlan,examples.data(),fr.data());
    fr[0] = fr[ir_size-1] = 0;
    for(size_t i = 1; i < ir_size-1; i++) examples[i] = std::abs(fr[i]);
    training = examples;
    p.plot_x(examples.data(),examples.size(),"butterworth");
    Matrix e = matrix_new(1,ir_size,examples);
    Matrix t = matrix_new(1,ir_size,training);
    
    NeuralNetwork net;
    Layer * input = new Layer(INPUT,ir_size,LINEAR);    
    Layer * hidden1= new Layer(HIDDEN,ir_size,ATAN);
    Layer * hidden2= new Layer(HIDDEN,ir_size,ATAN);
    Layer * hidden3= new Layer(HIDDEN,ir_size,ATAN);
    Layer * hidden4= new Layer(HIDDEN,1,ATAN);
    Layer * hidden4= new Layer(HIDDEN,1,ATAN);
    Layer * hidden5= new Layer(HIDDEN,1,ATAN);
    Layer * hidden6= new Layer(HIDDEN,1,ATAN);
    Layer * hidden7= new Layer(HIDDEN,1,ATAN);
    Layer * hidden8= new Layer(HIDDEN,1,ATAN);
    BiquadLayer * output= new BiquadLayer(&net,ir_size,LINEAR);
    net.addLayer(input);    
    net.addLayer(hidden1);
    net.addLayer(hidden2);
    net.addLayer(hidden3);    
    net.addLayer(hidden4);
    net.addLayer(hidden5);
    net.addLayer(hidden6);
    net.addLayer(hidden7);
    net.addLayer(hidden8);
    net.addLayer(output);    
    net.connect();
    net.loss_widget = 1e-6;

    ParameterSet p(e,t,1000,1);
    p.batch_size=1;
    p.learning_rate = 0.001;
    p.momentum_factor = 0.001;
    p.regularization_strength = 0;
    p.search_time=1000;
    p.verbose = true;
    p.shuffle = false;
    p.optimizer = GD_OPTIMIZER;    
    
    std::cout << "Electric-Frankenstein 9000" << std::endl;
    net.train(p,NONSTOCHASTIC);
    net.ForwardPass(e);        
    output->Print();
}
int main(int argc, char * argv[]) {     
    p.setstyle("lines");   
    Biquad();      
    sleep(15);
}



