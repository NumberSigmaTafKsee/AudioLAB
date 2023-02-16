
#pragma once 

#include <cmath>

namespace Analog::Filters::MoogFilters
{


    template<typename T>
    struct TMoogFilter : public FilterProcessor
    {    
        TMoogFilter();
        ~TMoogFilter();

        void init();
        void calc();
        T process(T x);
                
        T getCutoff();
        void setCutoff(T c);
        T getRes();
        void setRes(T r);
    
        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
        };
        void setPort(int port, DspFloatType v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                setCutoff(v);
                break;
            case PORT_RESONANCE:
                setRes(v);
                break;		
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*process(I);
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out);
        
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DspFloatType * out) {
            ProcessSIMD(n,out,out);
        }
        
        T cutoff;
        T res;
        T fs;
        T y1,y2,y3,y4;
        T oldx;
        T oldy1,oldy2,oldy3;
        T x;
        T r;
        T p;
        T k;
    };

    template<typename T>
    TMoogFilter<T>::TMoogFilter() : FilterProcessor()
    {
        fs=44100.0;
        init();
    }

    template<typename T>
    TMoogFilter<T>::~TMoogFilter()
    {
    }

    template<typename T>
    void TMoogFilter<T>::init()
    {
        // initialize values
        y1=y2=y3=y4=oldx=oldy1=oldy2=oldy3=0;
        calc();
    };

    template<typename T>
    void TMoogFilter<T>::calc()
    {
        T f = (cutoff+cutoff) / fs; //[0 - 1]
        p=f*(1.8f-0.8f*f);
        k=p+p-1.f;

        T t=(1.f-p)*1.386249f;
        T t2=12.f+t*t;
        r = res*(t2+6.f*t)/(t2-6.f*t);
    };

    template<typename T>
    T TMoogFilter<T>::process(T input)
    {
        // process input
        x = bpsigmoid(input) - r*y4;

        //T q = p;
        //p = p+p;
        //Four cascaded onepole filters (bilinear transform)
        y1= x*p +  oldx*p - k*y1;    
        y2=y1*p + oldy1*p - k*y2;    
        y3=y2*p + oldy2*p - k*y3;    
        y4=y3*p + oldy3*p - k*y4;
        //p = q;
        //Clipper band limited sigmoid
        y4-=((y4*y4*y4)/6.f);

        oldx = x; oldy1 = y1; oldy2 = y2; oldy3 = y3;
        return y4;
    }

    template<typename T>
    void TMoogFilter<T>::ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
    {
        Undenormal denormal;
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            // process input
            x = in[i] - r*y4;

            //T q = p;
            //p = p+p;
            //Four cascaded onepole filters (bilinear transform)
            y1= x*p +  oldx*p - k*y1;    
            y2=y1*p + oldy1*p - k*y2;    
            y3=y2*p + oldy2*p - k*y3;    
            y4=y3*p + oldy3*p - k*y4;
            //p = q;
            //Clipper band limited sigmoid
            y4-=((y4*y4*y4)/6.f);

            oldx = x; oldy1 = y1; oldy2 = y2; oldy3 = y3;
            out[i] = y4;
        }
    }

    template<typename T>
    T TMoogFilter<T>::getCutoff(){ return cutoff; }

    template<typename T>
    void TMoogFilter<T>::setCutoff(T c){ cutoff=c; calc(); }

    template<typename T>
    T TMoogFilter<T>::getRes(){ return res; }

    template<typename T>
    void TMoogFilter<T>::setRes(T r) { res=r; calc(); }


///////////////////////////////////////////////////////////////////////////////////////////
// MusicDSP Kaka
///////////////////////////////////////////////////////////////////////////////////////////
struct AnotherMoogFilter1 : public FilterProcessor
{
//Init
// cutoff = cutoff freq in Hz
//fs = sampling frequency //(e.g. 44100Hz)
//res = resonance [0 - 1] //(minimum - maximum)

    DspFloatType f,fs,k,p,scale,r,y1,y2,y3,y4,oldx,oldy1,oldy2,oldy3;
    DspFloatType cutoff,Q;
    DspFloatType x;

    AnotherMoogFilter1(DspFloatType sampleRate, DspFloatType cutoff, DspFloatType resonance) : FilterProcessor() {
                
        coefficients(sampleRate,cutoff,resonance);
        x=y1=y2=y3=y4=oldx=oldy1=oldy2=oldy3=0;
    }

    void coefficients(DspFloatType sampleRate,DspFloatType frequency, DspFloatType resonance) 
    {
        fs = sampleRate;
        cutoff = frequency;
        Q = resonance;

        f = 2 * cutoff / fs; //[0 - 1]
        k = 3.6*f - 1.6*f*f -1; //(Empirical tuning)
        p = (k+1)*0.5;

        // resonance sucks 
        scale = std::exp((1-p)*1.386249);
        r = resonance*scale;        
        //DspFloatType t=(1.f-p)*1.386249f;
        //DspFloatType t2=12.f+t*t;
        //r = Q*(t2+6.f*t)/(t2-6.f*t);
    }
    void setCutoff(DspFloatType c) {        
        c = clamp(c,0,fs/2);
        coefficients(fs,c,Q);
    }
    void setResonance(DspFloatType res) {
        res = clamp(res,0,1);
        coefficients(fs,cutoff,res);
    }
	enum {
		PORT_CUTOFF,
		PORT_RESONANCE,
	};
	void setPort(int port, DspFloatType v) {
		switch(port)
		{
			case PORT_CUTOFF: setCutoff(v); break;
			case PORT_RESONANCE: setResonance(v); break;
		}
	}
    DspFloatType Tick(DspFloatType input, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
    {
        Undenormal denormal;
        DspFloatType c = cutoff;
        DspFloatType res = Q;
        coefficients(fs,c + 0.5*X*c,Q + 0.5*Y*Q);
        //Loop
        //--Inverted feed back for corner peaking
        x = input - r*y4;                
        
        //Four cascaded onepole filters (bilinear transform)
        y1=x*p + oldx*p - k*y1;        
        y2=y1*p+oldy1*p - k*y2;        
        y3=y2*p+oldy2*p - k*y3;        
        y4=y3*p+oldy3*p - k*y4;        

        coefficients(fs,c,res);

        //Clipper band limited sigmoid
        y4 = y4 - (y4*y4*y4)/6;        
        oldx  = x;
        oldy1 = y1;
        oldy2 = y2;
        oldy3 = y3;

        return A*y4;
    }
	void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
		Undenormal denormal;
		#pragma omp simd
		for(size_t i = 0; i < n; i++) {
			DspFloatType c = cutoff;
			DspFloatType res = Q;
			
			//Loop
			//--Inverted feed back for corner peaking
			x = in[i] - r*y4;                
			
			//Four cascaded onepole filters (bilinear transform)
			y1=x*p + oldx*p - k*y1;        
			y2=y1*p+oldy1*p - k*y2;        
			y3=y2*p+oldy2*p - k*y3;        
			y4=y3*p+oldy3*p - k*y4;        
			

			//Clipper band limited sigmoid
			y4 = y4 - (y4*y4*y4)/6;        
			oldx  = x;
			oldy1 = y1;
			oldy2 = y2;
			oldy3 = y3;
			out[i] = y3;
		}
	}
};

struct AnotherMoogFilter2 : public FilterProcessor
{
    // Moog 24 dB/oct resonant lowpass VCF
    // References: CSound source code, Stilson/Smith CCRMA paper.
    // Modified by paul.kellett@maxim.abel.co.uk July 2000

    DspFloatType f, p, q;             //filter coefficients
    DspFloatType b0, b1, b2, b3, b4;  //filter buffers (beware denormals!)
    DspFloatType t1, t2;              //temporary buffers
    DspFloatType fs,fc,res;

    // Set coefficients given frequency & resonance [0.0...1.0]
    AnotherMoogFilter2(DspFloatType sr, DspFloatType cutoff, DspFloatType r) : FilterProcessor()
    {
        fs = sr;
        fc = cutoff/sr;
        res = r;
        calc();
        b0=b1=b2=b3=b4=0;
    }
    void calc()
    {
        q = 1.0f - fc;
        p = fc + 0.8f * fc * q;
        f = p + p - 1.0f;
        q = res * (1.0f + 0.5f * q * (1.0f - q + 5.6f * q * q));
    }
    void setCutoff(DspFloatType f) { fc = f/fs; }
    void setResonance(DspFloatType r) { res = r; }
	enum {
		PORT_CUTOFF,
		PORT_RESONANCE,
	};
	void setPort(int port, DspFloatType v) {
		switch(port)
		{
			case PORT_CUTOFF: setCutoff(v); break;
			case PORT_RESONANCE: setResonance(v); break;
		}
	}
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X = 0, DspFloatType Y=0)
    {
        Undenormal denormals;
        calc();
        DspFloatType in = I - q*b4;       
        t1 = b1; //std::tanh(b1);  
        b1 = (in + b0) * p - b1 * f;
        t2 = b2; //std::tanh(b2);  
        b2 = (b1 + t1) * p - b2 * f;
        t1 = b3; //std::tanh(b3); 
        b3 = (b2 + t2) * p - b3 * f;
        b4 = (b3 + t1) * p - b4 * f;
        b4 = b4 - b4 * b4 * b4 * 0.166667f;
        b0 = in;
        return b4;
    }
	void ProcessSIMD(size_t n, DspFloatType * input, DspFloatType * out) {
		Undenormal denormal;
		#pragma omp simd
		for(size_t i = 0; i < n; i++) {			
			calc();
			DspFloatType in = input[i] - q*b4;       
			t1 = b1; //std::tanh(b1);  
			b1 = (in + b0) * p - b1 * f;
			t2 = b2; //std::tanh(b2);  
			b2 = (b1 + t1) * p - b2 * f;
			t1 = b3; //std::tanh(b3); 
			b3 = (b2 + t2) * p - b3 * f;
			b4 = (b3 + t1) * p - b4 * f;
			b4 = b4 - b4 * b4 * b4 * 0.166667f;
			b0 = in;
			out[i] = b4;
		}
	}
};	

struct MoogVCF2 : public FilterProcessor
{
    //Init
    DspFloatType fc;
    DspFloatType fs;
    DspFloatType res;
    DspFloatType out1,out2,out3,out4;
    DspFloatType in1,in2,in3,in4;
    
    MoogVCF2(DspFloatType sr, DspFloatType Fc, DspFloatType R) : FilterProcessor()
    {
        fs = sr;
        fc = Fc/sr;
        res= R;
        out1=out2=out3=out4=0;
        in1=in2=in3=in4=0;
    }
    void setCutoff(DspFloatType f) { fc = f/fs; }
    void setResonance(DspFloatType r) { res = r; }
	enum {
		PORT_CUTOFF,
		PORT_RESONANCE,
	};
	void setPort(int port, DspFloatType v) {
		switch(port)
		{
			case PORT_CUTOFF: setCutoff(v); break;
			case PORT_RESONANCE: setResonance(v); break;
		}
	}
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
        DspFloatType f = fc * 1.16;
        DspFloatType fb = res * (1.0 - 0.15 * f * f);
        DspFloatType input = I;
        input -= out4 * fb;
        input *= 0.35013 * (f*f)*(f*f);
        out1 = input + 0.3 * in1 + (1 - f) * out1; // Pole 1
        in1  = input;
        out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
        in2  = out1;
        out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
        in3  = out2;
        out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
        in4  = out3;
        return out4;
    }
	void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
		Undenormal denormal;
		#pragma omp simd
		for(size_t i = 0; i < n; i++) {			
			DspFloatType f = fc * 1.16;
			DspFloatType fb = res * (1.0 - 0.15 * f * f);
			DspFloatType input = in[i];
			input -= out4 * fb;
			input *= 0.35013 * (f*f)*(f*f);
			out1 = input + 0.3 * in1 + (1 - f) * out1; // Pole 1
			in1  = input;
			out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
			in2  = out1;
			out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
			in3  = out2;
			out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
			in4  = out3;
			out[i] =  out4;
		}
	}
};


constexpr DspFloatType gaintable[199] = { 0.999969, 0.990082, 0.980347, 0.970764, 0.961304, 0.951996, 0.94281, 0.933777, 0.924866, 0.916077, 0.90741, 0.898865, 0.890442, 0.882141 , 0.873962, 0.865906, 0.857941, 0.850067, 0.842346, 0.834686, 0.827148, 0.819733, 0.812378, 0.805145, 0.798004, 0.790955, 0.783997, 0.77713, 0.770355, 0.763672, 0.75708 , 0.75058, 0.744141, 0.737793, 0.731537, 0.725342, 0.719238, 0.713196, 0.707245, 0.701355, 0.695557, 0.689819, 0.684174, 0.678558, 0.673035, 0.667572, 0.66217, 0.65686, 0.651581, 0.646393, 0.641235, 0.636169, 0.631134, 0.62619, 0.621277, 0.616425, 0.611633, 0.606903, 0.602234, 0.597626, 0.593048, 0.588531, 0.584045, 0.579651, 0.575287 , 0.570953, 0.566681, 0.562469, 0.558289, 0.554169, 0.550079, 0.546051, 0.542053, 0.538116, 0.53421, 0.530334, 0.52652, 0.522736, 0.518982, 0.515289, 0.511627, 0.507996 , 0.504425, 0.500885, 0.497375, 0.493896, 0.490448, 0.487061, 0.483704, 0.480377, 0.477081, 0.473816, 0.470581, 0.467377, 0.464203, 0.46109, 0.457977, 0.454926, 0.451874, 0.448883, 0.445892, 0.442932, 0.440033, 0.437134, 0.434265, 0.431427, 0.428619, 0.425842, 0.423096, 0.42038, 0.417664, 0.415009, 0.412354, 0.409729, 0.407135, 0.404572, 0.402008, 0.399506, 0.397003, 0.394501, 0.392059, 0.389618, 0.387207, 0.384827, 0.382477, 0.380127, 0.377808, 0.375488, 0.37323, 0.370972, 0.368713, 0.366516, 0.364319, 0.362122, 0.359985, 0.357849, 0.355713, 0.353607, 0.351532,0.349457, 0.347412, 0.345398, 0.343384, 0.34137, 0.339417, 0.337463, 0.33551, 0.333588, 0.331665, 0.329773, 0.327911, 0.32605, 0.324188, 0.322357, 0.320557,0.318756, 0.316986, 0.315216, 0.313446, 0.311707, 0.309998, 0.308289, 0.30658, 0.304901, 0.303223, 0.301575, 0.299927, 0.298309, 0.296692, 0.295074, 0.293488, 0.291931, 0.290375, 0.288818, 0.287262, 0.285736, 0.284241, 0.282715, 0.28125, 0.279755, 0.27829, 0.276825, 0.275391, 0.273956, 0.272552, 0.271118, 0.269745, 0.268341, 0.266968, 0.265594, 0.264252, 0.262909, 0.261566, 0.260223, 0.258911, 0.257599, 0.256317, 0.255035, 0.25375 };

struct StilsonMoog2 : public FilterProcessor
{
    
    inline DspFloatType crossfade( DspFloatType amount, DspFloatType a, DspFloatType b ) {
        return (1-amount)*a + amount*b;
    }

    DspFloatType fc,fs,Q,p;
    DspFloatType cutoff,resonance;
    DspFloatType lowpass,highpass,bandpass,lastX;
    DspFloatType state[4], output; //should be gl obal scope / preserved between calls
    DspFloatType pre_gain,post_gain;

    StilsonMoog2(DspFloatType Fc, DspFloatType R, DspFloatType Fs) : FilterProcessor() {        
        fs = Fs;        
        cutoff = Fc;
        resonance = R;
        pre_gain = 2.0f;
        post_gain = 3.0f;
        setCutoff(Fc);
        setResonance(R);
        memset(&state[0],0,4*sizeof(DspFloatType));
        lowpass = highpass = bandpass = lastX = output = 0;
    }
    void setResonance(DspFloatType resonance)        
    {
        DspFloatType ix, ixfrac;
        int ixint;                
        ix = p * 99;
        ixint = floor( ix );        
        ixfrac = ix - ixint;        
        this->resonance = resonance;
        Q = resonance * crossfade( ixfrac, gaintable[ ixint + 99 ], gaintable[ ixint + 100 ] );        
    }
    void setCutoff(DspFloatType frequency) 
    {
        //code for setting pole coefficient based on frequency        
        cutoff = clamp(frequency,0,fs/2);
        fc = 2 * frequency / fs;
        if(fc < 0.005) fc = 0.005;
        DspFloatType x2 = fc*fc;
        DspFloatType x3 = fc*x2;
        p = -0.69346 * x3 - 0.59515 * x2 + 3.2937 * fc - 1.0072; //cubic fit by DFL, not 100% accurate but better than nothing...
    }        
	enum {
		PORT_CUTOFF,
		PORT_RESONANCE,
	};
	void setPort(int port, DspFloatType v) {
		switch(port)
		{
			case PORT_CUTOFF: setCutoff(v); break;
			case PORT_RESONANCE: setResonance(v); break;
		}
	}
    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
    {        
        int i,pole;
        DspFloatType temp, input;
        Undenormal denormal;
        
        input  = std::tanh(pre_gain*I);
        output = 0.25 * ( input - output ); //negative feedback
        output = clamp(output,-1,1);

        for( pole = 0; pole < 4; pole++) {
                temp = state[pole];
                output = output + p * (output - temp);
                state[pole] = output;
                output = output + temp;
                //if(std::fabs(output) < 1e-6) output=0;
        }        
        lowpass = output;
        highpass = input - output;
        bandpass = 3 * state[2] - lowpass; //got this one from paul kellet
        output = lowpass;
        output *= Q;  //scale the feedback
        lastX = I;
        return std::tanh(post_gain*output);
    }
	void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
		Undenormal denormal;
		int i,pole;
		DspFloatType temp, input;			
		#pragma omp simd
		for(size_t i = 0; i < n; i++) {						
			input  = std::tanh(pre_gain*in[i]);
			output = 0.25 * ( input - output ); //negative feedback
			output = clamp(output,-1,1);
			for( pole = 0; pole < 4; pole++) {
					temp = state[pole];
					output = output + p * (output - temp);
					state[pole] = output;
					output = output + temp;
					//if(std::fabs(output) < 1e-6) output=0;
			}        
			lowpass = output;
			highpass = input - output;
			bandpass = 3 * state[2] - lowpass; //got this one from paul kellet
			output = lowpass;
			output *= Q;  //scale the feedback
			lastX = input;
			out[i] = std::tanh(post_gain*output);
		}
	}
};

struct MoogLike : public FilterProcessor
{

    enum {
        LOWPASS,
        HIGHPASS    
    }   
    type = LOWPASS;

    DspFloatType coef[9];
    DspFloatType d[4];
    DspFloatType omega; //peak freq
    DspFloatType g;     //peak mag

    DspFloatType  fs,res;
    DspFloatType  _in,_out;

    // calculating coefficients:

    DspFloatType k,p,q,a;
    DspFloatType a0,a1,a2,a3,a4;
    
    MoogLike(DspFloatType Fs, DspFloatType Fc, DspFloatType Q, DspFloatType G) : FilterProcessor()
    {
        omega = Fc;
        q  = Q;
        fs = Fs;
        g  = G;
        k =p=q=a=a0=a1=a2=a3=a4=0;
    }

    void SetCoefficients(DspFloatType Fc, DspFloatType R)
    {
        omega = Fc;
        q     = R;
        k=(4.0*g-3.0)/(g+1.0);
        p=1.0-0.25*k;
        p*=p;
        
        if(type == LOWPASS) {
            // LP:
            a=1.0/(std::tan(0.5*omega)*(1.0+p));
            p=1.0+a;
            q=1.0-a;

            a0= 1.0/(k+p*p*p*p);
            a1= 4.0*(k+p*p*p*q);
            a2= 6.0*(k+p*p*q*q);
            a3= 4.0*(k+p*q*q*q);
            a4= (k+q*q*q*q);
            p = a0*(k+1.0);

            coef[0]=p;
            coef[1]=4.0*p;
            coef[2]=6.0*p;
            coef[3]=4.0*p;
            coef[4]=p;
            coef[5]=-a1*a0;
            coef[6]=-a2*a0;
            coef[7]=-a3*a0;
            coef[8]=-a4*a0;
        }
        else {
            // or HP:
            a=std::tan(0.5*omega)/(1.0+p);
            p=a+1.0;
            q=a-1.0;

            a0=1.0/(p*p*p*p+k);
            a1=4.0*(p*p*p*q-k);
            a2=6.0*(p*p*q*q+k);
            a3=4.0*(p*q*q*q-k);
            a4=    (q*q*q*q+k);
            p=a0*(k+1.0);

            coef[0]=p;
            coef[1]=-4.0*p;
            coef[2]=6.0*p;
            coef[3]=-4.0*p;
            coef[4]=p;
            coef[5]=-a1*a0;
            coef[6]=-a2*a0;
            coef[7]=-a3*a0;
            coef[8]=-a4*a0;
        }
    }
	enum {
		PORT_CUTOFF,
		PORT_RESONANCE,
		PORT_LPMODE,
		PORT_HPMODE,
	};
	void setPort(int port, DspFloatType v) {
		switch(port)
		{
			case PORT_CUTOFF: SetCoefficients(v,q); break;
			case PORT_RESONANCE: SetCoefficients(omega,v); break;
			case PORT_LPMODE: type = LOWPASS; break;
			case PORT_HPMODE: type = HIGHPASS; break;
		}
	}
    DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
    {
        Undenormal denormal;
                
        _in = I;
        // per sample:
        _out=coef[0]*_in+d[0];
        d[0]=coef[1]*_in+coef[5]*_out+d[1];
        d[1]=coef[2]*_in+coef[6]*_out+d[2];
        d[2]=coef[3]*_in+coef[7]*_out+d[3];
        d[3]=coef[4]*_in+coef[8]*_out;
        return _out;
    }
	void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
		Undenormal denormal;		
		#pragma omp simd
		for(size_t i = 0; i < n; i++) {						
			_in = in[i];
			// per sample:
			_out=coef[0]*_in+d[0];
			d[0]=coef[1]*_in+coef[5]*_out+d[1];
			d[1]=coef[2]*_in+coef[6]*_out+d[2];
			d[2]=coef[3]*_in+coef[7]*_out+d[3];
			d[3]=coef[4]*_in+coef[8]*_out;
			out[i] = _out;
		}
	}
};

}
