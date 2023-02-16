///////////////////////////////////////////////////////////////////////////
// Amplifiers
// Clippers and Distortion
// Discrete Summation Formulas
// Wave Digital Filters
///////////////////////////////////////////////////////////////////////////

float pre_gain  = 1.0f;
float post_gain = 1.0f;

inline float amp_clamp(float x, float a, float b) {
    return x < a? a: x > b? b : x;
}
inline float preamp(float x) {
    return pre_gain * x;
}
inline float postamp(float x) {
    return post_gain * x;
}

inline float tanh_normal(float x, float K=10,float r = 1.0f) {
    return post_gain*std::tanh(pre_gain*K*x) / std::tanh(r);
}
inline float positive_signal(float x) {
    return (x+1)/2;
}
inline float negative_signal(float x) {
    return (x-1)/2;
}
inline float sigmoid(float x, float K=10) {
    return 1.0f / (1.0f + std::exp(-K*x));
}
inline float sigmoid_minus(float x) {
    return -sigmoid(x);
}
inline float bpsigmoid(float x) {
    return (2*sigmoid(x))-1.0f;
}
inline float full_rectify(float x) {
    return amp_clamp(x*std::abs(x),0,1);
}
inline float half_rectify(float x) {
    return amp_clamp(x,0,1);
}
inline float modulated_signals(float a, float b) {
    return a*b;
}
inline float circular_modulated_signals(float a, float b) {
    return std::fmod(a,b);
}
inline float positive_modulated_signals(float a, float b) {
    return positive_signal(a)*positive_signal(b);
}
inline float negative_modulated_signals(float a, float b) {
    return negative_signal(a)*negative_signal(b);
}


class Cubic 
{
public:
  //! Default constructor.
  Cubic( void ) : a1_(0.5), a2_(0.5), a3_(0.5), gain_(1.0), threshold_(1.0) {};

  //! Set the a1 coefficient value.
  void setA1( float a1 ) { a1_ = a1; };

  //! Set the a2 coefficient value.
  void setA2( float a2 )  { a2_ = a2; };

  //! Set the a3 coefficient value.
  void setA3( float a3 )  { a3_ = a3; };

  //! Set the gain value.
  void setGain( float gain ) { gain_ = gain; };

  //! Set the threshold value.
  void setThreshold( float threshold ) { threshold_ = threshold; };

  //! Input one sample to the function and return one output.
  float Tick( float input, float A = 1, float X = -1, float Y = 1 );

  //! Take a channel of the StkFrames object as inputs to the function and replace with corresponding outputs.
  /*!
    The StkFrames argument reference is returned.  The \c channel
    argument must be less than the number of channels in the
    StkFrames argument (the first channel is specified by 0).
    However, range checking is only performed if _STK_DEBUG_ is
    defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  //StkFrames& tick( StkFrames& frames, unsigned int channel = 0 );

  //! Take a channel of the \c iFrames object as inputs to the function and write outputs to the \c oFrames object.
  /*!
    The \c iFrames object reference is returned.  Each channel
    argument must be less than the number of channels in the
    corresponding StkFrames argument (the first channel is specified
    by 0).  However, range checking is only performed if _STK_DEBUG_
    is defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  //StkFrames& tick( StkFrames& iFrames, StkFrames &oFrames, unsigned int iChannel = 0, unsigned int oChannel = 0 );

protected:

  float a1_;
  float a2_;
  float a3_;
  float gain_; 
  float threshold_;
  float lastFrame_;
};

inline float Cubic :: Tick( float input, float A, float X, float Y)
{
  A *= pre_gain;;
  float inSquared = input * input;
  float inCubed = inSquared * input;

  lastFrame_ = gain_ * (a1_ * input + a2_ * inSquared + a3_ * inCubed);

  // Apply threshold if we are out of range.
  if ( fabs( lastFrame_ ) > threshold_ ) {
    lastFrame_ = ( lastFrame_ < 0 ? -threshold_ : threshold_ );
  }
  if(lastFrame_ < X) lastFrame_ = X;
  if(lastFrame_ < Y) lastFrame_ = Y;
  return lastFrame_;
}


struct Clipper
{
    std::function<float (float,float,float,float)> clip;

    Clipper() {
        clip = [](float I, float A, float X, float Y) { return pre_gain*I; };
    }
    Clipper(std::function<float (float,float,float,float)> c) : clip(c) {}

    float Tick(float I, float A =1, float X =1, float Y = 1) {        
        return amp_clamp(post_gain*clip(I,A,X,Y),X,Y);
    }    
};

struct Amplifier
{
    float G;
    float bias;
    float prv;
    float min = -1;
    float max = 1;

    std::function<float (float,float,float,float)> preclip;
    std::function<float (float,float,float,float)> postclip;

    Amplifier(float Gain = 1.0f, float b = 0.0f)
    {
        G = Gain;
        bias = b;
        preclip = [](float I, float A, float X, float Y) { return I; };
        postclip= [](float I, float A, float X, float Y) { return I; };
    }
    
    float Integrator(float in) {
        float r = in + prv;
        prv = in;
        return r;
    }
    float Differencer(float in) {
        float r = in - prv;
        prv = in;
        return r;
    }
    void SetBias(float b) {
        bias = b;
    }
    float Tick(float I, float A = 1, float X = 0, float Y = 0, float B=0)
    {
        float x = G*preclip(pre_gain*(I + (B+bias)),A,X,Y);        
        return amp_clamp((G*postclip(post_gain*x,A,X,Y)),min,max);
    }
};



inline float sigmoidDistortionFunction(float x, float gain,
                                        float max, float dc) {                                                        
    return max * gain * x / sqrt(1 + (gain * std::pow(gain * x, 2))) + dc;
}

inline float asymmetricSigmoidDistortionFunction(float x) 
{
    // Cutoff for chopping top
    static float cutoff = 0.05;
    static float slope = 0.1;
    static float gain = 20;
    static float max = 0.3;
    static float dc = 0;
    // Calculate constant to add to linear region to make it join up with the
    // sigmoid function
    static float b = sigmoidDistortionFunction(x, gain, max, dc) - slope * cutoff;
    if (x > cutoff) {
        return slope * x + b;
    } else {
        return sigmoidDistortionFunction(x, gain, max, dc);
    }
}
inline float assymetric_sigmoid(float I, float A = 1, float X = -1, float Y = 1) {        
    float x = asymmetricSigmoidDistortionFunction(A*pre_gain*I);
    return amp_clamp(post_gain*x,X,Y);
}
inline float asymmetricSigmoidDistortionFunction2(float x) {
    // Cutoff for chopping top
    static float cutoff = 0.05;
    static float gain = 20;
    static float max = 0.3;
    static float dc = 0;
    
    if (x > cutoff) {
        return sigmoidDistortionFunction(sigmoidDistortionFunction(x, gain, max, dc), gain * 2, max, dc);
    } else {
        return sigmoidDistortionFunction(x, gain, max, dc);
    }
}
inline float assymetric_sigmoid2(float I, float A = 1, float X = -1, float Y = 1) {
    float x = asymmetricSigmoidDistortionFunction2(A*pre_gain*I);
    return amp_clamp(post_gain*x,X,Y);
}

inline float distortionFunction(float x) {
    if (x < -0.08905) {
        // Assume x >= -1
        // Therefore, this first interval is actually -1 <= x < -0.08905
        return -(3 / 4) * (1 - std::pow((1 - (-x - 0.032847)), 12) +
                        (1 / 3) * (-x - 0.032847)) +
            0.01;
    } else if (x < 0.320018) {
        return -6.153 * std::pow(x, 2) + 3.9375 * x;
    } else {
        // Assume x <= 1
        // Therefore, this last interval is actually 0.320018 <= x <= 1
        return 0.630035;
    }
}

inline float distortion_function(float I, float A=1,float X = -1, float Y=1) 
{
    float x = distortionFunction(A*pre_gain*I);
    return amp_clamp(post_gain*x,X,Y);
}


float cubic_distortion(float in, float A = 1, float X = -1, float Y = 1)
{
    A *= pre_gain;
    float r = A*(in - (1.0f/3.0f)*std::pow(in,3.0f));    
    return amp_clamp(post_gain*r,X,Y);
}
float asin_distortion(float in, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float r = (2.f / M_PI) * std::asin(in * A);
    if(r < X) r = X;
    if(r > Y) r = Y;
    return amp_clamp(post_gain*r,X,Y);
}
float acos_distortion(float in, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float r = (2.f / M_PI) * std::acos(in * A);
    if(r < X) r = X;
    if(r > Y) r = Y;
    return amp_clamp(post_gain*r,X,Y);
}
float atan_distortion(float in, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float r = (2.f / M_PI) * std::atan(in * A);
    if(r < X) r = X;
    if(r > Y) r = Y;
    return amp_clamp(post_gain*r,X,Y);
}
float asinh_distortion(float in, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float r = (2.f / M_PI) * std::asinh(in * A);    
    return amp_clamp(A*r,X,Y);    
}
float acosh_distortion(float in, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float r = (2.f / M_PI) * std::acosh(in * A);
    if(r < X) r = X;
    if(r > Y) r = Y;
    return amp_clamp(post_gain*r,X,Y);
}
float atanh_distortion(float in, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float r = (2.f / M_PI) * std::atanh(in * A);
    return amp_clamp(post_gain*r,X,Y);
}
float exp_distortion(float x, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
    float sign = x < 0? -1.0f:1.0f;
    float r = sign * (1.f - std::exp(-std::fabs(x * A)));
    return amp_clamp(post_gain*r,X,Y);
}
float dc_distortion(float x, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;;
    float dc = (float)rand() / (float)(RAND_MAX);
    float r = 0;
    if(x < 0) r = cubic_distortion(x - X*dc,A,X,Y);
    else r = cubic_distortion(x+Y*dc,A,X,Y);
    return amp_clamp(post_gain*r,X,Y);
}
float bipolar_distortion(float x, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;;
    float r = 0;
    if(x > 0) r = atan_distortion(x,A,X,Y);
    else r = cubic_distortion(x,A,X,Y);
    return amp_clamp(post_gain*r,X,Y);
}
float quadratic_distortion(float x, float A=1, float X = -1, float Y=1)
{    
    A *= pre_gain;
    float r = x;
    if(x >= 0 && x < M_PI/2) r = atan_distortion(x,A,X,Y);
    else if(x >= (M_PI/2) && x < M_PI) r = atan_distortion(x,A,X,Y);
    else if(x >= M_PI && x < (3*M_PI/4)) r = exp_distortion(x,A,X,Y);
    else r= exp_distortion(x,A,X,Y);
    return amp_clamp(post_gain*r,X,Y);
}
float quadratic2_distortion(float x, float A=1, float X = -1, float Y=1)
{    
    A *= pre_gain;
    float r = x;
    if(x >= 0 && x < M_PI/2) r= atan_distortion(x,A,X,Y);
    else if(x >= (M_PI/2) && x < M_PI) r= cubic_distortion(x,A,X,Y);
    else if(x >= M_PI && x < (3*M_PI/4)) r= atan_distortion(x,A,X,Y);
    else r= cubic_distortion(x,A,X,Y);
    return amp_clamp(post_gain*r,X,Y);
}
float quadratic3_distortion(float x, float A=1, float X = -1, float Y=1)
{    
    A *= pre_gain;
    float r = x;
    if(x >= 0 && x < M_PI/2) r=cubic_distortion(x,A,X,Y);
    else if(x >= (M_PI/2) && x < M_PI) r=cubic_distortion(x,A,X,Y);
    else if(x >= M_PI && x < (3*M_PI/4)) r=atan_distortion(x,A,X,Y);
    else r= atan_distortion(x,A,X,Y);
    return amp_clamp(post_gain*r,X,Y);
}
float parametric_clip(float input, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;
	float softClip = (2.f / M_PI) * std::atan(input * A);
    float blend = input * (X*0.5 + softClip * Y*0.5);
    return amp_clamp(post_gain*blend,X,Y);
}


float arcTanDistortion (float input, float A=1, float X = -1, float Y = 1)
{
    A *= pre_gain;
    float gain = A + 1.0f;    
    float out = 2.0f / M_PI * std::atan(gain * input);    
    out = out / std::log(gain);
    return amp_clamp(post_gain*out,X,Y);
}

float softClipper (float input, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;;
    float newInput = input * A;
    float out = 0.0;
    
    if (newInput >= 1.0f)
        out = 1.0f;
    else if ((newInput > -1) && (newInput < 1))
        out = (3.0f / 2.0f) * (newInput - (std::pow(newInput, 3.0f) / 3.0f));
    else if (newInput <= -1)
        out = -1.0f;
    
    return amp_clamp(post_gain*out,X,Y);    
}

float errorf(float x, float K = 10, float X =-1, float Y=1) {
    float r = std::erf(pre_gain*K*x);
    r *= post_gain;
    return amp_clamp(r,X,Y);
}

float sigmoid (float input, float A=1, float X = -1, float Y=1)
{
    A *= pre_gain;;
    float gain = gain + 1.0f;    
    float out = (2.0f * (1.0f / (1.0f + std::exp(-gain * input)))) - 1;    
    out = (out) / (std::log(gain));
    return amp_clamp(post_gain*out,X,Y);
}

float hardclip(float input, float A=1, float X=-1, float Y=1) {

    input *= A*pre_gain;
    return amp_clamp(post_gain*input,X,Y);
}

float hyperbolicTangent (float input, float gain=1, float X = -1, float Y=1)
{
    gain = (gain + 1.0f);
    float out = (std::tanh(gain * pre_gain * input)) / (std::tanh(gain));
    return amp_clamp(post_gain*out,X,Y);
}

float diodeClipping (float input, float gain=1, float X=-1, float Y=1)
{
//    gain = gain + 1.0f;
    
    float diodeClippingAlgorithm = std::exp((0.1f * pre_gain * input) / (0.0253f * 1.68f)) - 1.0f;
    float out = 2 / M_PI * std::atan(diodeClippingAlgorithm * (gain * 16));
    return amp_clamp(post_gain*out,X,Y);
}

float fuzzExponential (float input, float gain=1, float X =1, float Y=1)
{
    gain *= pre_gain;
    float newInput = input * gain;
    float out;
    
    //Soft clipping
    if (newInput < 0.0f)
        out = -1.0f *  (1.0f - std::exp(-abs(newInput)));
    else
        out = 1.0f * (1.0f - std::exp(-abs(newInput)));
 
    //Half Wave Rectifier
    out = 0.5f * (out + abs(out));
    return amp_clamp(post_gain*out,X,Y);

}

float pieceWiseOverdrive (float input, float gain, float X =-1, float Y=1)
{
    gain = (gain + 1.0f);
    float newInput = pre_gain*input * (gain) ;
    float out = 0.0f;
    
    if (abs(newInput) <= 1.0f / 3.0f)
        out = 2.0f * newInput;
    else if (abs(newInput) > 2.0f / 3.0f)
    {
        if (newInput > 0.0f)
            out = newInput;
        if (newInput < 0.0f)
            out = -newInput;
    } else
    {
        if (newInput > 0.0f)
            out = (3.0f - std::pow((2.0f - newInput * 3.0f), 2.0f)) / 3.0f;
        if (newInput < 0.0f)
            out = -(3.0f - std::pow((2.0f - newInput * 3.0f), 2.0f)) / 3.0f;
    }
    
    out = (out / std::log(gain + 1.0f));
    return amp_clamp(post_gain*out,X,Y);    
}

float tube (float input, float gain, float X =-1, float Y=1)
{
    gain = (gain + 1.0f);
    float Q = -1.5f; //more negative = more linear
    float distortion = 5; //higher number = higher distortion
    float out;
    
    float newInput = pre_gain * input * (gain / 10);
    
    if (Q == 0)
    {
        out = newInput / (1 - std::exp(-distortion * newInput));
        if (newInput == Q)
        {
            out = 1 / distortion;
        }
    } else
    {
        out = ((newInput - Q) / (1 - std::exp(-distortion * (newInput - Q)))) + (Q / (1 - std::exp(distortion * Q)));
        if (newInput == Q)
        {
            out = (1 / distortion) + (Q / (1 - std::exp(distortion * Q)));
        }
    }
    
    return amp_clamp(post_gain*out,X,Y);
}

float arraya (float input, float gain, float X =-1, float Y=1)
{
    gain = (gain + 1.0f);
    auto newInput = pre_gain * input;
    
    //Arraya

    auto out = ((3.0f * newInput) / 2.0f) * (1.0f - (std::pow(newInput, 2.0f) / 3.0f));
    
//    Fuzz Exponential
    if (out < 0.0f)
        out = 1.0f * ((1.0f - std::exp(abs(out)) / (std::exp(1.0f) - 1.0f)));
    else
        out = -1.0f * ((1.0f - std::exp(abs(out)) / (std::exp(1.0f) - 1.0f)));
    
    //Exponential 2
//    out = (std::exp(1.0f) - std::exp(1.0f - out)) / (std::exp(1.0f) - 1.0f);
    
//    out = 0.5f * (out + abs(out));
//    out = abs(out);
    
    if (gain >= 10.0f)
        out = out * (gain / 100.0f);
    else
        out = out * (0.1f);
    
    //Arraya
    out = ((3.0f * out) / 2.0f) * (1.0f - (std::pow(out, 2.0f) / 3.0f));
    return amp_clamp(post_gain*out,X,Y);
}

float gallo (float input, float gain, float X =-1, float Y=1)
{
    gain = (gain + 1.0f);
    float a = -0.01f;
    float b = 0.7f;
    float k1 = std::pow(a, 2.0f);
    float k2 = 1 + (2 * a);
    float k3 = std::pow(b, 2.0f);
    float k4 = 1 - (2 * b);
    float out_1 = 0.0f;
    
    auto newInput = pre_gain * input * gain;
    
    if (newInput < a)
        out_1 = (k1 + newInput) / (k2 - newInput);
    if (newInput >= a && newInput <= b)
        out_1 = newInput;
    if (newInput > b)
        out_1 = (newInput - k3) / (newInput + k4);
    
    return amp_clamp(post_gain*out_1,X,Y);
}

float doubleSoftClipper (float input, float gain, float X =-1, float Y=1)
{
    gain = (gain + 1.0f);
    auto slope = 2.0f;
    auto upperLim = 0.8f;
    auto lowerLim = -1.0f;
    auto upperSkew = 1.0f;
    auto lowerSkew = 1.0f;
    auto xOffFactor = 0.0f;
    auto out = 0.0f;
    
    auto xOff = (1.0f / slope) * std::pow(slope, xOffFactor);
    
    input *= (gain / 10.0f);
    
    if (input > 0.0f)
    {
        input = (input - xOff) * upperSkew;
        
        if (input >= 1.0f / slope)
            out = upperLim;
        else if (input <= -1.0f / slope)
            out = 0.0f;
        else
            out = (3.0f / 2.0f) * upperLim * (slope * input - std::pow(slope * input, 3.0f) / 3.0f) / 2.0f + (upperLim / 2.0f);
    } else
    {
        input = (input + xOff) * lowerSkew;
        
        if (input >= 1.0f / slope)
            out = 0.0f;
        else if (input <= -1.0f / slope)
            out = lowerLim;
        else
            out = (3.0f / 2.0f) * -lowerLim * (slope * input - std::pow(slope * input, 3.0f) / 3.0f) / 2.0f + (lowerLim / 2.0f);
    }
    if(out < X) out = X;
    if(out > Y) out = Y;
    return out;
}

float crush (float input, float gain, float X =-1, float Y=1)
{
    gain = pre_gain*(gain + 1.0f);
    auto out = 0.0f;
    
    gain /= 100.0f;
    
    float dry = input;
    float wet = round(input * std::pow(2, gain));
    out = (wet + dry)  * asin(gain) + dry;
    return amp_clamp(post_gain*out,X,Y);
}

float tuboid (float input, float gain, float X =-1, float Y=1)
{
    gain = pre_gain*(gain + 1.0f);
    auto ktp = 1.0f;
    auto ktn = 3.0f;
    auto sfn = 0.0f;
    
    auto threshPos = 0.3f;
    auto threshNeg = -0.7f;
    
    auto out = 0.0f;
    
    gain /= 10.0f;
    
    auto so = input * gain;
    
    if (so >= threshPos)
        sfn = ktp * std::pow(so - threshPos, 3.0f);
    else if (so <= threshNeg)
        sfn = -ktn * abs(std::pow(so - threshNeg, 3.0f));
    else
        sfn = 0.0f;
    
    so = (input - sfn) * gain;
    out = so;    
    return amp_clamp(post_gain*out,X,Y);
}

float pakarinen_Yeh (float input, float gain, float X =-1, float Y=1)
{
    gain = pre_gain*(gain + 1.0f);
    auto out = 0.0f;
    
    gain /= 100.0f;
    
    auto x = input * gain;
    
    if ((x >= 0.320018f) && (x <= 1.0f))
        out = 0.630035f;
    else if ((x >= -0.08905f) && (x < 0.320018))
        out = (-6.153f * std::pow(x, 2.0f)) + (3.9375f * x);
    else if ((x >= -1.0f) && (x < -0.08905f))
        out = (-0.75f * (1.0f - std::pow(1.0f - (abs(x) - 0.029f), 12.0f) + (0.333f * (abs(x) - 0.029f)))) + 0.01f;
    else
        out = -0.9818f;
    
    out *= 1.5f;
    return amp_clamp(post_gain*out,X,Y);
}


/*
template<typename T>
T logiclip (T x) noexcept
{
    return 2.0f / (1.0f + JMath::exp (-2.0f * x)) - 1.0f;
}

template<typename T>
T hardclip (T x) noexcept
{
    return sgn (x) * std::fminf (std::fabsf(x), 1.0f);
}

template<typename T>
T tanclip (T x) noexcept
{
    float soft = 0.0f;
    return JMath::tanh ((1.0f - 0.5f * soft) * x);
}

template<typename T>
T quintic (T x) noexcept
{
    if (std::fabsf (x) < 1.25f)
    {
        return x - (256.0f / 3125.0f) * std::powf (x, 5.0f);
    } else
    {
        return sgn (x) * 1.0f;
    }
}

template<typename T>
T cubicBasic (T x) noexcept
{
    if (std::fabsf (x) < 1.5f)
    {
        return x - (4.0f / 27.0f) * std::powf (x, 3.0f);
    } else
    {
        return sgn (x) * 1.0f;
    }
}

template<typename T>
T algClip (T x) noexcept
{
    float soft = 0.0f;
    return x / std::sqrtf ((1.0f + 2.0f * soft + std::powf (x, 2.0f)));
}

template<typename T>
T arcClip (T x) noexcept
{
    float soft = 0.0f;
    return (2.0f / juce::MathConstants<T>::pi) * std::atanf ((1.6f - soft * 0.6f) * x);
}

template<typename T>
T sinclip (T x) noexcept
{
    if (std::fabsf (x) < juce::MathConstants<T>::pi)
    {
        return JMath::sin (x);
    }
    else
    {
        return sgn (x) * 1.0f;
    }
}
inline float FuzzCtrTable(const float x)
{
    static auto gen = std::minstd_rand(2112);
    static const float b = 20;

    static auto dist = std::uniform_real_distribution<float>(-1.0, 1.0);

    auto g = exp(-x * x * b);
    auto xadj = x + g * dist(gen);
    return xadj;
}

template <int scale> float FuzzTable(const float x)
{
    static auto gen = std::minstd_rand(2112);
    static const float range = 0.1 * scale;
    static auto dist = std::uniform_real_distribution<float>(-range, range);

    auto xadj = x * (1 - range) + dist(gen);
    return xadj;
}
static float_t absf(float_t x) {
    return (x >= 0.0 ? x : -x);
}

// cube function
static float_t cubef(float_t x) {
    return (x * x * x);
}

// use this to process audio via the rectification algorithm
static float_t Rectify(float_t sample) {
    return ((1 - RectifierThreshold) * sample) + (absf(sample) * RectifierThreshold);
};

// hard clip of input signal
static float_t HardClip(float_t sample, float_t thresh) {
    return 0.5 * (absf(sample + thresh) - absf(sample - thresh));
};

// cubic soft clip function
static float_t SoftCubicClip(float_t sample, float_t thresh) {
    float_t threshInv = 1 / thresh;
    return threshInv * ((thresh * 1.5 * HardClip(sample, thresh)) -
        (0.5 * cubef(HardClip(sample, thresh)) * threshInv));
};

// use this to process audio via the SoftCubicClip algorithm
static float_t SoftCubic(float_t sample) {
    return (invsqrt2 / 3) * (SoftCubicClip(sample, CubicSoftClipThreshold) +
        (CubicHarmonicBalance * SoftCubicClip(absf(sample), CubicSoftClipThreshold)));
};

// soft clip function with adjustable knee
static float_t SKClip(float_t sample, float_t knee) {
    return sample / (knee * absf(sample) + 1.0);
};

// use this to process audio via the SKClip algorithm
static float_t SoftKnee(float_t sample) {
    return 0.5 * (SKClip(sample, SoftClipKnee) + ((SoftClipKnee / 2.0) * SKClip(absf(sample), SoftClipKnee)));
};

// use this to process audio via the leaky integrator algorithm
static float_t LeakyInt(float_t sample, float_t previous_sample) {
    if (sample > previous_sample) {
       return invsqrt2 * (((1.0 - TcRise) * sample) + (TcRise * previous_sample));
    }
      else {
       return invsqrt2 * (((1.0 - TcFall) * sample) + (TcFall * previousfloat linearScale( float in, float min, float max )
{
	float ret;
	if ( min == 0.0f && max == 0.0f )
	{
		ret = 0.0f;
	}
	else if ( min > max )
	{
		ret = min - ( in * ( min - max ) );
	}
	else
	{
		ret = min + ( in * ( max - min ) );
	}
	return ret;
}

float linearDescale( float in, float min, float max )
{
	float ret;
	if ( min == 0.0f && max == 0.0f )
	{
		ret = 0.0f;
	}
	else if ( min > max )
	{
		ret = ( min - in ) / ( min - max );
	}
	else
	{
		ret = ( in - min ) / ( max - min );
	}
	return ret;
}

float expoScale( float in, float min, float max )
{
	// negative log makes no sense...
	if ( min < 0.0f || max < 0.0f )
	{
		return 0.0f;
	}

	// not handling min > max (inverse) case yet

	// figure out how many "octaves" (doublings) it takes to get from min to
	// max
	// we only have log & log10 so we have to do change of base
	// note this uses + instead of * so we can handle min == 0
	float octaves = log( max - min + 1 ) / log( 2.0f );
	return ( min - 1 ) + pow( 2.0f, in * octaves );
}

float expoDescale( float in, float min, float max )
{
	// see above
	if ( min < 0.0f || max < 0.0f )
	{
		return 0.0f;
	}

	// again, not handling min > max (inverse) case yet
	
	// note this was derived by simply inverting the previous function
	float log2 = log( 2.0f );
	return ( log( in - min + 1 ) / log2 ) / ( log( max - min + 1 ) / log2 );
}

float floorScale( float in, float min, float max )
{
	if ( min > max )
	{
		return ceil( linearScale( in, min, max ) );
	}
	else
	{
		return floor( linearScale( in, min, max ) );
	}
}

float expoShape( float in, float amount )
{
	if ( in == 0.0f )
		return in;

	float flip = in < 0.0f ? -1.0f : 1.0f;

	return pow( in * flip, amount ) * flip;
}

float softClipShape( float in, float amount )
{
	return in / ( 1 + fabs( in ) );
}

float sineShape( float in, float amount )
{
	return in * cos( in * amount );
}

float chebyshevShape( float in, float amount )
{
	return chebyshevRec( in, (int)amount );
}

float chebyshevRec( float in, int depth )
{
	if ( depth == 0 )
	{
		return 1.0f;
	}

	// lastval represents C(k-1)float Distortion::softClip(float sample)
{
    if (sample < -1.f) {
        return -softClipThreshold; //2/3
    }
    else if (sample > 1.f) {
        return softClipThreshold;
    }
    else {
        return sample - ((sample * sample * sample) / 3.f);
    }
}



// Arctangent nonlinearity
float Distortion::arctangent(float sample, float alpha)
{
    // f(x) = (2 / PI) * arctan(alpha * x[n]), where alpha >> 1 (drive param)
    return (2.f / PI)* atan(alpha * sample);
}

// Hard-clipping nonlinearity
float Distortion::hardClip(float sample)
{
    if (sample < -1.f) {
        return -1.f;
    }
    else if (sample > 1.f) {
        return 1.f;
    }
    else {
        return sample;
    }
}

// Square law series expansion
float Distortion::squareLaw(float sample, float alpha)
{
    return sample + alpha * sample * sample;
}

float Distortion::cubicWaveShaper(float sample)
{
    return 1.5f * sample - 0.5f * sample * sample * sample;
}

// Foldback nonlinearity, input range: (-inf, inf)
float Distortion::foldback(float sample)
{
    // Threshold should be > 0.f
    if (sample > controls.threshold || sample < -controls.threshold) {
        sample = fabs(fabs(fmod(sample - controls.threshold,
                                controls.threshold * 4))
                      - controls.threshold * 2) - controls.threshold;
    }
    return sample;
}

// A nonlinearity by Partice Tarrabia and Bram de Jong
float Distortion::waveShaper1(float sample, float alpha)
{
    const float k = 2.f * alpha / (1.f - alpha);
    return (1.f + k) * sample / (1.f + k * fabs(sample));
}

// A nonlinearity by Jon Watte
float Distortion::waveShaper2(float sample, float alpha)
{
    const float z = PI * alpha;
    const float s = 1.f / sin(z);
    const float b = 1.f / alpha;
    
    if (sample > b) {
        return 1.f;
    }
    else {
        return sin(z * sample) * s;
    }
}

// A nonlinearity by Bram de Jong, input range: [-1, 1]
float Distortion::waveShaper3(float sample, float alpha)
{
    // original design requires sample be positive
    // alpha: [0, 1]
    bool isNegative = false;
    float output = sample;
    if (sample < 0.f) {
        isNegative = true;
        output = -output;
    }
    
    if (output > alpha) {
        output = alpha + (output - alpha)
            / (1.f + powf(((output - alpha)) / (1.f - alpha), 2.f));
    }
    if (output > 1.f) {
        output = (alpha + 1.f) / 2.f;
    }
    
    if (isNegative) {
        output = -output;
    }
    
    return output;
}


float Distortion::gloubiBoulga(float sample)
{
    const double x = sample * 0.686306;
    const double a = 1 + exp(sqrt(fabs(x)) * -0.75);
    return (exp(x) - exp(-x * a)) / (exp(x) + exp(-x));
}

// Approximation based on description in gloubiBoulga
float Distortion::gloubiApprox(float sample)
{
    return sample - (0.15f * sample * sample) - (0.15f * sample * sample * sample);
}
	float lastVal = 1.0f;
	float out = in;
	float temp;

	// depth=1 is the base case
	for( int i = 1; i < depth; i++ )
	{
		temp = out;
		out = ( 2.0f * in * out ) - lastVal;
		lastVal = temp;
	}
	return out;
}_sample));
    }
};



inline float FuzzEdgeTable(const float x)
{
    static auto gen = std::minstd_rand(2112);
    static auto dist = std::uniform_real_distribution<float>(-1.0, 1.0);

    auto g = x * x * x * x;
    auto xadj = 0.85 * x + 0.15 * g * dist(gen);
    return xadj;
}

template<typename T>
T limitclip (T x) noexcept
{
    return juce::jlimit (-0.1f, 0.1f, x);
}

switch ((int)paramDistortionType.getTargetValue()) {                    
        //hard clipping
        case distortionTypeHardClipping: {
            float threshold = 0.5f;
            if (in > threshold)
                out = threshold;
            else if (in < -threshold)
                out = -threshold;
            else
                out = in;
            out *= 0.5f;
            break;
        }
            
        //soft clipping
        case distortionTypeSoftClipping: {
            float threshold1 = 1.0f / 3.0f;
            float threshold2 = 2.0f / 3.0f;
            if (in > threshold2)
                out = 1.0f;
            else if (in > threshold1)
                out = 1.0f - powf (2.0f - 3.0f * in, 2.0f) / 3.0f;
            else if (in < -threshold2)
                out = -1.0f;
            else if (in < -threshold1)
                out = -1.0f + powf (2.0f + 3.0f * in, 2.0f) / 3.0f;
            else
                out = 2.0f * in;
            out *= 0.5f;
            break;
        }
            
        //exponential
        case distortionTypeExponential: {
            if (in > 0.0f)
                out = 1.0f - expf (-in);
            else
                out = -1.0f + expf (in);
            out *= 0.05f;
            break;
        }
        
        //full-wave rectifier
        case distortionTypeFullWaveRectifier: {
            out = fabsf (in);
            break;
        }
            
        //half-wave rectifier
        case distortionTypeHalfWaveRectifier: {
            if (in > 0.0f)
                out = in;
            else
                out = 0.0f;
            break;
        }
         
        //fold-back
        case distortionTypeFoldBack: {
            float fblevel1 = 0.3f;
            float fblevel2 = fblevel1 * 2.0f;
            if (in > fblevel1)
                out = (fblevel2 - out) * 0.05f;
            else if (in < (-fblevel1)) {
                out = (-fblevel2 - out) * 0.05f;
            }
            break;
        }

        //squarer
        case distortionTypeSquarer: {
            out = powf(in, 2.0f);
            break;
        }
            
        //chebyshev 4th order
        case distortionTypeChebT4: {
            out = (8.0f * powf(in, 4.0f)) - (8.0f * powf(in, 2.0f)) - 1.0f;
            out = out * 0.1f;
            break;
        }
            
        //bit crusher
        case distortionTypeBitCrusher: {
            if (cpt == 0.0f)
                out = in;
            elseconst auto softClip = piDivisor * std::atanf(input * juce::Decibels::decibelsToGain(driveDB));
    
    auto blend = input * (1.0 - mix) + softClip * mix;
    
    blend *= juce::Decibels::decibelsToGain(outputDB);
    
                out = 0.0f;
            cpt = (cpt + 1) % undFac;
            break;
        }
            
                //slew limiter
                case distortionTypeSlewLimiter: {
                    if (in > out)
                        out = jmin (in, out + slewRise);
                    else
                        out = jmax (in, out - slewFall);
                    break;
                }
            } 
            
            
    const float eta = 2.f;
    const float Is = 1.e-6;
    const float Vt = 26.e-3;

    float Fs;
    float Ts;

    float C2;
    float R2;
    float x2;

    float R5;
    
    float potLev;

    float Vd;
    float thr;

    // Grouped Resistances
    float G;
                                 
    float Clipping::processSample(float Vi) {
    float b = 1.f; // for dampening
    float fd = -Vi / R2 + Is * sinh(Vd / (eta*Vt)) + G * Vd - x2;
    for (int i = 0; i < 50 && abs(fd) > thr; ++i) {
        float fdd = (Is / (eta*Vt)) * cosh(Vd / (eta*Vt)) + G;
        float Vnew = Vd - b * fd / fdd;
        float fn = -Vi / R2 + Is * sinh(Vnew / (eta*Vt)) + G * Vnew - x2;
        if (abs(fn) < abs(fd)) {
            Vd = Vnew;
            b = 1.f;
        }
        else {
            b *= 0.5f;
        }

        fd = -Vi / R2 + Is * sinh(Vd / (eta*Vt)) + G * Vd - x2;
    }
    
    x2 = 2 * Vd / R2 - x2;
    return  potLev * Vd;
}
void Distortion::updateCoefficients() {
    Ts = 1.f / Fs;
    R1 = Ts / (2.f*C1);
    updateGroupedResistances();
}

void Distortion::updateGroupedResistances() {
    G = 1.f / (R1 + R3 + Rp);
    Gb = (R3 + Rp) * G;
    Gi = 1.f + R4 * G;
    Gx1 = R1 * R4 * G;
}
float Distortion::processSample(float Vi) {
    float Vb = Gb * Vi - R1 * Gb*x1;
    float Vr1 = Vi - Vb;
    float Vo = Gi * Vi - Gx1 * x1;
    if (Vo > 4.5f) {
        Vo = 4.5f;
    } else if (Vo < -4.5f) {
        Vo = -4.5f;
    }
    x1 = (2.f * Vr1 / R1) - x1;
    return Vo;
}

struct table1d { // 1-dimensional function table
    float low;
    float high;
    float istep;
    int size;
    float data[];
};

template <int tab_size>
struct table1d_imp {
    float low;
    float high;
    float istep;
    int size;
    float data[tab_size];
    operator table1d&() const { return *(table1d*)this; }
};
static table1d_imp<200> tube_table __rt_data = {
	0,0.880571,55,200, {
	0.000000000000,0.016750418760,0.033500837521,0.050251256281,0.067001675042,
	0.083752093802,0.100502512563,0.117252931323,0.134003350084,0.150753768844,
	0.167504187605,0.184254606365,0.201005025126,0.217755443886,0.234505862647,
	0.251256281407,0.268006700168,0.284757118928,0.301507537688,0.318257956449,
	0.335008375209,0.351758793970,0.368509212730,0.385259631491,0.402010050251,
	0.418760469012,0.435510887772,0.452261306533,0.469011725293,0.485762144054,
	0.502512562814,0.519262981575,0.536013400335,0.552763819095,0.569514237856,
	0.586264656616,0.593999365915,0.599005700058,0.603333878092,0.607267980657,
	0.610930403272,0.614388585430,0.617684834322,0.620847956424,0.623898713479,
	0.626852710073,0.629722059392,0.632516408598,0.635243602944,0.637910133907,
	0.640521451968,0.643082191206,0.645596334502,0.648067337630,0.650498224157,
	0.652891659201,0.655250007565,0.657575380145,0.659869671394,0.662134589886,
	0.664371683478,0.666582360199,0.668767905735,0.670929498165,0.673068220463,
	0.675185071166,0.677280973536,0.679356783465,0.681413296334,0.683451252980,
	0.685471344926,0.687474218971,0.689460481240,0.691430700768,0.693385412682,
	0.695325121033,0.697250301337,0.699161402840,0.701058850561,0.702943047129,
	0.704814374447,0.706673195193,0.708519854182,0.710354679615,0.712177984199,
	0.713990066189,0.715791210325,0.717581688703,0.719361761564,0.721131678031,
	0.722891676777,0.724641986646,0.726382827231,0.728114409400,0.729836935788,
	0.731550601256,0.733255593312,0.734952092506,0.736640272797,0.738320301896,
	0.739992341580,0.741656547999,0.743313071948,0.744962059131,0.746603650406,
	0.748237982013,0.749865185791,0.751485389382,0.753098716418,0.754705286704,
	0.756305216386,0.757898618110,0.759485601172,0.761066271663,0.762640732600,
	0.764209084056,0.765771423280,0.767327844807,0.768878440572,0.770423300008,
	0.771962510148,0.773496155712,0.775024319196,0.776547080961,0.778064519305,
	0.779576710544,0.781083729079,0.782585647471,0.784082536500,0.785574465234,
	0.787061501083,0.788543709858,0.790021155826,0.791493901762,0.792962008997,
	0.794425537468,0.795884545762,0.797339091160,0.798789229679,0.800235016113,
	0.801676504067,0.803113746001,0.804546793259,0.805975696107,0.807400503763,
	0.808821264430,0.810238025325,0.811650832709,0.813059731915,0.814464767373,
	0.815865982638,0.817263420412,0.818657122571,0.820047130187,0.821433483549,
	0.822816222186,0.824195384887,0.825571009720,0.826943134054,0.828311794573,
	0.829677027299,0.831038867606,0.832397350237,0.833752509322,0.835104378390,
	0.836452990391,0.837798377701,0.839140572145,0.840479605004,0.841815507032,
	0.843148308470,0.844478039054,0.845804728029,0.847128404163,0.848449095754,
	0.849766830645,0.851081636231,0.852393539473,0.853702566904,0.855008744674,
	0.856312098429,0.857612653515,0.858910434856,0.860205466994,0.861497774101,
	0.862787379985,0.864074308096,0.865358581537,0.866640223071,0.867919255124,
	0.869195699798,0.870469578875,0.871740913821,0.873009725798,0.874276035666,
	0.875539863991,0.876801231052,0.878060156841,0.879316661078,0.880570763209
	}
};

double always_inline tubeclip(double x) {
    double f = fabs(x);
    f = f * tube_table.istep;
    int i = static_cast<int>(f);
    if (i < 0) {
        f = tube_table.data[0];
    } else if (i >= tube_table.size-1) {
        f = tube_table.data[tube_table.size-1];
    } else {
	f -= i;
	f = tube_table.data[i]*(1-f) + tube_table.data[i+1]*f;
    }
    return copysign(f, x);
}

float DistortionAudioProcessor::hardClipping(const float& _in)
{
    float out;
    float threshold = 0.5f;
    
    if (_in > threshold)
        out = threshold;
    else if (_in < -threshold)
        out = -threshold;
    else
        out = _in;
    
    return out;
}

float DistortionAudioProcessor::softClipping(const float& _in)
{
    float out;
    float threshold1 = 1.0f / 3.0f;
    float threshold2 = 2.0f / 3.0f;
    
    if (_in > threshold2)
        out = 1.0f;
    else if (_in > threshold1)
        out = 1.0f - powf (2.0f - 3.0f * _in, 2.0f) / 3.0f;
    else if (_in < -threshold2)
        out = -1.0f;
    else if (_in < -threshold1)
        out = -1.0f + powf (2.0f + 3.0f * _in, 2.0f) / 3.0f;
    else
        out = 2.0f * _in;
        out *= 0.5f;
    
    return out;
}

float DistortionAudioProcessor::exponential(const float& _in)
{
    float out;
    
    if (_in > 0.0f)
        out = 1.0f - expf (-_in);
    else
        out = -1.0f + expf (_in);
    
    return out;
}

float DistortionAudioProcessor::fullWaveRectifier(const float& _in)
{
    float out;
    
    out = fabsf (_in);
    
    return out;

}

float DistortionAudioProcessor::halfWaveRectifier(const float& _in)
{
    float out;
    
    if (_in > 0.0f)
        out = _in;
    else
        out = 0.0f;
    
    return out;
    
}

float DistortionAudioProcessor::ArayaAndSuyama(const float &_in)
{
    float out;
    
    out = (3/2) * (_in) * (1 - pow(_in, 2)/3);
    out = (3/2) * (out) * (1 - pow(out, 2)/3);
    out = (3/2) * (out) * (1 - pow(out, 2)/3);

    return out;
    
}

float DistortionAudioProcessor::doidicSymmetric(const float& _in)
{
    float out;
    
    out =  ( (2*fabsf(_in))  - pow(_in, 2)) * copysignf(1, _in);
    
    return out;
}

float DistortionAudioProcessor::doidicAssymetric(const float& _in)
{
    float out;
    
    if (_in >= -1 && _in < -0.08905f) {
        out = -(0.75)*( 1 - (1 - pow(fabs(_in) - 0.032847, 12)) + 1/3*(fabs(_in) - 0.032847)) + 0.01;
    }
    else if (_in >= -0.08905f && _in < 0.320018)
    {
        out = -6.153 * pow(_in,2) + 3.9375 * _in;
    }
    else if (_in >= 0.320018 && _in <= 1)
    {
        out = 0.630035;
    }
    
    return out;
    
}
*/