#pragma once

#include <vector>

struct ringbuffer
{
    std::vector<float> buffer;
    int read_cursor,write_cursor;

    ringbuffer() = default;

    ringbuffer(size_t n) {
        buffer.resize(n);
        zeros(buffer);
        read_cursor = 0;
        write_cursor = 0;
    }
    size_t size() const { 
        return buffer.size();
    }
    void resize(size_t s) {
        buffer.resize(s);
        zeros(buffer);
        read_cursor  %= s;
        write_cursor %= s;
    }
    void calc_delay(float time)
    {
        write_cursor = 0;
        read_cursor  = (time * sampleRate);
        read_cursor  %= buffer.size();
    }
    void set_delay(size_t samples)
    {
        write_cursor = 0;
        read_cursor  = samples % buffer.size();
    }
    float read(int cursor)
    {
        return buffer[cursor];
    }
    void write(int cursor, float value)
    {
        buffer[cursor] = value;
    }
    float& operator[](int idx) {
        return buffer[idx];
    }
    float read() 
    {
        float x1 = buffer[read_cursor++];
        read_cursor %= buffer.size();
        float x2 = buffer[read_cursor];
        float frac = x1 - floor(x1);
        float out = x1 + frac * (x2-x1);        
        return out;
    }
    void write(float v)
    {
        buffer[write_cursor++] = v;
        write_cursor %= buffer.size();
    }
};


struct delayline
{
    ringbuffer delay;
    size_t read_cursor,write_cursor;
    float feedback;
    float mix;
    size_t delayLen;
    float delayTime;

    enum InterpType
    {
        None,
        NearestNeighbor,
        Linear,
        Cubic,
        Spline3,
        Spline5,
        Hermite1,
        Hermite2,
        Hermite3,
        Hermite4,
    }
    interpType;

    
    delayline() = default;
    delayline(float delay_time, float delay_size = 1) {
        delay.resize(sampleRate*delay_size);        
        feedback=0.5;
        mix = 0.5;
        delayTime = delay_time;        
        delayLen  = delay_time * sampleRate;        
        read_cursor  = 0;
        write_cursor = delayLen;
        interpType = Linear;
    }

    float& operator[](size_t i) {
        return delay[i % delay.size()];
    }
    
    void setDelaySize(float size) {
        delay.resize(size);        
    }
    void reset() {
        read_cursor  = 0;
        write_cursor = delayLen;   
    }
    void setDelayTime(float f) 
    {
        delayTime = f;
        delayLen  = f;
        write_cursor = (read_cursor + delayLen) % delayLen;
        read_cursor  = 0;
        write_cursor = delayLen;
    }
    void setFeedback(float f) {
        feedback = f;
    }
    void setMix(float m) {
        mix = std::fmod(m,1);
    }
    void resize(size_t n) {
        delay.resize(n);
    }
    virtual float Tick(float I, float A = 1, float X = 1, float Y = 1) {        
        size_t n = read_cursor;
        float output = delay.read(read_cursor);            
        float d1 = A*(I - Y*output*feedback);                
        float f= d1 - floor(d1); //(int)d1;        
        output = Interpolate(n,f);
        delay.write(write_cursor, I*Y + feedback*output);        
        write_cursor %= delayLen;
        read_cursor  %= delayLen;
        return mix*I + (1.0-mix)*output;
    }
    size_t size() { return delay.size(); }    

    float Interpolate(size_t n, float frac)
    {
        switch(interpType)
        {
            case None: return delay[n];
            case NearestNeighbor: return NearestNeighborInterpolate(n,frac);
            case Linear: return LinearInterpolate(n,frac);
            case Cubic: return CubicInterpolate(n,frac);
            case Hermite1: return Hermite1Interpolate(n,frac);
            case Hermite2: return Hermite2Interpolate(n,frac);
            case Hermite3: return Hermite3Interpolate(n,frac);
            case Hermite4: return Hermite4Interpolate(n,frac);
            case Spline3:  return Spline3Interpolate(n,frac);
            case Spline5:  return Spline5Interpolate(n,frac);
        }
        return delay[read_cursor];
    }
    float NearestNeighborInterpolate(size_t n,float frac)
    {
        int   x  = round(frac);
        float x1 = delay[ (n + x) % delayLen];
        return x1;
    }
    float LinearInterpolate(size_t n,float frac)
    {            
        float x1 = delay[n];
        float x2 = delay[ (n+1) % delayLen];
        //float frac = x1 - (int)x1;
        return x1 + ((x2-x1)*frac);
    }    
    // just cubic stuff from musicdsp
    float CubicInterpolate(size_t inpos,float finpos)
    {            
        float xm1 = delay[(inpos - 1) % delayLen];
        float x0 = delay[inpos + 0];
        float x1  = delay[(inpos + 1) % delayLen];
        float x2  = delay[(inpos + 2) % delayLen];
        float a = (3 * (x0-x1) - xm1 + x2) / 2;
        float b = 2*x1 + xm1 - (5*x0 + x2) / 2;
        float c = (x1 - xm1) / 2;
        return (((a * finpos) + b) * finpos + c) * finpos + x0;
    }
    // just another kind of cubials it might be the same kakaloke really
    inline float Hermite1Interpolate(size_t inpos, float x)
    {            
        float y0 = delay[(inpos - 1) % delayLen];
        float y1 = delay[inpos + 0];
        float y2  = delay[(inpos + 1) % delayLen];
        float y3  = delay[(inpos + 2) % delayLen];
        // 4-point, 3rd-order Hermite (x-form)
        float c0 = y1;
        float c1 = 0.5f * (y2 - y0);
        float c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
        float c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

        return ((c3 * x + c2) * x + c1) * x + c0;
    }    
    inline float Hermite2Interpolate(size_t inpos, float x)
    {            
        float y0 = delay[(inpos - 1) % delayLen];
        float y1 = delay[inpos + 0];
        float y2  = delay[(inpos + 1) % delayLen];
        float y3  = delay[(inpos + 2) % delayLen];
        // 4-point, 3rd-order Hermite (x-form)
        float c0 = y1;
        float c1 = 0.5f * (y2 - y0);
        float c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
        float c2 = y0 - y1 + c1 - c3;

        return ((c3 * x + c2) * x + c1) * x + c0;
    }    
    inline float Hermite3Interpolate(size_t inpos, float x)
    {                
            float y0 = delay[(inpos - 1) % delayLen];
            float y1 = delay[inpos + 0];
            float y2  = delay[(inpos + 1) % delayLen];
            float y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            float c0 = y1;
            float c1 = 0.5f * (y2 - y0);
            float y0my1 = y0 - y1;
            float c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
            float c2 = y0my1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
    }    
    inline float Hermite4Interpolate(size_t inpos, float frac_pos)
    {            
        float xm1 = delay[(inpos - 1) % delayLen];
        float x0 = delay[inpos + 0];
        float x1  = delay[(inpos + 1) % delayLen];
        float x2  = delay[(inpos + 2) % delayLen];
        const float    c     = (x1 - xm1) * 0.5f;
        const float    v     = x0 - x1;
        const float    w     = c + v;
        const float    a     = w + v + (x2 - x0) * 0.5f;
        const float    b_neg = w + a;

        return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
    }
    float Spline3Interpolate(size_t inpos, float x)
    {            
        float L1 = delay[(inpos-1)%delayLen];
        float L0 = delay[inpos];
        float H0 = delay[(inpos + 1) % delayLen];
        float H1 = delay[(inpos + 2) % delayLen];
        return L0 + .5f* x*(H0-L1 + x*(H0 + L0*(-2) + L1 + x*( (H0 - L0)*9 + (L1 - H1)*3 + x*((L0 - H0)*15 + (H1 -  L1)*5 +  x*((H0 - L0)*6 + (L1 - H1)*2 )))));
    }
    float Spline5Interpolate(size_t inpos, float x)
    {
        /* 5-point spline*/
        int nearest_sample = delayLen;

        float p0=delay[(nearest_sample-2) % delayLen];
        float p1=delay[(nearest_sample-1) % delayLen];
        float p2=delay[nearest_sample];
        float p3=delay[(nearest_sample+1) % delayLen];
        float p4=delay[(nearest_sample+2) % delayLen];
        float p5=delay[(nearest_sample+3) % delayLen];

        return p2 + 0.04166666666*x*((p3-p1)*16.0+(p0-p4)*2.0
        + x *((p3+p1)*16.0-p0-p2*30.0- p4
        + x *(p3*66.0-p2*70.0-p4*33.0+p1*39.0+ p5*7.0- p0*9.0
        + x *( p2*126.0-p3*124.0+p4*61.0-p1*64.0- p5*12.0+p0*13.0
        + x *((p3-p2)*50.0+(p1-p4)*25.0+(p5-p0)*5.0)))));
    }
};



struct multitapdelayline 
{
    std::vector<size_t> tap_reads;
    size_t taps;
    size_t write_cursor, read_cursor;
    size_t delayLen;
    ringbuffer delay;        
    float feedback;
    float mix;

    enum InterpType
    {
        None,
        NearestNeighbor,
        Linear,
        Cubic,
        Spline3,
        Spline5,
        Hermite1,
        Hermite2,
        Hermite3,
        Hermite4,
    }
    interpType;

    multitapdelayline() = default;
    multitapdelayline (float delayTime, float delayScale = 1.0)
    {            
        delay.resize(delayScale*sampleRate);                
        feedback=0.5;
        mix = 0.5;             
        delayLen   = sampleRate*delayTime;            
        interpType = Linear;
    }
    void setDelaySize(size_t size) {
        delay.resize(size);
    }
    void setDelayTime(float time) {
        delayLen = sampleRate*time;
    }
    void addTap(float t) {
        size_t d = t*sampleRate;            
        tap_reads.push_back(d);                    
        taps++;
    }        
    float Tick(float I, float A=1, float X=1, float Y=1) 
    {       
        float output = 0;
        float sum    = 0;            
        size_t read;
        size_t len;
        for(size_t i = 0; i < taps; i++)
        {                            
            read = tap_reads[i];                                
            output = delay.read(read);                        
            float f  = output - floor(output); //(int)output;                                        
            float x1 = output;
            read = read  % delayLen;
            float x2 = delay[read];
            output = x1 + (f*(x2-x1));
            tap_reads[i] = read;
            sum += output;                
        }         
        if( taps > 0) sum /= taps;        
            
        float x1 = 0;
        size_t r = read_cursor;
        x1 = delay.read(read_cursor);  
        output = x1;
        float f  = output - (int)output;                                                                
        output = 0.5*(output + Interpolate(r,f));
        read_cursor = read_cursor % delayLen;
        delay.write(write_cursor, I + Y*feedback*output);        
        write_cursor %= delayLen;
        return mix*I + (1.0-mix)*(0.5*(sum+output));            
    }
    float Interpolate(size_t n, float frac)
    {
        switch(interpType)
        {
            case None: return delay[n];
            case NearestNeighbor: return NearestNeighborInterpolate(n,frac);
            case Linear: return LinearInterpolate(n,frac);
            case Cubic: return CubicInterpolate(n,frac);
            case Hermite1: return Hermite1Interpolate(n,frac);
            case Hermite2: return Hermite2Interpolate(n,frac);
            case Hermite3: return Hermite3Interpolate(n,frac);
            case Hermite4: return Hermite4Interpolate(n,frac);
            case Spline3:  return Spline3Interpolate(n,frac);
            case Spline5:  return Spline5Interpolate(n,frac);
        }
        return delay[read_cursor];
    }
    float NearestNeighborInterpolate(size_t n,float frac)
    {
        int   x  = round(frac);
        float x1 = delay[ (n + x) % delayLen];
        return x1;
    }
    float LinearInterpolate(size_t n,float frac)
    {            
        float x1 = delay[n];
        float x2 = delay[ (n+1) % delayLen];
        //float frac = x1 - (int)x1;
        return x1 + ((x2-x1)*frac);
    }    
    // just cubic stuff from musicdsp
    float CubicInterpolate(size_t inpos,float finpos)
    {            
        float xm1 = delay[(inpos - 1) % delayLen];
        float x0 = delay[inpos + 0];
        float x1  = delay[(inpos + 1) % delayLen];
        float x2  = delay[(inpos + 2) % delayLen];
        float a = (3 * (x0-x1) - xm1 + x2) / 2;
        float b = 2*x1 + xm1 - (5*x0 + x2) / 2;
        float c = (x1 - xm1) / 2;
        return (((a * finpos) + b) * finpos + c) * finpos + x0;
    }
    // just another kind of cubials it might be the same kakaloke really
    inline float Hermite1Interpolate(size_t inpos, float x)
    {            
        float y0 = delay[(inpos - 1) % delayLen];
        float y1 = delay[inpos + 0];
        float y2  = delay[(inpos + 1) % delayLen];
        float y3  = delay[(inpos + 2) % delayLen];
        // 4-point, 3rd-order Hermite (x-form)
        float c0 = y1;
        float c1 = 0.5f * (y2 - y0);
        float c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
        float c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

        return ((c3 * x + c2) * x + c1) * x + c0;
    }    
    inline float Hermite2Interpolate(size_t inpos, float x)
    {            
        float y0 = delay[(inpos - 1) % delayLen];
        float y1 = delay[inpos + 0];
        float y2  = delay[(inpos + 1) % delayLen];
        float y3  = delay[(inpos + 2) % delayLen];
        // 4-point, 3rd-order Hermite (x-form)
        float c0 = y1;
        float c1 = 0.5f * (y2 - y0);
        float c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
        float c2 = y0 - y1 + c1 - c3;

        return ((c3 * x + c2) * x + c1) * x + c0;
    }    
    inline float Hermite3Interpolate(size_t inpos, float x)
    {                
            float y0 = delay[(inpos - 1) % delayLen];
            float y1 = delay[inpos + 0];
            float y2  = delay[(inpos + 1) % delayLen];
            float y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            float c0 = y1;
            float c1 = 0.5f * (y2 - y0);
            float y0my1 = y0 - y1;
            float c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
            float c2 = y0my1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
    }    
    inline float Hermite4Interpolate(size_t inpos, float frac_pos)
    {            
        float xm1 = delay[(inpos - 1) % delayLen];
        float x0 = delay[inpos + 0];
        float x1  = delay[(inpos + 1) % delayLen];
        float x2  = delay[(inpos + 2) % delayLen];
        const float    c     = (x1 - xm1) * 0.5f;
        const float    v     = x0 - x1;
        const float    w     = c + v;
        const float    a     = w + v + (x2 - x0) * 0.5f;
        const float    b_neg = w + a;

        return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
    }
    float Spline3Interpolate(size_t inpos, float x)
    {            
        float L1 = delay[(inpos-1)%delayLen];
        float L0 = delay[inpos];
        float H0 = delay[(inpos + 1) % delayLen];
        float H1 = delay[(inpos + 2) % delayLen];
        return L0 + .5f* x*(H0-L1 + x*(H0 + L0*(-2) + L1 + x*( (H0 - L0)*9 + (L1 - H1)*3 + x*((L0 - H0)*15 + (H1 -  L1)*5 +  x*((H0 - L0)*6 + (L1 - H1)*2 )))));
    }
    float Spline5Interpolate(size_t inpos, float x)
    {
        /* 5-point spline*/
        int nearest_sample = delayLen;

        float p0=delay[(nearest_sample-2) % delayLen];
        float p1=delay[(nearest_sample-1) % delayLen];
        float p2=delay[nearest_sample];
        float p3=delay[(nearest_sample+1) % delayLen];
        float p4=delay[(nearest_sample+2) % delayLen];
        float p5=delay[(nearest_sample+3) % delayLen];

        return p2 + 0.04166666666*x*((p3-p1)*16.0+(p0-p4)*2.0
        + x *((p3+p1)*16.0-p0-p2*30.0- p4
        + x *(p3*66.0-p2*70.0-p4*33.0+p1*39.0+ p5*7.0- p0*9.0
        + x *( p2*126.0-p3*124.0+p4*61.0-p1*64.0- p5*12.0+p0*13.0
        + x *((p3-p2)*50.0+(p1-p4)*25.0+(p5-p0)*5.0)))));
    }
};  



struct combfilter
{
    delayline delay[2];
    float x1,y,y1;
    float gain[2];
    float delayTime[2];
    
    enum {
        X_index,
        Y_index,
    };

    combfilter(float _g1, float _g2, float _d1, float _d2)
    {
        gain[X_index] = _g1;
        gain[Y_index] = _g2;
        delayTime[0] = _d1;
        delayTime[1] = _d2;
        for(size_t i = 0; i < 2; i++)
        {
            delay[i].setDelaySize(44100);
            delay[i].setDelayTime(delayTime[i]);
        }       
        x1=y=y1=0;
    }    

    // X modulation * depth
    // Y modulation * depth
    float Tick(float I, float A=1, float X = 0, float Y=0)
    {
        float x = I;
        y = x + gain[X_index] * x1 - gain[Y_index] * y1;
        x1 = delay[X_index].Tick(x);
        y1 = delay[Y_index].Tick(y);
        return y;
    }
};


struct iircombfilter
{
    delayline delay;
    float g,x,y,y1;

    iircombfilter(float _g, float _d) 
    {
        delay.setDelaySize(44100);
        delay.setDelayTime(_d);
        g = _g;
        x = y = y1 = 0;
    }
    float Tick(float I, float A = 1, float X = 0, float Y= 0)
    {
        x = I;
        y = x - g*y1;
        y1= delay.Tick(y);
        return y;
    }
};


struct fircombfilter
{
    delayline delay;
    float g,x,x1,y;

    fircombfilter(float _g, float _d) 
    {
        delay.setDelaySize(44100);
        delay.setDelayTime(_d);
        g = _g;
        x = y = x1 = 0;
    }
    float Tick(float I, float A = 1, float X = 0, float Y= 0)
    {
        x = I;
        y = x + g*x1;
        x1= delay.Tick(x);
        return y;
    }
};



struct multitapcombfilter
{
    multitapdelayline delay[2];
    float x1,y,y1;
    float gain[2];
    float delayTime[2];
    
    enum {
        X_index,
        Y_index,
    };

    multitapcombfilter(float _g1, float _g2, float _d1, float _d2)
    {
        gain[X_index] = _g1;
        gain[Y_index] = _g2;
        delayTime[0] = _d1;
        delayTime[1] = _d2;
        for(size_t i = 0; i < 2; i++)
        {
            delay[i].setDelaySize(44100);
            delay[i].setDelayTime(delayTime[i]);
        }       
        x1=y=y1=0;
    }    
    void addTap(float t) {
        delay[0].addTap(t);
        delay[1].addTap(t);
    }
    // X modulation * depth
    // Y modulation * depth
    float Tick(float I, float A=1, float X = 0, float Y=0)
    {
        float x = I;
        y = x + gain[X_index] * x1 - gain[Y_index] * y1;
        x1 = delay[X_index].Tick(x);
        y1 = delay[Y_index].Tick(y);
        return y;
    }
};


struct multitapiircombfilter
{
    multitapdelayline delay;
    float g,x,y,y1;

    multitapiircombfilter(float _g, float _d) 
    {
        delay.setDelaySize(44100);
        delay.setDelayTime(_d);
        g = _g;
        x = y = y1 = 0;
    }
    void addTap(float t) {
        delay.addTap(t);
    }
    float Tick(float I, float A = 1, float X = 0, float Y= 0)
    {
        x = I;
        y = x - g*y1;
        y1= delay.Tick(y);
        return y;
    }
};


struct multitapfircombfilter
{
    multitapdelayline delay;
    float g,x,x1,y;

    multitapfircombfilter(float _g, float _d) 
    {
        delay.setDelaySize(44100);
        delay.setDelayTime(_d);
        g = _g;
        x = y = x1 = 0;
    }
    void addTap(float t) {
        delay.addTap(t);
    }
    float Tick(float I, float A = 1, float X = 0, float Y= 0)
    {
        x = I;
        y = x + g*x1;
        x1= delay.Tick(x);
        return y;
    }
};