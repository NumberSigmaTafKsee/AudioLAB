#pragma once


namespace kfr3
{
    template<typename T>
    struct DelayLine
    {
        kfr::univector<T> delay;
        size_t read_cursor,write_cursor;
        T feedback;
        T mix;
        size_t delayLen;
        T delayTime;

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

        DelayLine() =  default;

        DelayLine(T delay_time, T delay_size = 1) {
            delay.resize(sampleRate*delay_size);
            //memset(delay.data(),0,delay.size()*sizeof(T));
            zeros(delay);
            feedback=0.5;
            mix = 0.5;
            delayTime = delay_time;        
            delayLen  = delay_time * sampleRate;        
            read_cursor  = 0;
            write_cursor = delayLen;
            interpType = Linear;
        }
        T& operator[](size_t i) {
            return delay[i % delay.size()];
        }
        
        void setDelaySize(T size) {
            delay.resize(size);
            memset(delay.data(),0,delay.size()*sizeof(T));
        }
        void reset() {
            read_cursor  = 0;
            write_cursor = delayLen;   
        }
        void setDelayTime(T f) 
        {
            delayTime = f;
            delayLen  = f;
            write_cursor = (read_cursor + delayLen) % delayLen;
            read_cursor  = 0;
            write_cursor = delayLen;
        }
        void setFeedback(T f) {
            feedback = f;
        }
        void setMix(T m) {
            mix = std::fmod(m,1);
        }
        void resize(size_t n) {
            delay.resize(n);
        }
        virtual T Tick(T I, T A = 1, T X = 1, T Y = 1) {
            
            T output;
            size_t n = read_cursor;
            delay.ringbuf_read(read_cursor,output);            
            T d1 = A*(I - Y*output*feedback);                
            T f= d1 - floor(d1); //(int)d1;        
            output = Interpolate(n,f);
            delay.ringbuf_write(write_cursor, I*Y + feedback*output);        
            write_cursor %= delayLen;
            read_cursor %= delayLen;
            return mix*I + (1.0-mix)*output;
        }
        size_t size() { return delay.size(); }    

        T Interpolate(size_t n, T frac)
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
        T NearestNeighborInterpolate(size_t n,T frac)
        {
            int   x  = round(frac);
            T x1 = delay[ (n + x) % delayLen];
            return x1;
        }
        T LinearInterpolate(size_t n,T frac)
        {            
            T x1 = delay[n];
            T x2 = delay[ (n+1) % delayLen];
            //T frac = x1 - (int)x1;
            return x1 + ((x2-x1)*frac);
        }    
        // just cubic stuff from musicdsp
        T CubicInterpolate(size_t inpos,T finpos)
        {            
            T xm1 = delay[(inpos - 1) % delayLen];
            T x0 = delay[inpos + 0];
            T x1  = delay[(inpos + 1) % delayLen];
            T x2  = delay[(inpos + 2) % delayLen];
            T a = (3 * (x0-x1) - xm1 + x2) / 2;
            T b = 2*x1 + xm1 - (5*x0 + x2) / 2;
            T c = (x1 - xm1) / 2;
            return (((a * finpos) + b) * finpos + c) * finpos + x0;
        }
        // just another kind of cubials it might be the same kakaloke really
        inline T Hermite1Interpolate(size_t inpos, T x)
        {            
            T y0 = delay[(inpos - 1) % delayLen];
            T y1 = delay[inpos + 0];
            T y2  = delay[(inpos + 1) % delayLen];
            T y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
            T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);

            return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite2Interpolate(size_t inpos, T x)
        {            
            T y0 = delay[(inpos - 1) % delayLen];
            T y1 = delay[inpos + 0];
            T y2  = delay[(inpos + 1) % delayLen];
            T y3  = delay[(inpos + 2) % delayLen];
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
            T c2 = y0 - y1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite3Interpolate(size_t inpos, T x)
        {                
                T y0 = delay[(inpos - 1) % delayLen];
                T y1 = delay[inpos + 0];
                T y2  = delay[(inpos + 1) % delayLen];
                T y3  = delay[(inpos + 2) % delayLen];
                // 4-point, 3rd-order Hermite (x-form)
                T c0 = y1;
                T c1 = 0.5f * (y2 - y0);
                T y0my1 = y0 - y1;
                T c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
                T c2 = y0my1 + c1 - c3;

                return ((c3 * x + c2) * x + c1) * x + c0;
        }    
        inline T Hermite4Interpolate(size_t inpos, T frac_pos)
        {            
            T xm1 = delay[(inpos - 1) % delayLen];
            T x0 = delay[inpos + 0];
            T x1  = delay[(inpos + 1) % delayLen];
            T x2  = delay[(inpos + 2) % delayLen];
            const T    c     = (x1 - xm1) * 0.5f;
            const T    v     = x0 - x1;
            const T    w     = c + v;
            const T    a     = w + v + (x2 - x0) * 0.5f;
            const T    b_neg = w + a;

            return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
        }
        T Spline3Interpolate(size_t inpos, T x)
        {            
            T L1 = delay[(inpos-1)%delayLen];
            T L0 = delay[inpos];
            T H0 = delay[(inpos + 1) % delayLen];
            T H1 = delay[(inpos + 2) % delayLen];
            return L0 + .5f* x*(H0-L1 + x*(H0 + L0*(-2) + L1 + x*( (H0 - L0)*9 + (L1 - H1)*3 + x*((L0 - H0)*15 + (H1 -  L1)*5 +  x*((H0 - L0)*6 + (L1 - H1)*2 )))));
        }
        T Spline5Interpolate(size_t inpos, T x)
        {
            /* 5-point spline*/
            int nearest_sample = delayLen;

            T p0=delay[(nearest_sample-2) % delayLen];
            T p1=delay[(nearest_sample-1) % delayLen];
            T p2=delay[nearest_sample];
            T p3=delay[(nearest_sample+1) % delayLen];
            T p4=delay[(nearest_sample+2) % delayLen];
            T p5=delay[(nearest_sample+3) % delayLen];

            return p2 + 0.04166666666*x*((p3-p1)*16.0+(p0-p4)*2.0
            + x *((p3+p1)*16.0-p0-p2*30.0- p4
            + x *(p3*66.0-p2*70.0-p4*33.0+p1*39.0+ p5*7.0- p0*9.0
            + x *( p2*126.0-p3*124.0+p4*61.0-p1*64.0- p5*12.0+p0*13.0
            + x *((p3-p2)*50.0+(p1-p4)*25.0+(p5-p0)*5.0)))));
        }
    };
}
