#pragma once
#include "Allocator.hpp"

namespace FX
{
    ///////////////////////////////////////////////////////////////////////////////
    // Compressor
    ///////////////////////////////////////////////////////////////////////////////
    struct MonoCompressor : public MonoFXProcessor
    {
        bool autoTime=true;
        using sample_vector = std::vector<DspFloatType,Allocator::aligned_allocator<DspFloatType,64>>;
        sample_vector x_g, x_l, y_g, y_l, c; // input, output, control
        DspFloatType yL_prev=0;
        DspFloatType sampleRate=44100.0;

        DspFloatType compressionRatio = 10.0f;
        DspFloatType attackTimeMs = 15.0f, releaseTimeMs=100.0f;
        DspFloatType threshold=-24.0f; //dB
        DspFloatType makeUpGain=0.0f;
        DspFloatType currentGain = 0;

        DspFloatType ratio=10.0f,tauAttack=15.0f,tauRelease=100.0f,alphaAttack=15,alphaRelease=100;

        enum {
            PORT_RATIO,
            PORT_ATKTIME,
            PORT_RELTIME,
            PORT_THRESHOLD,
            PORT_MAKEUPGAIN,
            PORT_CURRENTGAIN,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_RATIO: compressionRatio = v; break;
                case PORT_ATKTIME: attackTimeMs = v; break;
                case PORT_RELTIME: releaseTimeMs = v; break;
                case PORT_THRESHOLD: threshold = v; break;
                case PORT_MAKEUPGAIN: makeUpGain = v; break;
                case PORT_CURRENTGAIN: currentGain = v; break;
                default: printf("No port %d\n",port);
            }
        }
        MonoCompressor() : MonoFXProcessor()
        {
			int bufferSize=2048;
			x_g.resize(bufferSize*2);
            x_l.resize(bufferSize*2);
            y_g.resize(bufferSize*2);
            y_l.resize(bufferSize*2);
            c.resize(bufferSize*2);
            memset(x_g.data(),0,x_g.size()*sizeof(DspFloatType));
            memset(y_g.data(),0,y_g.size()*sizeof(DspFloatType));
            memset(x_l.data(),0,x_l.size()*sizeof(DspFloatType));
            memset(y_l.data(),0,y_l.size()*sizeof(DspFloatType));
            memset(c.data(),0,c.size()*sizeof(DspFloatType));
        }
        void computeCompressionGain(size_t bufferSize, DspFloatType * buffer)
        {
            DspFloatType alphaAttack = exp(-1.0f / (0.001f * sampleRate * attackTimeMs));
            DspFloatType alphaRelease = exp(-1.0f / (0.001f * sampleRate * releaseTimeMs));
            DspFloatType yl_avg = 0.0f;
            #pragma omp simd aligned(buffer)
            for (int i = 0; i < bufferSize; ++i)
            {
                // Level detection- estimate level using peak detector
                if (fabs(buffer[i]) < 0.000001f) x_g[i] = -120;
                else x_g[i] = 20 * log10(fabs(buffer[i]));

                // Gain computer- static apply input/output curve
                if (x_g[i] >= threshold)
                    y_g[i] = threshold + (x_g[i] - threshold) / compressionRatio;
                else
                    y_g[i] = x_g[i];

                x_l[i] = x_g[i] - y_g[i];

                // Ballistics- smoothing of the gain 
                if (x_l[i] > yL_prev)
                    y_l[i] = alphaAttack * yL_prev + (1 - alphaAttack) * x_l[i];
                else
                    y_l[i] = alphaRelease * yL_prev + (1 - alphaRelease) * x_l[i];

                // Accumulate averaged gain for the whole buffer (for GUI display)
                yl_avg += y_l[i];

                // find control
                c[i] = pow(10.0f, (makeUpGain - y_l[i]) / 20.0f);
                yL_prev = y_l[i];
            }

            yl_avg /= bufferSize;
            currentGain = pow(10.0f, -yl_avg / 20.0f);    
        }
        void compressor(size_t n, DspFloatType * input)
        {
            alphaAttack = exp(-1/(0.001 * sampleRate * tauAttack));
            alphaRelease= exp(-1/(0.001 * sampleRate * tauRelease));
            #pragma omp simd aligned(input)
            for (int i = 0 ; i <  n; ++i)
            {
                //Level detection- estimate level using peak detector
                if (fabs(input[i]) < 0.000001) x_g[i] =-120;
                else x_g[i] =20*log10(fabs(input[i]));
                //Gain computer- static apply input/output curve
                if (x_g[i] >= threshold) y_g[i] = threshold+ (x_g[i] - threshold) / ratio;
                else y_g[i] = x_g[i];
                x_l[i] = x_g[i] - y_g[i];
                //Ballistics- smoothing of the gain 
                if (x_l[i]>yL_prev)  y_l[i]=alphaAttack * yL_prev+(1 - alphaAttack ) * x_l[i] ; 
                else				 y_l[i]=alphaRelease* yL_prev+(1 - alphaRelease) * x_l[i] ;
                //find control
                c[i] = pow(10,(makeUpGain - y_l[i])/20);
                yL_prev=y_l[i];
            }
        }
        void ProcessBlock (size_t n, DspFloatType * inputs, DspFloatType * outputs)
        {
            Undenormal denormal;
            if (threshold < 0.0f)
            {        
                // compression : calculates the control voltage
                //computeCompressionGain(n,input);
                compressor(n,inputs);
                // apply control voltage to the audio signal        
                #pragma omp simd aligned(inputs,outputs)
                for (int i = 0; i < n; ++i)
                {
                    DspFloatType cv = c[i];
                    outputs[i] = cv * inputs[i];            
                }        
            }
        }
    };

    struct StereoCompressor : public StereoFXProcessor
    {
		using sample_vector = std::vector<DspFloatType,Allocator::aligned_allocator<DspFloatType,64>>;
        sample_vector x_g, x_l,y_g, y_l,c;// input, output, control
        sample_vector inputBuffer;

        bool autoTime,compressorSwitch = true;
        DspFloatType ratio,threshold,makeUpGain,tauAttack,tauRelease,alphaAttack,alphaRelease,yL_prev;
        int nhost;    
        DspFloatType sampleRate=44100.0;
        
        enum {
            
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
            
            }
        }
        StereoCompressor() : StereoFXProcessor()
        {                        
            int bufferSize=2048;
            // Allocate a lot of dynamic memory here
            x_g.resize(bufferSize*2);
            x_l.resize(bufferSize*2);
            y_g.resize(bufferSize*2);
            y_l.resize(bufferSize*2);
            c.resize(bufferSize*2);
            inputBuffer.resize(bufferSize*2);        
            yL_prev=0;
            autoTime = false;
            compressorSwitch = true;
            ratio = 4;
            threshold=-6;
            resetAll();
        }
        void resetAll()
        {
            tauAttack=0;tauRelease = 0;
            alphaAttack=0;alphaRelease = 0;
            //threshold = 0;
            //ratio= 1;
            makeUpGain= 0;
            yL_prev=0;
            memset(x_g.data(),0,x_g.size()*sizeof(DspFloatType));
            memset(y_g.data(),0,y_g.size()*sizeof(DspFloatType));
            memset(x_l.data(),0,x_l.size()*sizeof(DspFloatType));
            memset(y_l.data(),0,y_l.size()*sizeof(DspFloatType));
            memset(c.data(),0,c.size()*sizeof(DspFloatType));
        }
        void compressor(size_t n, DspFloatType * buffer)
        {
            alphaAttack = exp(-1/(0.001 * sampleRate * tauAttack));
            alphaRelease= exp(-1/(0.001 * sampleRate * tauRelease));
            
            #pragma omp simd aligned(buffer)
            for (int i = 0 ; i < n; ++i)
            {
                //Level detection- estimate level using peak detector
                if (fabs(buffer[i]) < 0.000001) x_g[i] =-120;
                else x_g[i] =20*log10(fabs(buffer[i]));
                //Gain computer- static apply input/output curve
                if (x_g[i] >= threshold) y_g[i] = threshold+ (x_g[i] - threshold) / ratio;
                else y_g[i] = x_g[i];
                x_l[i] = x_g[i] - y_g[i];
                //Ballistics- smoothing of the gain 
                if (x_l[i]>yL_prev)  y_l[i]=alphaAttack * yL_prev+(1 - alphaAttack ) * x_l[i] ; 
                else				 y_l[i]=alphaRelease* yL_prev+(1 - alphaRelease) * x_l[i] ;
                //find control
                c[i] = pow(10,(makeUpGain - y_l[i])/20);
                yL_prev=y_l[i];
            }
        }
        void ProcessBlock (size_t numSamples, DspFloatType ** inputs, DspFloatType ** outputs)
        {
            if (compressorSwitch == true)
            {           
                memset(inputBuffer.data(),0,inputBuffer.size()*sizeof(DspFloatType));
                DspFloatType * inbuf = inputBuffer.data();
                #pragma omp simd aligned(inbuf,inputs,outputs)
                for(size_t i = 0; i < numSamples; i++)
                {
                    outputs[0][i] = inputs[0][i];
                    outputs[1][i] = inputs[1][i];
                    inbuf[i] = 0.5*(inputs[0][i]+inputs[1][i]);
                } 
                
                #pragma omp simd aligned(inbuf,inputs,outputs)
                for(size_t m = 0; m < 2; m++)
                if ( (threshold< 0) )
                {                
                    compressor(numSamples,inbuf);
                    // apply control voltage to the audio signal
                    for (int i = 0 ; i < numSamples ; ++i)
                    {
                        outputs[0][i] *= c[i];
                        outputs[1][i] *= c[i];                    
                        inbuf[i*2] = 0.5 * (outputs[0][i] + outputs[1][i]);
                    }            
                }
            }        
        }
        void ProcessInplace(size_t n, DspFloatType ** buffer) {
            ProcessBlock(n,buffer,buffer);
        }
    };
}
