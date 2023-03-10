#pragma once
#include "vcvfft.hpp"
#include "SndFile.hpp"

namespace DSP::Convolution
{
    struct PffftConvolver : public StereoFXProcessor
    {
        StereoConvolution * conv;
        
        
        PffftConvolver(const char * impulse_file, int BufferSize)    
        : StereoFXProcessor()
        {
            //
            sample_vector<float> lexiconL,lexiconR;
            SndFileReaderFloat sndfile(impulse_file);
            
            std::cout << sndfile.channels() << std::endl;
            std::cout << sndfile.frames() << std::endl;
            std::cout << sndfile.samplerate() << std::endl;
            std::vector<float> channelv(sndfile.frames()*sndfile.channels());
            
            sndfile >> channelv;

            lexiconL.resize(sndfile.frames());
            lexiconR.resize(sndfile.frames());

            for(size_t i = 0; i < sndfile.frames(); i++)
            {        
                lexiconL[i] = channelv[i*2];        
                lexiconR[i] = channelv[i*2+1];
            }

            //convolverL.init(BufferSize,lexiconL.data(),lexiconL.size());
            //convolverR.init(BufferSize,lexiconR.data(),lexiconR.size());
            conv = new StereoConvolution(lexiconL.data(), lexiconR.data(), lexiconL.size(),BufferSize);        
        }
        ~PffftConvolver()
        {
            if(conv) delete conv;
        }
        void ProcessBlock(size_t framesPerBuffer, float ** in, float ** out)
        {
            conv->Process(framesPerBuffer,in,out);        
        }
        void InplaceProcess(size_t n, float ** buffer)
        {
            ProcessBlock(n,buffer,buffer);
        }
    };
}