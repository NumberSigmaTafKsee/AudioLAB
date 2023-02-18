
#pragma once

template<typename Processor>
struct Voice
{
    Processor * proc;

    Voice(Processor * p) : proc(p) {

    }
    ~Voice() {
        delete proc;
    }

};

template<typename Processor>
struct Voices {
    std::vector<Voices<Processor>*> voices;
    std::vector<std::vector<DspFloatType>> temp;

    Voices(size_t nVoices) {
        voices.resize(nVoices);
        temp.resize(nVoices);
        for(size_t i = 0; i < nVoices; i++) {
            voices[i] = new Processor;
            tenp[i].resize(1024);
        }
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType *out) {
        temp.resize(n); 
        #pragma omp parallel for         
        for(size_t i = 0; i < voices.size(); i++ ) {
            voices[i]->proc->ProcessBlock(n,in,temp.data());
        }
        #pragma omp simd
        for(size_t i = 0; i < voices.size(); i++) {
            for(size_t j = 0; j < n; j++)
                out[j] += temp[i][j];
        }
        #pragma omp simd
        for(size_t i = 0; i < voices.size(); i++) {
            out[j] /= (DspFloatType)voices.size();
        }
    }
};