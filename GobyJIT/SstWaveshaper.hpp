#pragma once

struct SstWaveShaper
{
    sst::waveshapers::QuadWaveshaperState    waState;
    sst::waveshapers::QuadWaveshaperPtr      pWaveShaper;
    float Distortion = 0;
    
    SstWaveShaper(sst::waveshapers::WaveshaperType type = sst::waveshapers::WaveshaperType::wst_fuzz)
    {
        float R[sst::waveshapers::n_waveshaper_registers];
        sst::waveshapers::initializeWaveshaperRegister(type, R);
        for (int i = 0; i < sst::waveshapers::n_waveshaper_registers; ++i)
        {
            waState.R[i] = _mm_set1_ps(R[i]);
        }
        waState.init = _mm_cmpneq_ps(_mm_setzero_ps(), _mm_setzero_ps());

        pWaveShaper = sst::waveshapers::GetQuadWaveshaper(sst::waveshapers::WaveshaperType::wst_fuzz);
    }
    ~SstWaveShaper()
    {

    }

    void Process(size_t framesPerBuffer, float ** buffer)
    {
        __m128 gain = _mm_set1_ps(pow(10,(Distortion/20.0)));
        for(size_t i = 0; i < framesPerBuffer; i+=4)
        {    
            __m128 x = _mm_load_ps(&buffer[0][i]);                
            x = pWaveShaper(&waState,x,gain);        
            _mm_store_ps(&buffer[0][i],x);

            x = _mm_load_ps(&buffer[1][i]);
            x = pWaveShaper(&waState,x,gain);        
            _mm_store_ps(&buffer[1][i],x);
        }
        
    }
};
