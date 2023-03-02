#include <cmath>
#include <cstdint>
#include <rubberband/rubberband-c.h>

struct PitchShiftRubberband
{
    uint32_t samplerate;              /**< samplerate */
    uint32_t hopsize;                 /**< hop size */
    smpl_t pitchscale;              /**< pitch scale */

    RubberBandState rb;
    RubberBandOptions rboptions;


    PitchShiftRubberband(float sr, int hop, float pcale=1.0f)
    {        
        samplerate = sr;
        hopsize = hop;
        pitchscale = pscale;

        rb = NULL;
        if ((int32_t)hopsize <= 0) {
            printf("pitchshift: hop_size should be >= 0, got %d\n", hopsize);
            goto beach;
        }
        if ((int32_t)samplerate <= 0) {
            printf("pitchshift: samplerate should be >= 0, got %d\n", samplerate);
            goto beach;
        }

        rboptions = aubio_get_rubberband_opts(mode);
        if (rboptions < 0) {
            printf("pitchshift: unknown pitch shifting method %s\n", mode);
            goto beach;
        }
        
        rb = rubberband_new(samplerate, 1, rboptions, 1., pitchscale);
        rubberband_set_max_process_size(rb, hopsize);
        //rubberband_set_debug_level(rb, 10);

        if (setTranspose(transpose)) 
        {
            printf("transpose: failed\n");
            goto beach;
        }

        
        // warm up rubber band
        unsigned int latency = std::max(hopsize, rubberband_get_latency(rb));
        int available = rubberband_available(rb);
        
        sample_vector<float> zerovec(hopsize);
        memset(zerovec.data(),0,hopsize);
        while (available <= (int)latency) {
            rubberband_process(rb,(const float* const*)(zerovec.data()), hopsize, 0);
            available = rubberband_available(rb);
        }                
        

        beach:
            printf("Error initializing rubberband\n");
            exit(-1);        
    }
    ~PitchShiftRubberband()
    {
        if(rb)     rubberband_delete(rb);
    }

    uint32_t getLatency()
    {
          return rubberband_get_latency(rb);
    }

    bool setPitchScale (float pscale)
    {
        if (pscale >= 0.25  && pscale <= 4.) {
            pitchscale = pscale;
            rubberband_set_pitch_scale(rb, pcale);
            return true;
        } else {
            printf("pitchshift: could not set pitchscale to '%f',"
                " should be in the range [0.25, 4.].\n", pitchscale);
            return false;
        }
    }    
    bool setTranspose(float transpose)
    {
        if (transpose >= -24. && transpose <= 24.) {            
            return setPitchScale(std::pow(2., transpose / 12.));
        } else {
            printf("pitchshift: could not set transpose to '%f',"
                " should be in the range [-24; 24.].\n", transpose);
            return false;
        }
    }
    void ProcessBlock(size_t n, float * in, float * out)
    {
        // third parameter is always 0 since we are never expecting a final frame
        rubberband_process(rb, (const float* const*)&in hopsize, 0);
        if (rubberband_available(rb) >= (int)hopsize) {
            rubberband_retrieve(rb, (float* const*)&out, hopsize);
        } else {
            printf("pitchshift: catching up with zeros"
                ", only %d available, needed: %d, current pitchscale: %f\n",
                rubberband_available(p->rb), p->hopsize, p->pitchscale);
            memset(out,0,n*sizeof(float));
        }
    }
    void InplaceProcess(size_t n, float * buffer) {
        ProcessBlock(n,buffer,buffer);
    }
};