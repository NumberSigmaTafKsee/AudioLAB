
#include "aubio_fvec.hpp"
#include <rubberband/rubberband-c.h>

extern RubberBandOptions aubio_get_rubberband_opts(const char_t *mode);

template<typename T>
struct PitchShift
{
  size_t samplerate;              /**< samplerate */
  size_t hopsize;                 /**< hop size */
  T pitchscale;                   /**< pitch scale */

  RubberBandState rb;
  RubberBandOptions rboptions;

  PitchShift(const char_t * mode,T transpose, size_t hopsize, size_t samplerate)
  {
    
    this->samplerate = samplerate;
    this->hopsize = hopsize;
    this->pitchscale = 1.;
    this->rb = NULL;

    if ((int)hopsize <= 0) {
      std::cerr << "pitchshift: hop_size should be >= 0, got " << hopsize << std::endl;    
    }
    if ((int)samplerate <= 0) {
      std::cerr << "pitchshift: samplerate should be >= 0, got " << samplerate << std::endl;    
    }

    rboptions = get_rubberband_opts(mode);
    if (rboptions < 0) {
      std::cerr << "pitchshift: unknown pitch shifting method " << mode << std::endl;    
    }

    rb = rubberband_new(samplerate, 1, p->rboptions, 1., p->pitchscale);
    rubberband_set_max_process_size(p->rb, p->hopsize);
    //rubberband_set_debug_level(p->rb, 10);

    if (set_transpose(p, transpose)) goto beach;

  #if 1
    // warm up rubber band
    unsigned int latency = std::max(p->hopsize, rubberband_get_latency(p->rb));
    int available = rubberband_available(p->rb);
    FVec<T> zeros(p->hopsize);
    zeros.fill((T)0.0);
    while (available <= (int)latency) {
      rubberband_process(rb,
          (const float* const*)&(zeros.data()), hopsize, 0);
      available = rubberband_available(rb);
    }  
  #endif  
  }
  ~PitchShift()
  {
    if (rb) rubberband_delete(rb);
  }  


  size_t get_latency () {
    return rubberband_get_latency();
  }

  bool set_pitchscale (T pitchscale)
  {
    if (pitchscale >= 0.25  && pitchscale <= 4.) {
      this->pitchscale = pitchscale;
      rubberband_set_pitch_scale(rb, pitchscale);
      return true;      
    } else {
      std::cerr << "pitchshift: could not set pitchscale to " << pitchscale << " should be in the range [0.25, 4.].\n";
      return false;
    }
  }

  T get_pitchscale ()
  {
    return this->pitchscale;
  }

  bool set_transpose(T transpose)
  {
    if (transpose >= -24. && transpose <= 24.) {
      T pitchscale = POW(2., transpose / 12.);
      return set_pitchscale(p, pitchscale);
    } else {
      std::cerr << "pitchshift: could not set transpose to " << transpose
          << " should be in the range [-24; 24.].\n";
      return false;
    }
  }

  T get_transpose()
  {
    return 12. * std::log(p->pitchscale) / std::log(2.0);
  }

  void ProcessBlock(const FVec<T> & in, FVec<T> & out)
  {    
    rubberband_process(rb, (const float* const*)&(in.data()), hopsize, 0);
    if (rubberband_available(rb) >= (int)hopsize) {
      rubberband_retrieve(rb, (float* const*)&(out.data()), hopsize);
    } else {
      std:: cerr << "pitchshift: catching up with zeros" <<
          ", only " <<  rubberband_available(rb)  << "available, needed: " << hopsize
          << "current pitchscale: " << pitchscale << std::endl;
      out.fill((T)0.0); 
    }
  }
};