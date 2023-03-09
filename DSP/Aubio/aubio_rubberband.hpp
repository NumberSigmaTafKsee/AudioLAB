
#include "aubio_fvec.hpp"
#include <rubberband/rubberband-c.h>
#include <pthread.h>

#define MIN_STRETCH_RATIO 0.025
#define MAX_STRETCH_RATIO 40.
#define PATH_MAX 1024

char** split_str(const char* str, const char sep) {
  char** result = 0;
  uint32_t count = 0;
  char input[PATH_MAX];
  char* in_ptr = input;
  char* last_sep = 0;
  char delim[2]; delim[0] = sep; delim[1] = 0;

  strncpy(input, str, PATH_MAX);
  input[PATH_MAX - 1] = '\0';

  // count number of elements
  while (*in_ptr) {
    if (sep == *in_ptr) {
      count++;
      last_sep = in_ptr;
    }
    in_ptr++;
  }
  // add space for trailing token.
  count += last_sep < (input + strlen(input) - 1);
  count++;

  result = (char**)malloc(count * sizeof(char*));
  if (result) {
    uint32_t idx = 0;
    char* params = strtok(input, delim);
    while (params) {
      // make sure we don't got in the wild
      if (idx >= count)
        break;
      *(result + idx++) = strdup(params);
      params = strtok(0, delim);
    }
    // add null string at the end if needed
    if (idx < count - 1)
      *(result + idx) = 0;
  }
  return result;
}
RubberBandOptions get_rubberband_opts(const char *mode)
{
  RubberBandOptions rboptions = RubberBandOptionProcessRealTime;

  if ( strcmp(mode,"crispness:0") == 0 ) {
    rboptions |= RubberBandOptionTransientsSmooth;
    rboptions |= RubberBandOptionWindowLong;
    rboptions |= RubberBandOptionPhaseIndependent;
  } else if ( strcmp(mode, "crispness:1") == 0 ) {
    rboptions |= RubberBandOptionDetectorSoft;
    rboptions |= RubberBandOptionTransientsSmooth;
    rboptions |= RubberBandOptionWindowLong;
    rboptions |= RubberBandOptionPhaseIndependent;
  } else if ( strcmp(mode, "crispness:2") == 0 ) {
    rboptions |= RubberBandOptionTransientsSmooth;
    rboptions |= RubberBandOptionPhaseIndependent;
  } else if ( strcmp(mode, "crispness:3") == 0 ) {
    rboptions |= RubberBandOptionTransientsSmooth;
  } else if ( strcmp(mode, "crispness:4") == 0 ) {
    // same as "default"
  } else if ( strcmp(mode, "crispness:5") == 0 ) {
    rboptions |= RubberBandOptionTransientsCrisp;
  } else if ( strcmp(mode, "crispness:6") == 0 ) {
    rboptions |= RubberBandOptionTransientsCrisp;
    rboptions |= RubberBandOptionWindowShort;
    rboptions |= RubberBandOptionPhaseIndependent;
  } else if ( strcmp(mode, "default") == 0 ) {
    // nothing to do
  } else {
    // attempt to parse a list of options, separated with ','
    char **params = split_str(mode, ':');
    uint32_t i = 0;
    if (!params || !params[0]) {      
      rboptions = -1;
    }
    while (*(params + i) != NULL) {
      if ( strcmp(params[i], "ProcessOffline" ) == 0 )        {
             rboptions = RubberBandOptionProcessOffline;        
        rboptions = -1;
      }
      else if ( strcmp(params[i], "ProcessRealTime" ) == 0 )       rboptions |= RubberBandOptionProcessRealTime;
      else if ( strcmp(params[i], "StretchElastic" ) == 0 )        rboptions |= RubberBandOptionStretchElastic;
      else if ( strcmp(params[i], "StretchPrecise" ) == 0 )        rboptions |= RubberBandOptionStretchPrecise;
      else if ( strcmp(params[i], "TransientsCrisp" ) == 0 )       rboptions |= RubberBandOptionTransientsCrisp;
      else if ( strcmp(params[i], "TransientsMixed" ) == 0 )       rboptions |= RubberBandOptionTransientsMixed;
      else if ( strcmp(params[i], "TransientsSmooth" ) == 0 )      rboptions |= RubberBandOptionTransientsSmooth;
      else if ( strcmp(params[i], "DetectorCompound" ) == 0 )      rboptions |= RubberBandOptionDetectorCompound;
      else if ( strcmp(params[i], "DetectorPercussive" ) == 0 )    rboptions |= RubberBandOptionDetectorPercussive;
      else if ( strcmp(params[i], "DetectorSoft" ) == 0 )          rboptions |= RubberBandOptionDetectorSoft;
      else if ( strcmp(params[i], "PhaseLaminar" ) == 0 )          rboptions |= RubberBandOptionPhaseLaminar;
      else if ( strcmp(params[i], "PhaseIndependent" ) == 0 )      rboptions |= RubberBandOptionPhaseIndependent;
      else if ( strcmp(params[i], "ThreadingAuto" ) == 0 )         rboptions |= RubberBandOptionThreadingAuto;
      else if ( strcmp(params[i], "ThreadingNever" ) == 0 )        rboptions |= RubberBandOptionThreadingNever;
      else if ( strcmp(params[i], "ThreadingAlways" ) == 0 )       rboptions |= RubberBandOptionThreadingAlways;
      else if ( strcmp(params[i], "WindowStandard" ) == 0 )        rboptions |= RubberBandOptionWindowStandard;
      else if ( strcmp(params[i], "WindowShort" ) == 0 )           rboptions |= RubberBandOptionWindowShort;
      else if ( strcmp(params[i], "WindowLong" ) == 0 )            rboptions |= RubberBandOptionWindowLong;
      else if ( strcmp(params[i], "SmoothingOff" ) == 0 )          rboptions |= RubberBandOptionSmoothingOff;
      else if ( strcmp(params[i], "SmoothingOn" ) == 0 )           rboptions |= RubberBandOptionSmoothingOn;
      else if ( strcmp(params[i], "FormantShifted" ) == 0 )        rboptions |= RubberBandOptionFormantShifted;
      else if ( strcmp(params[i], "FormantPreserved" ) == 0 )      rboptions |= RubberBandOptionFormantPreserved;
      else if ( strcmp(params[i], "PitchHighSpeed" ) == 0 )        rboptions |= RubberBandOptionPitchHighSpeed;
      else if ( strcmp(params[i], "PitchHighQuality" ) == 0 )      rboptions |= RubberBandOptionPitchHighQuality;
      else if ( strcmp(params[i], "PitchHighConsistency" ) == 0 )  rboptions |= RubberBandOptionPitchHighConsistency;
      else if ( strcmp(params[i], "ChannelsApart" ) == 0 )         rboptions |= RubberBandOptionChannelsApart;
      else if ( strcmp(params[i], "ChannelsTogether" ) == 0 )      rboptions |= RubberBandOptionChannelsTogether;
      else {        
        rboptions = -1;
      }      
      i++;
    } 
    free(params);
  }
  return rboptions;
}

template<typename T>
struct PitchShift
{
  size_t samplerate;              /**< samplerate */
  size_t hopsize;                 /**< hop size */
  T pitchscale;                   /**< pitch scale */

  RubberBandState rb;
  RubberBandOptions rboptions;

  PitchShift(const char * mode,T transpose, size_t hopsize, size_t samplerate)
  {
    
    this->samplerate = samplerate;
    this->hopsize = hopsize;
    this->pitchscale = 1.;
    this->rb = NULL;

    rboptions = get_rubberband_opts(mode);
    if (rboptions < 0) {      
      throw std::runtime_error("unknown pitch shifting method");
    }

    rb = rubberband_new(samplerate, 1, rboptions, 1., pitchscale);
    rubberband_set_max_process_size(rb, hopsize);
    //rubberband_set_debug_level(p->rb, 10);

    setTranspose(transpose);

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


  size_t getLatency () {
    return rubberband_get_latency();
  }

  void setPitchScale (T pitchscale)
  {
    if (pitchscale >= 0.25  && pitchscale <= 4.) {
      this->pitchscale = pitchscale;
      rubberband_set_pitch_scale(rb, pitchscale);      
    } else {      
      throw std::runtime_error("could not set pitchscale");
    }
  }

  T getPitchScale ()
  {
    return this->pitchscale;
  }

  void setTranspose(T transpose)
  {
    if (transpose >= -24. && transpose <= 24.) {
      T pitchscale = std::pow(2., transpose / 12.);
      setPitchScale(p, pitchscale);
    } else {
      throw std::runtime_error("could not set transpose");
    }
  }

  T getTranspose()
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

template<typename T>
struct TimeStretch
{
  uint32_t samplerate;              /**< samplerate */
  uint32_t hopsize;                 /**< hop size */
  T stretchratio;            /**< time ratio */
  T pitchscale;              /**< pitch scale */

  RubberBandState rb;
  RubberBandOptions rboptions;

  TimeStretch(const char * mode, T stretchratio, uint32_t hopsize, uint32_t samplerate)
  {    
    this->hopsize = hopsize;
    this->pitchscale = 1.;
    
    if (stretchratio <= MAX_STRETCH_RATIO && stretchratio >= MIN_STRETCH_RATIO) {
      this->stretchratio = stretchratio;
    } else {
      throw std::runtime_error("stretchratio should be in the range")
    }

    this->rboptions = aubio_get_rubberband_opts(mode);
    if (p->rboptions < 0) {      
      throw std::runtime_error("unknown time stretching method");
    }

    this->rb = rubberband_new(samplerate, 1, rboptions, stretchratio, pitchscale);
    if (!this->rb) throw std::runtime_error("Error allocating rubberband");

    this->samplerate = samplerate;  
  }
  ~TimeStretch()
  {
    if(rb) rubberband_delete(rb);
  }
  uint32_t getSampleRate() const { return samplerate; }
  uint32_t getLatency() { return rubberband_get_latency(rb); }
  T getStretch() const { return stretchratio; }
  T getTranspose() const {   return 12. * std::log(p->pitchscale) / std::log(2.0); }
  
  void setStretch(const T& stretch)
  {
    if (!rb) {      
      throw std::runtime_error("rubberband not created");
    }
    if (stretch >= MIN_STRETCH_RATIO && stretch <= MAX_STRETCH_RATIO) {
      this->stretchratio = stretch;
      rubberband_set_time_ratio(b, 1./stretchratio);      
    } else {
      throw std::runtime_error("could not set stretch ratio");
    }
  }
  void setPitchScale (T pitchscale)
  {
    if (!rb) {
      throw std::runtime_error("rubberband not created");
    }
    if (pitchscale >= 0.0625  && pitchscale <= 4.) {
      p->pitchscale = pitchscale;
      rubberband_set_pitch_scale(p->rb, p->pitchscale);      
    } else {
      throw std::runtime_error("could not set pitch scale");
    }
  }
  void setTranspose(const T& transpose)
  {
    if (transpose >= -24. && transpose <= 24.) {
      smpl_t pitchscale = POW(2., transpose / 12.);
      setPitchScale(pitchscale);
    } else {
      throw std::runtime_error("could not set transpose");
    }
  }
};


