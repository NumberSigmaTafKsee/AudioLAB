/***************************************************/
/*! \class Delay
    \brief STK non-interpolating delay line class.
    This class implements a non-interpolating digital delay-line.  If
    the delay and maximum length are not specified during
    instantiation, a fixed maximum length of 4095 and a delay of zero
    is set.
    
    A non-interpolating delay line is typically used in fixed
    delay-length applications, such as for reverberation.
    by Perry R. Cook and Gary P. Scavone, 1995--2021.
*/
/***************************************************/

class Delay : public Filter
{
public:

  //! The default constructor creates a delay-line with maximum length of 4095 samples and zero delay.
  /*!
    An StkError will be thrown if the delay parameter is less than
    zero, the maximum delay parameter is less than one, or the delay
    parameter is greater than the maxDelay value.
   */
  Delay( unsigned long delay = 0, unsigned long maxDelay = 4095 );

  //! Class destructor.
  ~Delay();

  //! Get the maximum delay-line length.
  unsigned long getMaximumDelay( void ) { return inputs_.size() - 1; };

  //! Set the maximum delay-line length.
  /*!
    This method should generally only be used during initial setup
    of the delay line.  If it is used between calls to the tick()
    function, without a call to clear(), a signal discontinuity will
    likely occur.  If the current maximum length is greater than the
    new length, no memory allocation change is made.
  */
  void setMaximumDelay( unsigned long delay );

  //! Set the delay-line length.
  /*!
    The valid range for \e delay is from 0 to the maximum delay-line length.
  */
  void setDelay( unsigned long delay );

  //! Return the current delay-line length.
  unsigned long getDelay( void ) const { return delay_; };

  //! Return the value at \e tapDelay samples from the delay-line input.
  /*!
    The tap point is determined modulo the delay-line length and is
    relative to the last input value (i.e., a tapDelay of zero returns
    the last input value).
  */
  StkFloat tapOut( unsigned long tapDelay );

  //! Set the \e value at \e tapDelay samples from the delay-line input.
  void tapIn( StkFloat value, unsigned long tapDelay );

  //! Sum the provided \e value into the delay line at \e tapDelay samples from the input.
  /*!
    The new value is returned.  The tap point is determined modulo
    the delay-line length and is relative to the last input value
    (i.e., a tapDelay of zero sums into the last input value).
  */
  StkFloat addTo( StkFloat value, unsigned long tapDelay );

  //! Return the last computed output value.
  StkFloat lastOut( void ) const { return lastFrame_[0]; };

  //! Return the value that will be output by the next call to tick().
  /*!
    This method is valid only for delay settings greater than zero!
   */
  StkFloat nextOut( void ) { return inputs_[outPoint_]; };

  //! Calculate and return the signal energy in the delay-line.
  StkFloat energy( void ) const;

  //! Input one sample to the filter and return one output.
  StkFloat tick( StkFloat input );

  //! Take a channel of the StkFrames object as inputs to the filter and replace with corresponding outputs.
  /*!
    The StkFrames argument reference is returned.  The \c channel
    argument must be less than the number of channels in the
    StkFrames argument (the first channel is specified by 0).
    However, range checking is only performed if _STK_DEBUG_ is
    defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  StkFrames& tick( StkFrames& frames, unsigned int channel = 0 );

  //! Take a channel of the \c iFrames object as inputs to the filter and write outputs to the \c oFrames object.
  /*!
    The \c iFrames object reference is returned.  Each channel
    argument must be less than the number of channels in the
    corresponding StkFrames argument (the first channel is specified
    by 0).  However, range checking is only performed if _STK_DEBUG_
    is defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  StkFrames& tick( StkFrames& iFrames, StkFrames &oFrames, unsigned int iChannel = 0, unsigned int oChannel = 0 );

protected:

  unsigned long inPoint_;
  unsigned long outPoint_;
  unsigned long delay_;
};

inline StkFloat Delay :: tick( StkFloat input )
{
  inputs_[inPoint_++] = input * gain_;

  // Check for end condition
  if ( inPoint_ == inputs_.size() )
    inPoint_ = 0;

  // Read out next value
  lastFrame_[0] = inputs_[outPoint_++];

  if ( outPoint_ == inputs_.size() )
    outPoint_ = 0;

  return lastFrame_[0];
}

inline StkFrames& Delay :: tick( StkFrames& frames, unsigned int channel )
{
#if defined(_STK_DEBUG_)
  if ( channel >= frames.channels() ) {
    oStream_ << "Delay::tick(): channel and StkFrames arguments are incompatible!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }
#endif

  StkFloat *samples = &frames[channel];
  unsigned int hop = frames.channels();
  for ( unsigned int i=0; i<frames.frames(); i++, samples += hop ) {
    inputs_[inPoint_++] = *samples * gain_;
    if ( inPoint_ == inputs_.size() ) inPoint_ = 0;
    *samples = inputs_[outPoint_++];
    if ( outPoint_ == inputs_.size() ) outPoint_ = 0;
  }

  lastFrame_[0] = *(samples-hop);
  return frames;
}

inline StkFrames& Delay :: tick( StkFrames& iFrames, StkFrames& oFrames, unsigned int iChannel, unsigned int oChannel )
{
#if defined(_STK_DEBUG_)
  if ( iChannel >= iFrames.channels() || oChannel >= oFrames.channels() ) {
    oStream_ << "Delay::tick(): channel and StkFrames arguments are incompatible!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }
#endif

  StkFloat *iSamples = &iFrames[iChannel];
  StkFloat *oSamples = &oFrames[oChannel];
  unsigned int iHop = iFrames.channels(), oHop = oFrames.channels();
  for ( unsigned int i=0; i<iFrames.frames(); i++, iSamples += iHop, oSamples += oHop ) {
    inputs_[inPoint_++] = *iSamples * gain_;
    if ( inPoint_ == inputs_.size() ) inPoint_ = 0;
    *oSamples = inputs_[outPoint_++];
    if ( outPoint_ == inputs_.size() ) outPoint_ = 0;
  }

  lastFrame_[0] = *(oSamples-oHop);
  return iFrames;
}


/***************************************************/
/*! \class DelayA
    \brief STK allpass interpolating delay line class.
    This class implements a fractional-length digital delay-line using
    a first-order allpass filter.  If the delay and maximum length are
    not specified during instantiation, a fixed maximum length of 4095
    and a delay of 0.5 is set.
    An allpass filter has unity magnitude gain but variable phase
    delay properties, making it useful in achieving fractional delays
    without affecting a signal's frequency magnitude response.  In
    order to achieve a maximally flat phase delay response, the
    minimum delay possible in this implementation is limited to a
    value of 0.5.
    by Perry R. Cook and Gary P. Scavone, 1995--2021.
*/
/***************************************************/

class DelayA : public Filter
{
public:

  //! Default constructor creates a delay-line with maximum length of 4095 samples and delay = 0.5.
  /*!
    An StkError will be thrown if the delay parameter is less than
    zero, the maximum delay parameter is less than one, or the delay
    parameter is greater than the maxDelay value.
   */  
  DelayA( StkFloat delay = 0.5, unsigned long maxDelay = 4095 );

  //! Class destructor.
  ~DelayA();

  //! Clears all internal states of the delay line.
  void clear( void );

  //! Get the maximum delay-line length.
  unsigned long getMaximumDelay( void ) { return inputs_.size() - 1; };
  
  //! Set the maximum delay-line length.
  /*!
    This method should generally only be used during initial setup
    of the delay line.  If it is used between calls to the tick()
    function, without a call to clear(), a signal discontinuity will
    likely occur.  If the current maximum length is greater than the
    new length, no memory allocation change is made.
  */
  void setMaximumDelay( unsigned long delay );

  //! Set the delay-line length
  /*!
    The valid range for \e delay is from 0.5 to the maximum delay-line length.
  */
  void setDelay( StkFloat delay );

  //! Return the current delay-line length.
  StkFloat getDelay( void ) const { return delay_; };

  //! Return the value at \e tapDelay samples from the delay-line input.
  /*!
    The tap point is determined modulo the delay-line length and is
    relative to the last input value (i.e., a tapDelay of zero returns
    the last input value).
  */
  StkFloat tapOut( unsigned long tapDelay );

  //! Set the \e value at \e tapDelay samples from the delay-line input.
  void tapIn( StkFloat value, unsigned long tapDelay );

  //! Return the last computed output value.
  StkFloat lastOut( void ) const { return lastFrame_[0]; };

  //! Return the value which will be output by the next call to tick().
  /*!
    This method is valid only for delay settings greater than zero!
   */
  StkFloat nextOut( void );

  //! Input one sample to the filter and return one output.
  StkFloat tick( StkFloat input );

  //! Take a channel of the StkFrames object as inputs to the filter and replace with corresponding outputs.
  /*!
    The StkFrames argument reference is returned.  The \c channel
    argument must be less than the number of channels in the
    StkFrames argument (the first channel is specified by 0).
    However, range checking is only performed if _STK_DEBUG_ is
    defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  StkFrames& tick( StkFrames& frames, unsigned int channel = 0 );

  //! Take a channel of the \c iFrames object as inputs to the filter and write outputs to the \c oFrames object.
  /*!
    The \c iFrames object reference is returned.  Each channel
    argument must be less than the number of channels in the
    corresponding StkFrames argument (the first channel is specified
    by 0).  However, range checking is only performed if _STK_DEBUG_
    is defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  StkFrames& tick( StkFrames& iFrames, StkFrames &oFrames, unsigned int iChannel = 0, unsigned int oChannel = 0 );

protected:  

  unsigned long inPoint_;
  unsigned long outPoint_;
  StkFloat delay_;
  StkFloat alpha_;
  StkFloat coeff_;
  StkFloat apInput_;
  StkFloat nextOutput_;
  bool doNextOut_;
};

inline StkFloat DelayA :: nextOut( void )
{
  if ( doNextOut_ ) {
    // Do allpass interpolation delay.
    nextOutput_ = -coeff_ * lastFrame_[0];
    nextOutput_ += apInput_ + ( coeff_ * inputs_[outPoint_] );
    doNextOut_ = false;
  }

  return nextOutput_;
}

inline StkFloat DelayA :: tick( StkFloat input )
{
  inputs_[inPoint_++] = input * gain_;

  // Increment input pointer modulo length.
  if ( inPoint_ == inputs_.size() )
    inPoint_ = 0;

  lastFrame_[0] = nextOut();
  doNextOut_ = true;

  // Save the allpass input and increment modulo length.
  apInput_ = inputs_[outPoint_++];
  if ( outPoint_ == inputs_.size() )
    outPoint_ = 0;

  return lastFrame_[0];
}

inline StkFrames& DelayA :: tick( StkFrames& frames, unsigned int channel )
{
#if defined(_STK_DEBUG_)
  if ( channel >= frames.channels() ) {
    oStream_ << "DelayA::tick(): channel and StkFrames arguments are incompatible!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }
#endif

  StkFloat *samples = &frames[channel];
  unsigned int hop = frames.channels();
  for ( unsigned int i=0; i<frames.frames(); i++, samples += hop ) {
    inputs_[inPoint_++] = *samples * gain_;
    if ( inPoint_ == inputs_.size() ) inPoint_ = 0;
    *samples = nextOut();
    lastFrame_[0] = *samples;
    doNextOut_ = true;
    apInput_ = inputs_[outPoint_++];
    if ( outPoint_ == inputs_.size() ) outPoint_ = 0;
  }

  return frames;
}

inline StkFrames& DelayA :: tick( StkFrames& iFrames, StkFrames& oFrames, unsigned int iChannel, unsigned int oChannel )
{
#if defined(_STK_DEBUG_)
  if ( iChannel >= iFrames.channels() || oChannel >= oFrames.channels() ) {
    oStream_ << "DelayA::tick(): channel and StkFrames arguments are incompatible!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }
#endif

  StkFloat *iSamples = &iFrames[iChannel];
  StkFloat *oSamples = &oFrames[oChannel];
  unsigned int iHop = iFrames.channels(), oHop = oFrames.channels();
  for ( unsigned int i=0; i<iFrames.frames(); i++, iSamples += iHop, oSamples += oHop ) {
    inputs_[inPoint_++] = *iSamples * gain_;
    if ( inPoint_ == inputs_.size() ) inPoint_ = 0;
    *oSamples = nextOut();
    lastFrame_[0] = *oSamples;
    doNextOut_ = true;
    apInput_ = inputs_[outPoint_++];
    if ( outPoint_ == inputs_.size() ) outPoint_ = 0;
  }

  return iFrames;
}


/***************************************************/
/*! \class DelayL
    \brief STK linear interpolating delay line class.
    This class implements a fractional-length digital delay-line using
    first-order linear interpolation.  If the delay and maximum length
    are not specified during instantiation, a fixed maximum length of
    4095 and a delay of zero is set.
    Linear interpolation is an efficient technique for achieving
    fractional delay lengths, though it does introduce high-frequency
    signal attenuation to varying degrees depending on the fractional
    delay setting.  The use of higher order Lagrange interpolators can
    typically improve (minimize) this attenuation characteristic.
    by Perry R. Cook and Gary P. Scavone, 1995--2021.
*/
/***************************************************/

class DelayL : public Filter
{
public:

  //! Default constructor creates a delay-line with maximum length of 4095 samples and zero delay.
  /*!
    An StkError will be thrown if the delay parameter is less than
    zero, the maximum delay parameter is less than one, or the delay
    parameter is greater than the maxDelay value.
   */  
  DelayL( StkFloat delay = 0.0, unsigned long maxDelay = 4095 );

  //! Class destructor.
  ~DelayL();

  //! Get the maximum delay-line length.
  unsigned long getMaximumDelay( void ) { return inputs_.size() - 1; };

  //! Set the maximum delay-line length.
  /*!
    This method should generally only be used during initial setup
    of the delay line.  If it is used between calls to the tick()
    function, without a call to clear(), a signal discontinuity will
    likely occur.  If the current maximum length is greater than the
    new length, no memory allocation change is made.
  */
  void setMaximumDelay( unsigned long delay );

  //! Set the delay-line length.
  /*!
    The valid range for \e delay is from 0 to the maximum delay-line length.
  */
  void setDelay( StkFloat delay );

  //! Return the current delay-line length.
  StkFloat getDelay( void ) const { return delay_; };

  //! Return the value at \e tapDelay samples from the delay-line input.
  /*!
    The tap point is determined modulo the delay-line length and is
    relative to the last input value (i.e., a tapDelay of zero returns
    the last input value).
  */
  StkFloat tapOut( unsigned long tapDelay );

  //! Set the \e value at \e tapDelay samples from the delay-line input.
  void tapIn( StkFloat value, unsigned long tapDelay );

  //! Return the last computed output value.
  StkFloat lastOut( void ) const { return lastFrame_[0]; };

  //! Return the value which will be output by the next call to tick().
  /*!
    This method is valid only for delay settings greater than zero!
   */
  StkFloat nextOut( void );

  //! Input one sample to the filter and return one output.
  StkFloat tick( StkFloat input );

  //! Take a channel of the StkFrames object as inputs to the filter and replace with corresponding outputs.
  /*!
    The StkFrames argument reference is returned.  The \c channel
    argument must be less than the number of channels in the
    StkFrames argument (the first channel is specified by 0).
    However, range checking is only performed if _STK_DEBUG_ is
    defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  StkFrames& tick( StkFrames& frames, unsigned int channel = 0 );

  //! Take a channel of the \c iFrames object as inputs to the filter and write outputs to the \c oFrames object.
  /*!
    The \c iFrames object reference is returned.  Each channel
    argument must be less than the number of channels in the
    corresponding StkFrames argument (the first channel is specified
    by 0).  However, range checking is only performed if _STK_DEBUG_
    is defined during compilation, in which case an out-of-range value
    will trigger an StkError exception.
  */
  StkFrames& tick( StkFrames& iFrames, StkFrames &oFrames, unsigned int iChannel = 0, unsigned int oChannel = 0 );

 protected:

  unsigned long inPoint_;
  unsigned long outPoint_;
  StkFloat delay_;
  StkFloat alpha_;
  StkFloat omAlpha_;
  StkFloat nextOutput_;
  bool doNextOut_;
};

inline StkFloat DelayL :: nextOut( void )
{
  if ( doNextOut_ ) {
    // First 1/2 of interpolation
    nextOutput_ = inputs_[outPoint_] * omAlpha_;
    // Second 1/2 of interpolation
    if (outPoint_+1 < inputs_.size())
      nextOutput_ += inputs_[outPoint_+1] * alpha_;
    else
      nextOutput_ += inputs_[0] * alpha_;
    doNextOut_ = false;
  }

  return nextOutput_;
}

inline void DelayL :: setDelay( StkFloat delay )
{
  if ( delay + 1 > inputs_.size() ) { // The value is too big.
    oStream_ << "DelayL::setDelay: argument (" << delay << ") greater than  maximum!";
    handleError( StkError::WARNING ); return;
  }

  if (delay < 0 ) {
    oStream_ << "DelayL::setDelay: argument (" << delay << ") less than zero!";
    handleError( StkError::WARNING ); return;
  }

  StkFloat outPointer = inPoint_ - delay;  // read chases write
  delay_ = delay;

  while ( outPointer < 0 )
    outPointer += inputs_.size(); // modulo maximum length

  outPoint_ = (long) outPointer;   // integer part

  alpha_ = outPointer - outPoint_; // fractional part
  omAlpha_ = (StkFloat) 1.0 - alpha_;

  if ( outPoint_ == inputs_.size() ) outPoint_ = 0;
  doNextOut_ = true;
}

inline StkFloat DelayL :: tick( StkFloat input )
{
  inputs_[inPoint_++] = input * gain_;

  // Increment input pointer modulo length.
  if ( inPoint_ == inputs_.size() )
    inPoint_ = 0;

  lastFrame_[0] = nextOut();
  doNextOut_ = true;

  // Increment output pointer modulo length.
  if ( ++outPoint_ == inputs_.size() )
    outPoint_ = 0;

  return lastFrame_[0];
}

inline StkFrames& DelayL :: tick( StkFrames& frames, unsigned int channel )
{
#if defined(_STK_DEBUG_)
  if ( channel >= frames.channels() ) {
    oStream_ << "DelayL::tick(): channel and StkFrames arguments are incompatible!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }
#endif

  StkFloat *samples = &frames[channel];
  unsigned int hop = frames.channels();
  for ( unsigned int i=0; i<frames.frames(); i++, samples += hop ) {
    inputs_[inPoint_++] = *samples * gain_;
    if ( inPoint_ == inputs_.size() ) inPoint_ = 0;
    *samples = nextOut();
    doNextOut_ = true;
    if ( ++outPoint_ == inputs_.size() ) outPoint_ = 0;
  }

  lastFrame_[0] = *(samples-hop);
  return frames;
}

inline StkFrames& DelayL :: tick( StkFrames& iFrames, StkFrames& oFrames, unsigned int iChannel, unsigned int oChannel )
{
#if defined(_STK_DEBUG_)
  if ( iChannel >= iFrames.channels() || oChannel >= oFrames.channels() ) {
    oStream_ << "DelayL::tick(): channel and StkFrames arguments are incompatible!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }
#endif

  StkFloat *iSamples = &iFrames[iChannel];
  StkFloat *oSamples = &oFrames[oChannel];
  unsigned int iHop = iFrames.channels(), oHop = oFrames.channels();
  for ( unsigned int i=0; i<iFrames.frames(); i++, iSamples += iHop, oSamples += oHop ) {
    inputs_[inPoint_++] = *iSamples * gain_;
    if ( inPoint_ == inputs_.size() ) inPoint_ = 0;
    *oSamples = nextOut();
    doNextOut_ = true;
    if ( ++outPoint_ == inputs_.size() ) outPoint_ = 0;
  }

  lastFrame_[0] = *(oSamples-oHop);
  return iFrames;
}



Delay :: Delay( unsigned long delay, unsigned long maxDelay )
{
  // Writing before reading allows delays from 0 to length-1. 
  // If we want to allow a delay of maxDelay, we need a
  // delay-line of length = maxDelay+1.
  if ( delay > maxDelay ) {
    oStream_ << "Delay::Delay: maxDelay must be > than delay argument!\n";
    handleError( StkError::FUNCTION_ARGUMENT );
  }

  if ( ( maxDelay + 1 ) > inputs_.size() )
    inputs_.resize( maxDelay + 1, 1, 0.0 );

  inPoint_ = 0;
  this->setDelay( delay );
}

Delay :: ~Delay()
{
}

void Delay :: setMaximumDelay( unsigned long delay )
{
  if ( delay < inputs_.size() ) return;
  inputs_.resize( delay + 1, 1, 0.0 );
}

void Delay :: setDelay( unsigned long delay )
{
  if ( delay > inputs_.size() - 1 ) { // The value is too big.
    oStream_ << "Delay::setDelay: argument (" << delay << ") greater than maximum!\n";
    handleError( StkError::WARNING ); return;
  }

  // read chases write
  if ( inPoint_ >= delay ) outPoint_ = inPoint_ - delay;
  else outPoint_ = inputs_.size() + inPoint_ - delay;
  delay_ = delay;
}

StkFloat Delay :: energy( void ) const
{
  unsigned long i;
  StkFloat e = 0;
  if ( inPoint_ >= outPoint_ ) {
    for ( i=outPoint_; i<inPoint_; i++ ) {
      StkFloat t = inputs_[i];
      e += t*t;
    }
  } else {
    for ( i=outPoint_; i<inputs_.size(); i++ ) {
      StkFloat t = inputs_[i];
      e += t*t;
    }
    for ( i=0; i<inPoint_; i++ ) {
      StkFloat t = inputs_[i];
      e += t*t;
    }
  }
  return e;
}

StkFloat Delay :: tapOut( unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  return inputs_[tap];
}

void Delay :: tapIn( StkFloat value, unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  inputs_[tap] = value;
}

StkFloat Delay :: addTo( StkFloat value, unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  return inputs_[tap]+= value;
}


DelayA :: DelayA( StkFloat delay, unsigned long maxDelay )
{
  if ( delay < 0.5 ) {
    oStream_ << "DelayA::DelayA: delay must be >= 0.5!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }

  if ( delay > (StkFloat) maxDelay ) {
    oStream_ << "DelayA::DelayA: maxDelay must be > than delay argument!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }

  // Writing before reading allows delays from 0 to length-1. 
  if ( maxDelay + 1 > inputs_.size() )
    inputs_.resize( maxDelay + 1, 1, 0.0 );

  inPoint_ = 0;
  this->setDelay( delay );
  apInput_ = 0.0;
  doNextOut_ = true;
}

DelayA :: ~DelayA()
{
}

void DelayA :: clear()
{
  for ( unsigned int i=0; i<inputs_.size(); i++ )
    inputs_[i] = 0.0;
  lastFrame_[0] = 0.0;
  apInput_ = 0.0;
}

void DelayA :: setMaximumDelay( unsigned long delay )
{
  if ( delay < inputs_.size() ) return;
  inputs_.resize(delay + 1, 1, 0.0);
}

void DelayA :: setDelay( StkFloat delay )
{
  unsigned long length = inputs_.size();
  if ( delay + 1 > length ) { // The value is too big.
    oStream_ << "DelayA::setDelay: argument (" << delay << ") greater than maximum!";
    handleError( StkError::WARNING ); return;
  }

  if ( delay < 0.5 ) {
    oStream_ << "DelayA::setDelay: argument (" << delay << ") less than 0.5 not possible!";
    handleError( StkError::WARNING );
  }

  StkFloat outPointer = inPoint_ - delay + 1.0;     // outPoint chases inpoint
  delay_ = delay;

  while ( outPointer < 0 )
    outPointer += length;  // modulo maximum length

  outPoint_ = (long) outPointer;         // integer part
  if ( outPoint_ == length ) outPoint_ = 0;
  alpha_ = 1.0 + outPoint_ - outPointer; // fractional part

  if ( alpha_ < 0.5 ) {
    // The optimal range for alpha is about 0.5 - 1.5 in order to
    // achieve the flattest phase delay response.
    outPoint_ += 1;
    if ( outPoint_ >= length ) outPoint_ -= length;
    alpha_ += (StkFloat) 1.0;
  }

  coeff_ = (1.0 - alpha_) / (1.0 + alpha_);  // coefficient for allpass
}

StkFloat DelayA :: tapOut( unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  return inputs_[tap];
}

void DelayA :: tapIn( StkFloat value, unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  inputs_[tap] = value;
}


DelayL :: DelayL( StkFloat delay, unsigned long maxDelay )
{
  if ( delay < 0.0 ) {
    oStream_ << "DelayL::DelayL: delay must be >= 0.0!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }

  if ( delay > (StkFloat) maxDelay ) {
    oStream_ << "DelayL::DelayL: maxDelay must be > than delay argument!";
    handleError( StkError::FUNCTION_ARGUMENT );
  }

  // Writing before reading allows delays from 0 to length-1. 
  if ( maxDelay + 1 > inputs_.size() )
    inputs_.resize( maxDelay + 1, 1, 0.0 );

  inPoint_ = 0;
  this->setDelay( delay );
  doNextOut_ = true;
}

DelayL :: ~DelayL()
{
}

void DelayL :: setMaximumDelay( unsigned long delay )
{
  if ( delay < inputs_.size() ) return;
  inputs_.resize(delay + 1, 1, 0.0);
}

StkFloat DelayL :: tapOut( unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  return inputs_[tap];
}

void DelayL :: tapIn( StkFloat value, unsigned long tapDelay )
{
  long tap = inPoint_ - tapDelay - 1;
  while ( tap < 0 ) // Check for wraparound.
    tap += inputs_.size();

  inputs_[tap] = value;
}

Echo :: Echo( unsigned long maximumDelay ) : Effect()
{
  this->setMaximumDelay( maximumDelay );
  delayLine_.setDelay( length_ >> 1 );
  effectMix_ = 0.5;
  this->clear();
}

void Echo :: clear( void )
{
  delayLine_.clear();
  lastFrame_[0] = 0.0;
}

void Echo :: setMaximumDelay( unsigned long delay )
{
  if ( delay == 0 ) {
    oStream_ << "Echo::setMaximumDelay: parameter cannot be zero!";
    handleError( StkError::WARNING ); return;
  }

  length_ = delay;
  delayLine_.setMaximumDelay( delay );
}

void Echo :: setDelay( unsigned long delay )
{
  if ( delay > length_ ) {
    oStream_ << "Echo::setDelay: parameter is greater than maximum delay length!";
    handleError( StkError::WARNING ); return;
  }

  delayLine_.setDelay( delay );
}