	

#define V_ENVELOPE_FOLLOWER_NUM_POINTS      2000
class vEnvelopeFollower :
{
      public:
            vEnvelopeFollower();
            virtual ~vEnvelopeFollower();
            inline void Calculate(float *b)
            {
                    envelopeVal -= *buff;
                    if (*b < 0)
                            envelopeVal += *buff = -*b;
                    else
                            envelopeVal += *buff = *b;
                    if (buff++ == bufferEnd)
                            buff = buffer;
            }
            void SetBufferSize(float value);
            void GetControlValue(){return envelopeVal / (float)bufferSize;}

    private:
            float   buffer[V_ENVELOPE_FOLLOWER_NUM_POINTS];
            float   *bufferEnd, *buff, envelopeVal;
            int     bufferSize;
      float val;
};

vEnvelopeFollower::vEnvelopeFollower()
{
    bufferEnd = buffer + V_ENVELOPE_FOLLOWER_NUM_POINTS-1;
    buff = buffer;
    val = 0;
    float *b = buffer;
    do
    {
            *b++ = 0;
    }while (b <= bufferEnd);
    bufferSize = V_ENVELOPE_FOLLOWER_NUM_POINTS;
    envelopeVal= 0;
}

vEnvelopeFollower::~vEnvelopeFollower()
{

}

void vEnvelopeFollower::SetBufferSize(float value)
{

    bufferEnd = buffer + (bufferSize = 100 + (int)(value * ((float)V_ENVELOPE_FOLLOWER_NUM_POINTS-102)));
    buff = buffer;
    float val =  envelopeVal / bufferSize;
    do
    {
            *buff++ = val;
    }while (buff <= bufferEnd);
    buff = buffer;
}
