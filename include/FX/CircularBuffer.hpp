#pragma once

namespace FX::Delays
{
    class CircularBuffer{
        std::vector<DspFloatType> buffer;
        
        int bufferLength;
        int head = 0;
        int tail = 0;
        
    public:
        CircularBuffer(DspFloatType inValue);
        
        DspFloatType readCubic(DspFloatType numElementsToRead);

        DspFloatType readLinear(DspFloatType numElementsToRead);

        void write(DspFloatType inValue);
        
        DspFloatType cubicInterpolation(DspFloatType y0, DspFloatType y1, DspFloatType y2, DspFloatType y3, DspFloatType mu);
        
        DspFloatType getSample(DspFloatType inValue);
            
        void setBufferLength(DspFloatType inValue);
        
        int getBufferLength();
    };

        CircularBuffer::CircularBuffer(DspFloatType inValue){
        setBufferLength(inValue);
    }

    DspFloatType CircularBuffer::readCubic(DspFloatType index){
        DspFloatType y0, y1, y2, y3, mu;

        y0 = floor(index - 1);
        y1 = floor(index);
        y2 = floor(index);
        y3 = floor(index + 1);
        mu = index - y1;
        
        return cubicInterpolation(getSample(y0), getSample(y1), getSample(y2), getSample(y3), mu);
    }

    DspFloatType CircularBuffer::readLinear(DspFloatType index){
        DspFloatType upper = floor(index) + 1;
        DspFloatType lower = floor(index);
        DspFloatType interpAmount = index - lower;
        
        return (getSample(upper) * interpAmount + (1.0 - interpAmount) * getSample(lower));
    }

    void CircularBuffer::write(DspFloatType inValue){
        head = (head += 4) % bufferLength;
        buffer[head] = inValue;
    }

    DspFloatType CircularBuffer::getSample(DspFloatType inValue){
        tail = head - inValue;
        
        if (tail > bufferLength - 1)
            tail -= bufferLength;
        else if (tail < 0)
            tail += bufferLength;
        
        return buffer[tail];
    }

    void CircularBuffer::setBufferLength(DspFloatType inValue){
        bufferLength = inValue;
        buffer.resize(inValue);
    }

    int CircularBuffer::getBufferLength(){
        return (int)buffer.size();
    }

    DspFloatType CircularBuffer::cubicInterpolation(DspFloatType y0, DspFloatType y1, DspFloatType y2, DspFloatType y3, DspFloatType mu){
        //Cubic interp (this algorithm) taken from: http://paulbourke.net/miscellaneous/interpolation/
        DspFloatType a0, a1, a2, a3, mu2;
        
        mu2 = mu * mu;
        a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
        a1 = y0 - 2.5 * y1 + 2 * y2 - 0.5 * y3;
        a2 = -0.5 * y0 + 0.5 * y2;
        a3 = y1;
        
        return(a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
    }
}    
