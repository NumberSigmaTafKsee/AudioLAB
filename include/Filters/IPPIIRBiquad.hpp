#pragma once

struct IPPIIRBiquad: public Casino::IPP::IIRBiquad<DspFloatType>
{
    IPPIIRBiquad() = default;
    IPPIIRBiquad(const BiquadSection &c) {
        setCoefficients(c);
    }
    IPPIIRBiquad(const BiquadSOS & sos) {
        setCoefficients(sos);
    }
    void setCoefficients(const BiquadSection & c)
    {
        DspFloatType buf[6];        
        buf[0] = c.z[0];
        buf[1] = c.z[1];
        buf[2] = c.z[2];
        buf[3] = 1.0;
        buf[4] = c.p[0];
        buf[5] = c.p[1];
        this->initCoefficients(blockSize,1,buf);
    }
    void setCoefficients(const BiquadSOS & sos)
    {        
        DspFloatType buf[6*sos.size()];
        int x = 0;
        for(size_t i = 0; i < sos.size(); i++)
        {    
            buf[x++] = sos[i].z[0];
            buf[x++] = sos[i].z[1];
            buf[x++] = sos[i].z[2];
            buf[x++] = 1.0;
            buf[x++] = sos[i].p[0];
            buf[x++] = sos[i].p[1];
        }
        this->initCoefficients(blockSize,sos.size(),buf);
    }
    void setCoefficients(const Filters::FilterCoefficients & c)
    {
        DspFloatType buf[6];        
        buf[0] = c.b[0];
        buf[1] = c.b[1];
        buf[2] = c.b[2];
        buf[3] = 1.0;
        buf[4] = c.a[0];
        buf[5] = c.a[1];
        this->initCoefficients(blockSize,1,buf);
    }    
    void setCoefficients(const std::vector<Filters::FilterCoefficients> & c)
    {
        DspFloatType buf[6*c.size()];
        int x = 0;
        for(size_t i = 0; i < c.size(); i++)
        {    
            buf[x++] = c[i].b[0];
            buf[x++] = c[i].b[1];
            buf[x++] = c[i].b[2];
            buf[x++] = 1.0;
            buf[x++] = c[i].a[0];
            buf[x++] = c[i].a[1];
        }
        this->initCoefficients(blockSize,c.size(),buf);
    }    
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
    {
        assert(n == this->len);
        this->Execute(in,out);
    }
};
