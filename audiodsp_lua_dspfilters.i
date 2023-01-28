%module DspFilters
%{
#include <iostream>
#include <vector>
#include <Undenormal.hpp>
#include "DspFilters/Dsp.h"
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_math.i"

namespace Dsp {

    typedef std::complex<double> complex_t;
    typedef std::pair<complex_t, complex_t> complex_pair_t;

    // A conjugate or real pair
    struct ComplexPair : complex_pair_t
    {
        ComplexPair ();
        explicit ComplexPair (const complex_t& c1);
        ComplexPair (const complex_t& c1,
                    const complex_t& c2);
        bool isConjugate () const;
        bool isReal () const;
        bool isMatchedPair () const;
        bool is_nan () const;
    };

    struct PoleZeroPair
    {
        ComplexPair poles;
        ComplexPair zeros;

        PoleZeroPair ();

        // single pole/zero
        PoleZeroPair (const complex_t& p, const complex_t& z);
        PoleZeroPair (const complex_t& p1, const complex_t& z1,
                        const complex_t& p2, const complex_t& z2);
        bool isSinglePole () const;
        bool is_nan () const;
    };

    // Identifies the general class of filter
    enum Kind
    {
        kindLowPass,
        kindHighPass,
        kindBandPass,
        kindBandStop,
        kindLowShelf,
        kindHighShelf,
        kindBandShelf,
        kindOther
    };
}


namespace Dsp {
    // Unique IDs to help identify parameters
    enum ParamID
    {
    idSampleRate,
    idFrequency,
    idQ,
    idBandwidth,
    idBandwidthHz,
    idGain,
    idSlope,
    idOrder,
    idRippleDb,
    idStopDb,
    idRolloff,

    idPoleRho,
    idPoleTheta,
    idZeroRho,
    idZeroTheta,

    idPoleReal,
    idZeroReal
    };

    enum
    {
    maxParameters = 8
    };

    struct Params
    {
    void clear ()
    {
        for (int i = 0; i < maxParameters; ++i)
        value[i] = 0;
    }

    double& operator[] (int index)
    {
        return value[index];
    }

    const double& operator[] (int index) const
    {
        return value[index];
    }

    %extend {
        double __getitem__(size_t i) { return (*$self)[i]; }
        void   __setitem__(size_t i, double v) { (*$self)[i] = v; }
    }
    double value[maxParameters];
    };

    //
    // Provides meta-information about a filter parameter
    // to achieve run-time introspection.
    //
    class ParamInfo
    {
    public:
        typedef double (ParamInfo::*toControlValue_t) (double) const;
        typedef double (ParamInfo::*toNativeValue_t) (double) const;
        typedef std::string (ParamInfo::*toString_t) (double) const;

        // dont use this one
        ParamInfo (); // throws std::logic_error

        ParamInfo (ParamID id,
                    const char* szLabel,
                    const char* szName,
                    double arg1,
                    double arg2,
                    double defaultNativeValue,
                    toControlValue_t toControlValue_proc,
                    toNativeValue_t toNativeValue_proc,
                    toString_t toString_proc);

        ParamID getId () const;
        const char* getLabel () const;
        const char* getName () const;
        double getDefaultValue () const;
        double toControlValue (double nativeValue) const;
        double toNativeValue (double controlValue) const;
        std::string toString (double nativeValue) const;
        double clamp (double nativeValue) const;
        double Int_toControlValue (double nativeValue) const;
        double Int_toNativeValue (double controlValue) const;

        double Real_toControlValue (double nativeValue) const;
        double Real_toNativeValue (double controlValue) const;

        double Log_toControlValue (double nativeValue) const;
        double Log_toNativeValue (double controlValue) const;

        double Pow2_toControlValue (double nativeValue) const;
        double Pow2_toNativeValue (double controlValue) const;

        std::string Int_toString (double nativeValue) const;
        std::string Hz_toString (double nativeValue) const;
        std::string Real_toString (double nativeValue) const;
        std::string Db_toString (double nativeValue) const;

        //
        // Creates the specified ParamInfo
        //

        static ParamInfo defaultSampleRateParam ();
        static ParamInfo defaultCutoffFrequencyParam ();
        static ParamInfo defaultCenterFrequencyParam ();
        static ParamInfo defaultQParam ();
        static ParamInfo defaultBandwidthParam ();
        static ParamInfo defaultBandwidthHzParam ();
        static ParamInfo defaultGainParam ();
        static ParamInfo defaultSlopeParam ();
        static ParamInfo defaultRippleDbParam ();
        static ParamInfo defaultStopDbParam ();
        static ParamInfo defaultRolloffParam ();
        static ParamInfo defaultPoleRhoParam ();
        static ParamInfo defaultPoleThetaParam ();
        static ParamInfo defaultZeroRhoParam ();
        static ParamInfo defaultZeroThetaParam ();
        static ParamInfo defaultPoleRealParam ();
        static ParamInfo defaultZeroRealParam ();
    };

    class Filter
    {
    public:
        virtual ~Filter();

        virtual Kind getKind () const = 0;

        virtual const std::string getName () const = 0;

        virtual int getNumParams () const = 0;  

        virtual ParamInfo getParamInfo (int index) const = 0;

        Params getDefaultParams() const;

        const Params& getParams() const;
        double getParam (int paramIndex) const;
        void setParam (int paramIndex, double nativeValue);
        int findParamId (int paramId);

        void setParamById (int paramId, double nativeValue);

        void setParams (const Params& parameters);
        void copyParamsFrom (Dsp::Filter const* other);

        virtual std::vector<PoleZeroPair> getPoleZeros() const = 0;        
        virtual complex_t response (double normalizedFrequency) const = 0;

        virtual int getNumChannels() = 0;
        virtual void reset () = 0;
        virtual void process (int numSamples, float* const* arrayOfChannels) = 0;
        virtual void process (int numSamples, double* const* arrayOfChannels) = 0;
    };
}

namespace Dsp {

    template <class DesignClass,
            int Channels,
            class StateType = DirectFormII>
    class SmoothedFilterDesign
    : public FilterDesign <DesignClass,
                        Channels,
                        StateType>
    {
    public:
        typedef FilterDesign <DesignClass, Channels, StateType> filter_type_t;

        SmoothedFilterDesign (int transitionSamples);

        template <typename Sample>
        void processBlock (int numSamples,
                            Sample* const* destChannelArray);
        void process (int numSamples, float* const* arrayOfChannels);
        void process (int numSamples, double* const* arrayOfChannels);

        %extend {
            void ProcessBlock(size_t n, float * in, float * out) {
                float *temp[2];
                memcpy(out,in,n*sizeof(float));
                temp[0] = out;
                temp[1] = out;            
                $self->process(n,temp);                                    
            }
            void ProcessBlock(size_t n, float ** in, float ** out) {
                float *temp[2];
                memcpy(out[0],in[0],n*sizeof(float));
                memcpy(out[1],in[1],n*sizeof(float));
                temp[0] = out[0];
                temp[1] = out[1];            
                $self->process(n,temp);                                
            }
            void ProcessBlock(size_t n, double * in, double * out) {
                double *temp[2];
                memcpy(out,in,n*sizeof(double));
                temp[0] = out;
                temp[1] = out;            
                $self->process(n,temp);                                    
            }
            void ProcessBlock(size_t n, double ** in, double ** out) {
                double *temp[2];
                memcpy(out[0],in[0],n*sizeof(double));
                memcpy(out[1],in[1],n*sizeof(double));
                temp[0] = out[0];
                temp[1] = out[1];            
                $self->process(n,temp);                                
            }
        }
        
        //Params getDefaultParams() const;

        const Params& getParams() const;
        double getParam (int paramIndex) const;
        void setParam (int paramIndex, double nativeValue);
        int findParamId (int paramId);

        void setParamById (int paramId, double nativeValue);

        void setParams (const Params& parameters);
        void copyParamsFrom (Dsp::Filter const* other);    
    };
}

%template (float_vector_stereo) std::vector<float*>;
%template (double_vector_stereo) std::vector<double*>;

%template (RBJLowPass) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::LowPass,1>;
%template (StereoRBJLowPass) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::LowPass,2>;

%template (RBJHighPass) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::HighPass,1>;
%template (StereoRBJHighPass) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::HighPass,2>;

%template (RBJBandPass1) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandPass1,1>;
%template (StereoRBJBandPass1) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandPass1,2>;

%template (RBJBandPass2) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandPass2,1>;
%template (StereoRBJBandPass2) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandPass2,2>;

%template (RBJBandStop) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandStop,1>;
%template (StereoRBJBandStop) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandStop,2>;

%template (RBJLowShelf) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::LowShelf,1>;
%template (StereoRBJLowShelf) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::LowShelf,2>;

%template (RBJHighShelf) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::HighShelf,1>;
%template (StereoRBJHighShelf) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::HighShelf,2>;

%template (RBJBandShelf) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandShelf,1>;
%template (StereoRBJBandShelf) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::BandShelf,2>;

%template (RBJAllPass) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::AllPass,1>;
%template (StereoRBJAllPass) Dsp::SmoothedFilterDesign<Dsp::RBJ::Design::AllPass,2>;

%template (ButterworthBandPass) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::BandPass<64>,1>;
%template (StereoButterworthBandPass) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::BandPass<64>,2>;

%template (ButterworthBandStop) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::BandStop<64>,1>;
%template (StereoButterworthBandStop) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::BandStop<64>,2>;

%template (ButterworthLowPass) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::LowPass<64>,1>;
%template (StereoButterworthLowPass) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::LowPass<64>,2>;

%template (ButterworthHighPass) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::HighPass<64>,1>;
%template (StereoButterworthHighPass) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::HighPass<64>,2>;

%template (ButterworthBandShelf) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::BandShelf<64>,1>;
%template (StereoButterworthBandShelf) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::BandShelf<64>,2>;

%template (ButterworthLowShelf) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::LowShelf<64>,1>;
%template (StereoButterworthLowShelf) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::LowShelf<64>,2>;

%template (ButterworthHighShelf) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::HighShelf<64>,1>;
%template (StereoButterworthHighSHelf) Dsp::SmoothedFilterDesign<Dsp::Butterworth::Design::HighShelf<64>,2>;



%template (ChebyshevIBandPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::BandPass<64>,1>;
%template (StereoChebyshevIBandPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::BandPass<64>,2>;

%template (ChebyshevIBandStop) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::BandStop<64>,1>;
%template (StereoChebyshevIBandStop) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::BandStop<64>,2>;

%template (ChebyshevILowPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::LowPass<64>,1>;
%template (StereoChebyshevILowPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::LowPass<64>,2>;

%template (ChebyshevIHighPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::HighPass<64>,1>;
%template (StereoChebyshevIHighPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::HighPass<64>,2>;

%template (ChebyshevIBandShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::BandShelf<64>,1>;
%template (StereoChebyshevIBandShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::BandShelf<64>,2>;

%template (ChebyshevILowShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::LowShelf<64>,1>;
%template (StereoChebyshevILowShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::LowShelf<64>,2>;

%template (ChebyshevIHighShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::HighShelf<64>,1>;
%template (StereoChebyshevIHighSHelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevI::Design::HighShelf<64>,2>;



%template (ChebyshevIIBandPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::BandPass<64>,1>;
%template (StereoChebyshevIIBandPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::BandPass<64>,2>;

%template (ChebyshevIIBandStop) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::BandStop<64>,1>;
%template (StereoChebyshevIIBandStop) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::BandStop<64>,2>;

%template (ChebyshevIILowPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::LowPass<64>,1>;
%template (StereoChebyshevIILowPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::LowPass<64>,2>;

%template (ChebyshevIIHighPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::HighPass<64>,1>;
%template (StereoChebyshevIIHighPass) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::HighPass<64>,2>;

%template (ChebyshevIIBandShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::BandShelf<64>,1>;
%template (StereoChebyshevIIBandShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::BandShelf<64>,2>;

%template (ChebyshevIILowShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::LowShelf<64>,1>;
%template (StereoChebyshevIILowShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::LowShelf<64>,2>;

%template (ChebyshevIIHighShelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::HighShelf<64>,1>;
%template (StereoChebyshevIIHighSHelf) Dsp::SmoothedFilterDesign<Dsp::ChebyshevII::Design::HighShelf<64>,2>;

%template (EllipticBandPass) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::BandPass<64>,1>;
%template (StereoEllipticBandPass) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::BandPass<64>,2>;

%template (EllipticBandStop) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::BandStop<64>,1>;
%template (StereoEllipticBandStop) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::BandStop<64>,2>;

%template (EllipticLowPass) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::LowPass<64>,1>;
%template (StereoEllipticLowPass) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::LowPass<64>,2>;

%template (EllipticHighPass) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::HighPass<64>,1>;
%template (StereoEllipticHighPass) Dsp::SmoothedFilterDesign<Dsp::Elliptic::Design::HighPass<64>,2>;


%template (BesselBandPass) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::BandPass<64>,1>;
%template (StereoBesselBandPass) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::BandPass<64>,2>;

%template (BesselBandStop) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::BandStop<64>,1>;
%template (StereoBesselBandStop) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::BandStop<64>,2>;

%template (BesselLowPass) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::LowPass<64>,1>;
%template (StereoBesselLowPass) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::LowPass<64>,2>;

%template (BesselHighPass) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::HighPass<64>,1>;
%template (StereoBesselHighPass) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::HighPass<64>,2>;

%template (BesselLowShelf) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::LowShelf<64>,1>;
%template (StereoBesselLowShelf) Dsp::SmoothedFilterDesign<Dsp::Bessel::Design::LowShelf<64>,2>;

%template (LegendreBandPass) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::BandPass<64>,1>;
%template (StereoLegendreBandPass) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::BandPass<64>,2>;

%template (LegendreBandStop) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::BandStop<64>,1>;
%template (StereoLegendreBandStop) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::BandStop<64>,2>;

%template (LegendreLowPass) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::LowPass<64>,1>;
%template (StereoLegendreLowPass) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::LowPass<64>,2>;

%template (LegendreHighPass) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::HighPass<64>,1>;
%template (StereoLegendreHighPass) Dsp::SmoothedFilterDesign<Dsp::Legendre::Design::HighPass<64>,2>;

%template (pole_zero_vector) std::vector<Dsp::PoleZeroPair>;


%template (complex_t)       std::complex<double>;
// need to wrap std_pair.i sometime
//%template (complex_pair_t)  std::pair<std::complex<double>,std::complex<double>>;    