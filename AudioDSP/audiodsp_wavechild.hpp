#pragma once

namespace WDF {
//==============================================================================
// ** 1-PORT ** (base class for every WDF classes)
//==============================================================================
template <typename T>
class OnePort
{
    public:
        OnePort (T R, String n = String::empty)
            : name (n), Rp (R), a (0), b (0), port(this) {}
        //----------------------------------------------------------------------
        virtual String name () const { return name.isEmpty() ? label() : _name; }
        virtual String label () const = 0;
        //----------------------------------------------------------------------
        virtual inline void incident (T wave) = 0;
        //----------------------------------------------------------------------
        virtual inline T reflected () = 0;
        //----------------------------------------------------------------------
        virtual T R () { return port->Rp; }       // Port resistance
        virtual T G () { return 1.0 / port->Rp; } // Port conductance (inv.Rp)
        //----------------------------------------------------------------------
        virtual void connect (OnePort<T>* other)
        {
            port = other;
            other->port = this;
        }
        //----------------------------------------------------------------------
        T voltage () // v
        {
            return (port->a + port->b) / 2.0;
        }
        //----------------------------------------------------------------------
        T current () // i
        {
            return (port->a - port->b) / (port->Rp + port->Rp);
        }
        //----------------------------------------------------------------------
    protected:
        //----------------------------------------------------------------------
        String _name; // Port name
        T Rp; // Port resistance
        //----------------------------------------------------------------------
        T a; // incident wave (incoming wave)
        T b; // reflected wave (outgoing wave)
        //----------------------------------------------------------------------
        OnePort<T>* port; // internal pointer (used for direct connect form)
        //----------------------------------------------------------------------
};
//==============================================================================
// ** 2-PORT **
//==============================================================================
template <typename T>
class TwoPort : public OnePort<T> // parent
{
    public:
        OnePort<T>* child;
        //----------------------------------------------------------------------
        TwoPort (String name = String::empty)
            : OnePort (1.0, name) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "2P"; }
        //----------------------------------------------------------------------
        virtual void connectParent (OnePort<T>* parent)
        {
            OnePort::connect (parent);
        }
        //----------------------------------------------------------------------
        virtual void connectChild (OnePort<T>* port) = 0;
        //----------------------------------------------------------------------
        virtual void connect (OnePort<T>* p, OnePort<T>* c)
        {
            OnePort::connect (p);
            connectChild (c);
        }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            child.port->a = child.port->reflected ();
            computeParentB ();
            return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave;
            computeChildB ();
            child.port->incident (child.port->b);
        }
        //----------------------------------------------------------------------
        virtual inline void computeChildB () = 0;
        virtual inline void computeParentB () = 0;
        //----------------------------------------------------------------------
};
//==============================================================================
// ** 3-PORT **
//==============================================================================
template <typename T>
class ThreePort : public OnePort<T> // adapted
{
    public:
        OnePort<T> *left, *right;
        //----------------------------------------------------------------------
        ThreePort (String name = String::empty)
            : left (nullptr), right (nullptr), OnePort (1.0, name) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "3P"; }
        //----------------------------------------------------------------------
        virtual void connect (OnePort<T>* left, OnePort<T>* right) = 0;
        //----------------------------------------------------------------------
        virtual inline T reflected () = 0;
        virtual inline void incident (T wave) = 0;
        //----------------------------------------------------------------------
};
//==============================================================================
// ** SERIE **
//==============================================================================
template <typename T>
class Serie : public ThreePort<T>
{
    public:
        Serie (String name = "--")
            : ThreePort (name)
        {}
        //----------------------------------------------------------------------
        virtual String label () const { return "--"; }
        //----------------------------------------------------------------------
        virtual void connect (OnePort<T>* l, OnePort<T>* r)
        {
            left = l; right = r;
            port->Rp = (left->R() + right->R());
        }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = -(left->port->reflected()
                     + right->port->reflected());
            return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            register T lrW = (wave + left->port->b + right->port->b);
             left->port->incident ( left->port->b - ( left->R()/port->R()) * lrW);
            right->port->incident (right->port->b - (right->R()/port->R()) * lrW);
            port->a = wave;
        }
        //----------------------------------------------------------------------
};
//==============================================================================
// ** PARALLEL **
//==============================================================================
template <typename T>
class Parallel : public ThreePort<T>
{
    public:
        Parallel (String name = "||")
            : ThreePort (name)
        {}
        //----------------------------------------------------------------------
        virtual String label () const { return "||"; }
        //----------------------------------------------------------------------
        virtual void connect (OnePort<T>* l, OnePort<T>* r)
        {
            left = l; right = r;
            port->Rp = (left->R() * right->R())
                     / (left->R() + right->R());
        }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            register T lrG = left->G() + right->G();
            port->b = ( left->G()/lrG) *  left->port->reflected() +
                      (right->G()/lrG) * right->port->reflected();
            return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            register T lrG = left->G() + right->G();
            register T lrW = (wave + left->port->b + right->port->b);
             left->port->incident ( left->port->b - ( left->G()/lrG) * lrW);
            right->port->incident (right->port->b - (right->G()/lrG) * lrW);
            port->a = wave;
        }
        //----------------------------------------------------------------------
};
//==============================================================================
// ** RESISTOR **
//==============================================================================
template <typename T>
class Resistor : public OnePort<T>
{
    public:
        Resistor (T R, String name = String::empty)
            : OnePort (R, name) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "R"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = 0; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave;
        }
        //----------------------------------------------------------------------
};
//==============================================================================
// ** CAPACITOR **
//==============================================================================
template <typename T>
class Capacitor : public OnePort<T>
{
    public:
        Capacitor (T C, T Fs, String name = String::empty)
            : OnePort (Fs/2.0*C, name), state (0) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "C"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = state; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave; state = port->a;
        }
        //----------------------------------------------------------------------
    private:
        T state;
        //----------------------------------------------------------------------
};
//==============================================================================
// ** INDUCTOR **
//==============================================================================
template <typename T>
class Inductor : public OnePort<T>
{
    public:
        Inductor (T L, T Fs, String name = String::empty)
            : OnePort (2.0*L/Fs, name), state (0) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "L"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = -state; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave; state = port->a;
        }
        //----------------------------------------------------------------------
    private:
        T state;
        //----------------------------------------------------------------------
};
//==============================================================================
// ** OPEN CIRCUIT **
//==============================================================================
template <typename T>
class OpenCircuit : public OnePort<T>
{
    public:
        OpenCircuit (T R, String name = String::empty)
            : OnePort (R, name) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "Oc"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = port->a; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave;
        }
        //----------------------------------------------------------------------
};
//==============================================================================
// ** SHORT CIRCUIT **
//==============================================================================
template <typename T>
class ShortCircuit : public OnePort<T>
{
    public:
        ShortCircuit (T R, String name = String::empty)
            : OnePort (R, name) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "Sc"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = -port->a; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave;
        }
        //----------------------------------------------------------------------
};
//==============================================================================
// ** VOLTAGE SOURCE **
//==============================================================================
template <typename T>
class VoltageSource : public OnePort<T>
{
    public:
        VoltageSource (T V, T R, String name = String::empty)
            : OnePort (R, name), Vs (V) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "Vs"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = -port->a + 2.0 * Vs; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave;
        }
        //----------------------------------------------------------------------
    private:
        T Vs;
        //----------------------------------------------------------------------
};
//==============================================================================
// ** CURRENT SOURCE **
//==============================================================================
template <typename T>
class CurrentSource : public OnePort<T>
{
    public:
        CurrentSource (T I, T R, String name = String::empty)
            : OnePort (R, name), Is (I) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "Is"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            port->b = port->a + 2.0 * R() * Is; return port->b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            port->a = wave;
        }
        //----------------------------------------------------------------------
    private:
        T Is;
        //----------------------------------------------------------------------
};
//==============================================================================
// ** IDEAL TRANSFORMER **
//==============================================================================
template <typename T>
class IdealTransformer : public TwoPort<T>
{
    public:
        IdealTransformer (T ratio, String name = String::empty)
            : TwoPort (1.0, name), N (ratio) {}
        //----------------------------------------------------------------------
        virtual String label () const { return "][" }
        //----------------------------------------------------------------------
        virtual void connectChild (OnePort<T>* port)
        {
            T Rs = port.R();
            port->Rp = Rs / (n*n);
            child.Rp = Rs;
            TwoPort::connectChild (port);
        }
        //----------------------------------------------------------------------
        virtual inline void computeChildB ()
        {
            child.port->b = port->a * (1.0/N);
        }
        //----------------------------------------------------------------------
        virtual inline void computeParentB ()
        {
            port->b = child.port->a * N;
        }
        //----------------------------------------------------------------------
    private:
        T N;
        //----------------------------------------------------------------------
};
//==============================================================================
/**
    EXTRA TEMPLATES
    ---------------
    Not direct WDF classes but interfaces for your non-linear blackboxes.
**/
//==============================================================================
// ** Newton/Raphson ** (implicit equation solver)
//==============================================================================
template <typename T>
class NewtonRaphson
{
    public:
        NewtonRaphson (T guess = 100.0) : xguess (guess) {}
        //----------------------------------------------------------------------
        inline T solve (int max_iter = 100, T epsilon = 1e-9)
        {
            T x = xguess;
            T err = 1e6;
            int iteration = 0;
            while (fabs(err) / fabs(x) > epsilon)
            {
                xguess = iterate (x);
                err = x - xguess;
                x = xguess;
                if (iteration > max_iter) break;
                ++iteration;
            }
            return x;
        }
        //----------------------------------------------------------------------
        inline T iterate (T x, T dx = 1e-6)
        {
            T F = evaluate (x);
            T xNew = x - dx*F / (evaluate (x + dx) - F);
            return xNew;
        }
        //----------------------------------------------------------------------
        virtual inline T evaluate (T x) = 0; // declare your implicit equation
        //----------------------------------------------------------------------
    private :
        T xguess;
        //----------------------------------------------------------------------
};
//==============================================================================
} // namespace WDF
//==============================================================================


namespace Wavechild670 {
//==============================================================================
//----------------------------------------------------------------------
// Level Time Constant 6-way switch parameters (time from 10dB limiting)
//----------------------------------------------------------------------
static const double ltc[6][6] =
{
    //------------------------------------------------------------------
    // CT    CU        CV        RT        RU       RV   | Release Time
    //------------------------------------------------------------------
    { 2e-6, 8e-6, 20e-6,  51.9e3,  10e9,  10e9 }, // 0.3s
    { 2e-6, 8e-6, 20e-6, 149.9e3,  10e9,  10e9 }, // 0.8s
    { 4e-6, 8e-6, 20e-6,   220e3,  10e9,  10e9 }, // 2.0s
    { 8e-6, 8e-6, 20e-6,   220e3,  10e9,  10e9 }, // 5.0s
    { 4e-6, 8e-6, 20e-6,   220e3, 100e3,  10e9 }, // 2.0s / 10.0s
    { 2e-6, 8e-6, 20e-6,   220e3, 100e3, 100e3 }  // 0.3s / 5.0s / 25.0s
    //------------------------------------------------------------------
};
//==============================================================================
template <typename T>
class LevelTimeConstant
{
    public:
        LevelTimeConstant (T Fs)
            : //----------------------------------------------------------------
              R1 (220e3, "RT"),
              R2 (  1e9, "RU"),
              R3 (  1e9, "RV"),
              //----------------------------------------------------------------
              C1 ( 2e-6, Fs, "CT"),
              C2 ( 8e-6, Fs, "CU"),
              C3 (20e-6, Fs, "CV"),
              //----------------------------------------------------------------
        {
            wiring ();
        }
        //----------------------------------------------------------------------
        void parameters (T Fs, const int index)
        {
            jassert(index >= 0 && index < 6);
            update (Fs, ltc[index][0], ltc[index][1], ltc[index][2],
                        ltc[index][3], ltc[index][4], ltc[index][5]);
        }
        //----------------------------------------------------------------------
        T process (T Iin) // Iin == current (current law apply)
        {
            root.incident (root.reflected() - (2.0*(Iin * root.R())));
            return C1.voltage();
        }
        //----------------------------------------------------------------------
    protected:
        WDF::Resistor<T>    R1, R2, R3;
        WDF::Capacitor<T>   C1, C2, C3;
        //----------------------------------------------------------------------
        WDF::Serie<T>       serie_A;
        WDF::Serie<T>       serie_B;
        WDF::Parallel<T>    paral_A;
        WDF::Parallel<T>    paral_B;
        WDF::Parallel<T>    root;
        //----------------------------------------------------------------------
        /**
                --------------------------
                |       |    |     |     |
                |       |    |     R2    R3
              root      R1   C1    |     |
                |       |    |     C2    C3
                |       |    |     |     |
                --------------------------
        **/
        //----------------------------------------------------------------------
        inline void wiring ()
        {
            paral_A.connect (&R1,       &C1);
            serie_A.connect (&R2,       &C2);
            serie_B.connect (&R3,       &C3);
            paral_B.connect (&serie_A,  &serie_B);
               root.connect (&paral_A,  &paral_B);
        }
        //----------------------------------------------------------------------
        void update (T Fs, T CT = 2e-6,  T CU = 8e-6, T CV = 20e-6,
                           T RT = 220e3, T RU = 1e9,  T RV = 1e9)
        {
            T hFs = FS*.5;
            //------------------------------------------------------------------
            C1.Rp = hFs*CT;
            C2.Rp = hFs*CU;
            C3.Rp = hFs*CV;
            //------------------------------------------------------------------
            R1.Rp = RT;
            R2.Rp = RU;
            R3.Rp = RV;
            //------------------------------------------------------------------
            wiring ();
        }
        //----------------------------------------------------------------------
};

//==============================================================================
template <typename T>
class NonIdealTransformer : public WDF::TwoPort<T>
{
    public:
        NonIdealTransformer (T Fs,
                             T Nt,
                             T Lp_, T Rp_,
                             T Lm_, T Rc_,
                             T Ls_, T Rs_, T Cw_,
                             String name = String::empty)
            : OnePort (input.R, name),
              //----------------------------------------------------------------
              // Components
              //----------------------------------------------------------------
              Lp (Lp_, Fs, "Lp"),
              Lm (Lm_, Fs, "Lm"),
              Ls (Ls_, Fs, "Ls"),
              Cw (Cw_, Fs, "Cw"),
              Rp (Rp_,     "Rp"),
              Rc (Rc_,     "Rc"),
              Rs (Rs_,     "Rs"),
              //----------------------------------------------------------------
              transfo (Nt, "][")
              //----------------------------------------------------------------
        {}
        //----------------------------------------------------------------------
        virtual String label () const { "][" };
        //----------------------------------------------------------------------
        virtual void connect (OnePort<T>* parent, OnePort<T>* child)
        {
            root.connect (parent); wiring (child);
        }
        //----------------------------------------------------------------------
        virtual void connectParent (OnePort<T>* parent)
        {
            root.connect (parent);
        }
        //----------------------------------------------------------------------
        virtual void connectChild (OnePort<T>* child)
        {
            wiring (child);
        }
        //----------------------------------------------------------------------
        virtual inline void computeChildB ()
        {
        }
        //----------------------------------------------------------------------
        virtual inline void computeParentB ()
        {
        }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            return root.reflected ();
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            root.incident (wave);
        }
        //----------------------------------------------------------------------
        virtual T R () { return root.R (); }
        //----------------------------------------------------------------------
        inline T Vout () { return Cw.voltage(); }
        //----------------------------------------------------------------------
    protected:
        //----------------------------------------------------------------------
        WDF::IdealTransformer<T>    transfo;
        //----------------------------------------------------------------------
        WDF::Inductor<T>            Lp, Lm, Ls;
        WDF::Resistor<T>            Rp, Rc, Rs;
        WDF::Capacitor<T>           Cw;
        //----------------------------------------------------------------------
        WDF::Serie<T>               root;
        //----------------------------------------------------------------------
        WDF::Serie<T>               serie_A;
        WDF::Serie<T>               serie_B;
        WDF::Serie<T>               serie_C;
        WDF::Parallel<T>            paral_A;
        WDF::Parallel<T>            paral_B;
        WDF::Parallel<T>            paral_C;
        //----------------------------------------------------------------------
        inline void wiring (OnePort<T>* child)
        {
            serie_A.connect (&Lp,       &Rp);
            serie_B.connect (&Rs,       &Ls);
            paral_A.connect (&Lm,       &Rc);
            //------------------------------------------------------------------
            paral_B.connect (child,     &Cw);
            serie_C.connect (&paral_B,  &serie_B);
            //------------------------------------------------------------------
            transfo.connectChild (&serie_C);
            //------------------------------------------------------------------
            paral_C.connect (&transfo,  &paral_A);
               root.connect (&serie_A,  &paral_C);
        }
        //----------------------------------------------------------------------
};
//==============================================================================
template <typename T>
class InputCoupledTransformer : public WDF::OnePort<T>
{
    public:
        InputCoupledTransformer ()
            : OnePort (1.0)
        {}
        //----------------------------------------------------------------------
        virtual String label () const { "][" };
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            a = serie.reflected ();
            b = -a; // short circuit rules
            return transformer.Vout();
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            Vin.Vs = wave;
            serie.incident (b);
        }
        //----------------------------------------------------------------------
    protected:
        WDF::VoltageSource<T>       Vin;
        WDF::Resistor<T>            Rload;
        WDF::Resistor<T>            Rterm;
        //----------------------------------------------------------------------
	NonIdealTransformer<T>	    transformer;
        //----------------------------------------------------------------------
        WDF::Serie<T>               serie;
        WDF::Parallel<T>            paral;
        //----------------------------------------------------------------------
        inline void wiring ()
        {
            transformer.connectChild (&Rload);
            paral.connect (&transformer, &Rterm);
            serie.connect (&paral,       &Vin);
        }
        //----------------------------------------------------------------------
};
//==============================================================================
//
//                  SIGNAL AMPLIFIER TRANSFORMERS
//                  -----------------------------
//
//               INPUT (Tx10)             OUTPUT (Tx20)
//      --------------------------------------------------------
//      Rc       10 kOhms                 10 kOhms
//      Lm       35.7 H                   35.7 H
//      Rp       10 Ohms                  5 Ohms
//      Lp       4 mH                     100 uH
//      Rs       50 Ohms                  50 Ohms
//      Ls       1 mH                     400 uH
//      Cw       210 pF                   1 pF
//      --------------------------------------------------------
//      Ratio    1:9                      9:1
//      --------------------------------------------------------
//
//      T101    - Fairchild signal amp input    (mono)
//      T102    - Fairchild signal amp input    (stereo)
//      --------------------------------------------------------
//      T201    - Fairchild signal amp output   (mono)
//      T202    - Fairchild signal amp output   (stereo)
//
//      --------------------------------------------------------
//
//==============================================================================
/**
template <typename T>
class Tx10 : public Transformer<T> { public:
      Tx10 (T Fs, String name = String::empty) :
      Transformer (Fs, 9.0, 600.0, 1360.0, 10e3, 35.7,
        10.0, 4e-3, 50.0, 1e-3, 1000e3, 210e-12, name)
{} };
//==============================================================================
template <typename T>
class Tx20 : public Transformer<T> { public:
      Tx20 (T Fs, String name = String::empty) :
      Transformer (Fs, 1.0/9.0, 600.0, 1360.0, 10e3, 35.7,
        5.0, 100e-6, 50.0, 400e-6, 1000e3, 1e-12, name)
{} };
*/
//==============================================================================
//
//                SIDECHAIN AMPLIFIER TRANSFORMERS
//                --------------------------------
//
//               INPUT (T10x)             OUTPUT (T20x)
//      --------------------------------------------------------
//      Rc       10 kOhms                 10 kOhms
//      Lm       35.7 H                   35.7 H
//      Rp       10 Ohms                  5 Ohms
//      Lp       4 mH                     100 uH
//      Rs       50 Ohms                  50 Ohms
//      Ls       1 mH                     400 uH
//      Cw       210 pF                   1 pF
//      --------------------------------------------------------
//      Ratio    1:9                      9:1
//      --------------------------------------------------------
//
//      T101    - Fairchild signal amp input    (mono)
//      T102    - Fairchild signal amp input    (stereo)
//      --------------------------------------------------------
//      T201    - Fairchild signal amp output   (mono)
//      T202    - Fairchild signal amp output   (stereo)
//
//      --------------------------------------------------------
//
//==============================================================================
/**
template <typename T>
class Tx30 : public Transformer<T> { public:
      Tx30 (T Fs, String name = String::empty) :
      Transformer (Fs, 9.0, 600.0, 1360.0, 10e3, 35.7,
        10.0, 4e-3, 50.0, 1e-3, 1000e3, 210e-12, name)
{} };
*/
//==============================================================================
} // namespace Wavechild670
//==============================================================================

//==============================================================================
// TRANSFORMERS
//==============================================================================
// Fairchild 670 Signal amp input           T101/201    600 ohms/50k ohms. Ratio 1+1:9+9.
//------------------------------------------------------------------------------
// Fairchild 670 Signal amp output          T102/202    600 ohms/60k ct/ Ratio 9+9:1+1
//------------------------------------------------------------------------------
// Fairchild 670 Control amp input          T103/203    600 ohms/170k ohms. Ratio 17+17:1+1
//------------------------------------------------------------------------------
// Fairchild 670 Control amp output         T104/204    600 ohms/10k ct/ Ratio 4:1
//                                                      Feedback winding ratio 9.5:1
//------------------------------------------------------------------------------
// Fairchild 670 Bias supply Transformer.   T301        50/60 Hz
//               Secondary windings:                    375-0-375 V at 200 mA,
//                                                          6.3V CT at 5A,
//                                                          5.0V at 2A
//------------------------------------------------------------------------------
// Fairchild 670 Mains Transformer.         T302        50/60 Hz
//               Secondary windings:                    26.8V (for Selenium bridge rectifier)
//                                              tapped at 24V (for Silicon bridge rectifier) at 200 mA.
//------------------------------------------------------------------------------
// Fairchild 670 Heater Transformer.        T303/304    50/60 Hz
//               Secondary windings:                    6.3 V CT at 3 A
//                                                      6.3 at 2.6A
//------------------------------------------------------------------------------
// Fairchild 670 Bias supply Choke.         L301        71 Ohms 10H at 200 mA
//------------------------------------------------------------------------------
// Fairchild 670 Bias supply Choke.         L302        85 Ohms  5H at 200 mA
//==============================================================================

//==============================================================================

//==============================================================================
template <typename T>
class SidechainAmplifier
{
    public:
        SidechainAmplifier (T Fs)
            : AC (0.5), DC (0.1)
        {}
        //----------------------------------------------------------------------
        void parameters (T ACThreshold, T DCThreshold)
        {
            DC = 12.2 * (DCThreshold + 0.1);
            AC = 0.5 * ACThreshold * ACThreshold;
        }
        //----------------------------------------------------------------------
        // Fairchild 670 Class-B Sidechain Amplifier model
	//----------------------------------------------------------------------
        virtual inline T process (T Vsc, T VlevelCap)
        {
	    //------------------------------------------------------------------
	    // AC Threshold Input Transformer
	    //------------------------------------------------------------------
            Vpot = AC * transformer.process (Vsc);
	    //------------------------------------------------------------------
            // DC Threshold Vsc Stage, 12AX7 amplifier
	    //------------------------------------------------------------------
            Vs1 = -6.0 * ((log(1.0 + exp( Vpot - DC)))
                        - (log(1.0 + exp(-Vpot - DC))));
	    //------------------------------------------------------------------
            // Drive stage, 12BH7 + 6973 amplifier stages
	    //------------------------------------------------------------------
	    Vdiff = fabs (hardclip (8.4 * Vs1, -100.0, 100.0)) - VlevelCap;
	    //------------------------------------------------------------------
            // The nominal output current through the bridge rectifier
            // is calculated using a diode model in series with a resistance.
	    //------------------------------------------------------------------
	    Inom = 0.000375 * log(1.0 + exp(((10.0 * Vdiff) / 0.6) - 10.0)) * 0.0125;
	    //------------------------------------------------------------------
	    // One side-saturation (does not saturate negatives)
	    //------------------------------------------------------------------
	    return Inom - 0.05 * log(1.0 + exp(((10.0 * Inom) / 0.5) - 10.0));
        }
	//----------------------------------------------------------------------
        inline T hardclip (T x, T min, T max) { return (x < min) ? min
                                                     : (x > max) ? max
                                                     :             x; }
        //----------------------------------------------------------------------
    protected:
        T DC, AC;
        T Vpot, Vs1, Vdiff, Inom;
        //----------------------------------------------------------------------
};

//==============================================================================
template <typename T>
class TubeStage : public WDF::OnePort<T>, WDF::NewtonRaphson<T>
{
    public:
        TubeStage (T Fs)
            : //----------------------------------------------------------------
              // Components
              //----------------------------------------------------------------
              Rout (600.0,      "Rout"),      // signal output
              Rsc (1000.0,      "Rsc"),       // sidechain input
              //----------------------------------------------------------------
              Ck (2.0*4e-6, Fs, "2C1"),       // cathode capacitor (twice)
              Vk (-3.1,  705.0, "Vbal R11"),  // cathode (balance)
              Vp (240.0,  33.0, "240V R12"),  // plate (power supply)
              //----------------------------------------------------------------
              Vgk (0.0), Iak (0.0), lVk (0.0) // lVk = last Vk (cathode voltage)
              //----------------------------------------------------------------
        {
        }
        //----------------------------------------------------------------------
        virtual String label () const { return "Tube"; }
        //----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            a = root.reflected (); return a;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T wave)
        {
            b = nonlinear (wave);
            root.incident (b);
            lVk = Vk.voltage (); // z-1
        }
        //----------------------------------------------------------------------
        inline T process (T Vgate)
        {
            reflected ();
            incident (Vgate);
            return transformer.Vout();
        }
        //----------------------------------------------------------------------
        // Fairchild 670 Class-A Signal Amplifier (with Push/Pull topology)
        //----------------------------------------------------------------------
        inline void wiring (WDF::OnePort<T>* coupled)
        {
            paral_O.connect (&Rout,    &Rsc);
            transfo.connect (&serie_T, &paral_O);
            serie_T.connect (&transfo, &Vp);
            //------------------------------------------------------------------
            serie_K.connect (&Ck,       coupled);
            paral_K.connect (&Vk,      &serie_K);
            //------------------------------------------------------------------
               root.connect (&serie_T, &paral_K);
        }
        //----------------------------------------------------------------------
    protected:
        //----------------------------------------------------------------------
        WDF::VoltageSource<T>   Vk;
        WDF::VoltageSource<T>   Vp;
        WDF::Capacitor<T>       Ck;
        WDF::Resistor<T>        Rout;
        WDF::Resistor<T>        Rsc;
        //----------------------------------------------------------------------
        NonIdealTransformer<T>  transfo;
        //----------------------------------------------------------------------
        WDF::Serie<T>           root;
        WDF::Serie<T>           serie_K;
        WDF::Serie<T>           serie_T;
        WDF::Parallel<T>        paral_O;
        WDF::Parallel<T>        paral_K;
        //----------------------------------------------------------------------
    private:
        //----------------------------------------------------------------------
        // nonlinearity is isolated at trunk (wdf requirement)
        //----------------------------------------------------------------------
        inline T nonlinear (T Vgate)
        {
            Vgk = Vgate - lVk;              // grid-cathode voltage
            //------------------------------------------------------------------
            Iak = 0.0;                      // computed by solve()
            T Vak = solve ();               // Newton/Raphson iterations
            //------------------------------------------------------------------
            return Vak - root.R()*Iak;      // estimate reflected
        }
        //----------------------------------------------------------------------
        // implicit equation will be evaluate by Newton/Raphson solver
        //----------------------------------------------------------------------
        virtual inline T evaluate (T Vak)
        {
            Iak = Ia (Vgk, Vak) * NTI;
            return Vak + R*Iak - a;         // [ Vak + R*Iak - a = 0 ]
        }
        //----------------------------------------------------------------------
        // GE 6386 Remote Cutoff Triode
        //----------------------------------------------------------------------
        // The model parameters were calculated using Levenberg-Marquardt least
        // squares estimation and hand tuning to fit the 6386 characteristics
        // as given in the General Electric 6386 datasheet.
        // by Peter Raffensperger (2012)
        //----------------------------------------------------------------------
        inline T Ia (T Vgk, T Vak)          // Ia = anode current (in amps)
        {
            if (Vak < 0.0) Vak = 0.0;
            if (Vgk > 0.0) Vgk = 0.0;
            //------------------------------------------------------------------
            return (3.981e-8 * pow(Vak, 2.383))
                 / (pow((0.5 - 0.1*Vgk), 1.8)
                 * (0.5 + exp((-0.03922*Vak)
                 - (0.2*Vgk))));
        }
        //----------------------------------------------------------------------
    private:
        T Vgk, Iak, lVk;
        //----------------------------------------------------------------------
};

template <typename T>
class SidechainAmplifier
{
    public:
        SidechainAmplifier (T Fs)
            : AC (0.5), DC (0.1)
        {}
        //----------------------------------------------------------------------
        void parameters (T ACThreshold, T DCThreshold)
        {
            DC = 12.2 * (DCThreshold + 0.1);
            AC = 0.5 * ACThreshold * ACThreshold;
        }
        //----------------------------------------------------------------------
        // Fairchild 670 Class-B Sidechain Amplifier model
	//----------------------------------------------------------------------
        virtual inline T process (T Vsc, T VlevelCap)
        {
	    //------------------------------------------------------------------
	    // AC Threshold Input Transformer
	    //------------------------------------------------------------------
            Vpot = AC * transformer.process (Vsc);
	    //------------------------------------------------------------------
            // DC Threshold Vsc Stage, 12AX7 amplifier
	    //------------------------------------------------------------------
            Vs1 = -6.0 * ((log(1.0 + exp( Vpot - DC)))
                        - (log(1.0 + exp(-Vpot - DC))));
	    //------------------------------------------------------------------
            // Drive stage, 12BH7 + 6973 amplifier stages
	    //------------------------------------------------------------------
	    Vdiff = fabs (hardclip (8.4 * Vs1, -100.0, 100.0)) - VlevelCap;
	    //------------------------------------------------------------------
            // The nominal output current through the bridge rectifier
            // is calculated using a diode model in series with a resistance.
	    //------------------------------------------------------------------
	    Inom = 0.000375 * log(1.0 + exp(((10.0 * Vdiff) / 0.6) - 10.0)) * 0.0125;
	    //------------------------------------------------------------------
	    // One side-saturation (does not saturate negatives)
	    //------------------------------------------------------------------
	    return Inom - 0.05 * log(1.0 + exp(((10.0 * Inom) / 0.5) - 10.0));
        }
	//----------------------------------------------------------------------
        inline T hardclip (T x, T min, T max) { return (x < min) ? min
                                                     : (x > max) ? max
                                                     :             x; }
        //----------------------------------------------------------------------
    protected:
        T DC, AC;
        T Vpot, Vs1, Vdiff, Inom;
        //----------------------------------------------------------------------
};

//==============================================================================
template <typename T> struct UnitDelay { T a, b; };
//==============================================================================
template <typename T>
class BidirectionnalUnitDelay
{
    public:
        //----------------------------------------------------------------------
        void process ()
        {
            unit1.b = unit2.a;
            unit2.b = unit1.a;
        }
        //----------------------------------------------------------------------
        UnitDelay<T> unit1;
        UnitDelay<T> unit2;
        //----------------------------------------------------------------------
};
//==============================================================================
template <typename T>
class TransformerInputCircuit
{
	public:
		TransformerInputCircuit ()
		{
		}
		//----------------------------------------------------------------------
};
//==============================================================================
template <typename T>
class SignalAmplifier : public WDF::OnePort<T>
{
    public:
        SignalAmplifier (T Fs)
            : //----------------------------------------------------------------
              WDF::OnePort<T> (1.0),
		cathodeTocathode (new BidirectionnalUnitDelay<T>()),
		transformer (new InputCoupledTransformer<T>()),
		push (new TubeStage<T>(Fs)),
		pull (new TubeStage<T>(Fs)),
		VgateBias (-7.2)
              //----------------------------------------------------------------
        {
	    push->wiring (pull);
	    pull->wiring (push);
	}
        //----------------------------------------------------------------------
        virtual String label () const { return "Amp"; }
	//----------------------------------------------------------------------
        virtual inline T process (T Vin, T VlevelCap)
        {
            T Vgate = transformer->process (Vin);
            T VoutPush = push->process (VgateBias - VlevelCap + Vgate);
            T VoutPull = pull->process (VgateBias - VlevelCap + Vgate);
            cathodeTocathode->process ();
            return VoutPush - VoutPull;
        }
	//----------------------------------------------------------------------
        virtual inline T reflected ()
        {
            b = 0.0; return b;
        }
        //----------------------------------------------------------------------
        virtual inline void incident (T value)
        {
	    a = value;
        }
        //----------------------------------------------------------------------
    protected:
        //----------------------------------------------------------------------
        // Fairchild 670 Class-A Signal Amplifier model
        //----------------------------------------------------------------------
        ScopedPointer<BidirectionnalUnitDelay<T>> cathodeTocathode;
        ScopedPointer<InputCoupledTransformer<T>> transformer;
        ScopedPointer<TubeStage<T>> push; // GE 6386
        ScopedPointer<TubeStage<T>> pull; // GE 6386
	T VgateBias;
        //----------------------------------------------------------------------
};
//==============================================================================

//==============================================================================
#define SQRT_2 sqrt(2.0)
//==============================================================================
template <typename T>
class StereoProcessor
{
    public:
        StereoProcessor ()
            : Fs (44100.0),       gain (1.0),
              //-------------------------
                   A (0.0),          B (0.0),
                capA (0.0),       capB (0.0),
              levelA (1.0),     levelB (1.0),
          thresholdA (1.0), thresholdB (1.0),
                 tcA (2),          tcB (2),
              //-------------------------
              hardclipout (true),
                 feedback (false),
                  midside (false),
                   linked (true)
        {}
        //----------------------------------------------------------------------
        void init (T sampleRate)
        {
            Fs = sampleRate;
            //------------------------------------------------------------------
            signalAmpA = new SignalAmplifier<double> (Fs);
            signalAmpB = new SignalAmplifier<double> (Fs);
            //------------------------------------------------------------------
            sidechainAmpA = new SidechainAmplifier<double> (Fs);
            sidechainAmpB = new SidechainAmplifier<double> (Fs);
            //------------------------------------------------------------------
            timeConstantA = new LevelTimeConstant<double> (Fs);
            timeConstantB = new LevelTimeConstant<double> (Fs);
            //------------------------------------------------------------------
            timeConstantA->parameters (Fs, tcA);
            timeConstantB->parameters (Fs, tcB);
            //------------------------------------------------------------------
            capA = A = 0.0;
            capB = B = 0.0;
            //------------------------------------------------------------------
            warmup ();
        }
        //----------------------------------------------------------------------
        void parameters (const int tA, const int tB)
        {
            tcA = tA; timeConstantA->parameters (Fs, tcA);
            tcB = tB; timeConstantB->parameters (Fs, tcB);
        }
        //----------------------------------------------------------------------
        inline T sidechain (T VscA, T VscB)
        {
            T IscA = sidechainAmpA->process (VscA, capA);
            T IscB = sidechainAmpB->process (VscB, capB);

            if (linked)
            {
                T IscT = (IscA + IscB) * 0.5;
                T Ax = timeConstantA->process (IscT);
                T Bx = timeConstantB->process (IscT);
                capA =
                capB = (Ax + Bx) * 0.5;
            }
            else
            {
                capA = timeConstantA->process (IscA);
                capB = timeConstantA->process (IscB);
            }
        }
        //----------------------------------------------------------------------
        inline T hardclip (T x, T min, T max) { return (x < min) ? min
                                                     : (x > max) ? max
                                                     :             x; }
        //----------------------------------------------------------------------
        inline void process (float *left, float *right)
        {
            A = (midside) ? (T)((left[0]  + left[0] ) / SQRT_2) : left[0];
            B = (midside) ? (T)((right[0] - right[0]) / SQRT_2) : right[0];

            A *= levelA;
            B *= levelB;

            if (!feedback) sidechain (A, B);

            A = signalAmpA->process (A, capA);
            B = signalAmpB->process (B, capB);

            if ( feedback) sidechain (A, B);

            A = (midside) ? (A + B) / SQRT_2 : A;
            B = (midside) ? (A - B) / SQRT_2 : B;

            A *= gain;
            B *= gain;

            A = (hardclipout) ? hardclip(A, -1.0, 1.0) : A;
            B = (hardclipout) ? hardclip(B, -1.0, 1.0) : B;

            left[0]  = (float)A;
            right[0] = (float)B;
        }
        //----------------------------------------------------------------------
        void warmup (T timeInSec = 0.5)
        {
            long i, samples = (long) (timeInSec*Fs)/2;
            i = 0; for (; i < samples; ++i) { signalAmpA->process (0.0, capA);
                                              signalAmpB->process (0.0, capB); }
            i = 0; for (; i < samples; ++i) { T VscA = signalAmpA->process (0.0, capA);
                                              T VscB = signalAmpB->process (0.0, capB);
                                              sidechain (VscA, VscB); }
        }
        //----------------------------------------------------------------------
        T Fs; // samplerate
        //----------------------------------------------------------------------
        int tcA, tcB;
        bool hardclipout, midside, linked, feedback;
        T A, B, capA, capB, levelA, levelB, thresholdA, thresholdB, gain;
        //----------------------------------------------------------------------
        ScopedPointer<SignalAmplifier<T>>    signalAmpA;
        ScopedPointer<SignalAmplifier<T>>    signalAmpB;
        //----------------------------------------------------------------------
        ScopedPointer<LevelTimeConstant<T>>  timeConstantA;
        ScopedPointer<LevelTimeConstant<T>>  timeConstantB;
        //----------------------------------------------------------------------
        ScopedPointer<SidechainAmplifier<T>> sidechainAmpA;
        ScopedPointer<SidechainAmplifier<T>> sidechainAmpB;
        //----------------------------------------------------------------------
};
//==============================================================================
#undef SQRT_2
//==============================================================================


//==============================================================================
} // namespace Wavechild670
//==============================================================================