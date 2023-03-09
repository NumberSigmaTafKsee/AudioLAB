#pragma once

namespace Loudia
{
        class Meddis {
    protected:
        // Internal fixed constants
        static const Real M;
        static const Real A;
        static const Real B;
        static const Real g;
        static const Real y;
        static const Real l;
        static const Real r;
        static const Real x;
        static const Real h;

        // Internal parameters
        Real _sampleRate;
        int _channels;
        bool _substractSpont;
        
        Real dt;
        Real gdt;
        Real ydt;
        Real ldt;
        Real rdt;
        Real xdt;

        // Internal variables
        MatrixXR kt;
        MatrixXR spont;

        MatrixXR c;
        MatrixXR q;
        MatrixXR w;

    public:
        Meddis(Real sampleRate, int channels, bool substractSpont = true);

        ~Meddis();

        void setup();

        void process(const MatrixXR& data, MatrixXR* output);

        int channels() const;

        Real sampleRate() const;

        void reset();
    };

    // Internal fixed constants
    const Real Meddis::M = 1.0;
    const Real Meddis::A = 5.0;
    const Real Meddis::B = 300.0;
    const Real Meddis::g = 2000.0;
    const Real Meddis::y = 5.05;
    const Real Meddis::l = 2500.0;
    const Real Meddis::r = 6580.0;
    const Real Meddis::x = 66.31;
    const Real Meddis::h = 50000.0;

    Meddis::Meddis(Real sampleRate, int channels, bool substractSpont) : 
    _sampleRate( sampleRate ),
    _channels( channels ),
    _substractSpont( substractSpont ),
    kt( 1, channels ), 
    spont( 1, channels ), 
    c( 1, channels ), 
    q( 1, channels ), 
    w( 1, channels ) 
    {
    LOUDIA_DEBUG("MEDDIS: Constructor sampleRate:" << sampleRate << 
            ", channels:" << channels << 
            ", substractSpont:" << substractSpont);
    
    setup();
    }

    Meddis::~Meddis() {}

    void Meddis::setup(){
    LOUDIA_DEBUG("MEDDIS: Setting up...");

    // Configure the internal constants
    dt = 1./_sampleRate;
    gdt = g*dt;
    ydt = y*dt;
    ldt = l*dt;
    rdt = r*dt;
    xdt = x*dt;

    // Prepare the buffers
    reset();
    LOUDIA_DEBUG("MEDDIS: Set up finished.");
    }

    void Meddis::process(const MatrixXR& samples, MatrixXR* output){
    // Process will be called with a matrix where columns will be channels
    // and rows will be the time axis
    MatrixXR  row, limitedSt, replenish, eject, loss, reuptake, reprocess;

    (*output).resize(samples.rows(), _channels);

    for (int i = 0; i < samples.rows(); ++i) {
        row = samples.row(i);
    
        limitedSt = (row.array() + A).clipUnder();

        kt = (limitedSt * gdt).array() / (limitedSt.array() + B);

        replenish = ydt * ((-q).array() + M).clipUnder();
        eject = kt.array() * q.array();
        loss = ldt * c;
        reuptake = rdt * c;
        reprocess = xdt * w;
        
        q += replenish - eject + reprocess;
        c += eject - loss - reuptake;
        w += reuptake - reprocess;

        // Now iterate through each time slice of the data.  Use the
        // max function to implement the "if (0>" test.
        (*output).row(i) = h * c;
        
        if(_substractSpont){
        (*output).row(i) = ((*output).row(i) - spont.row(0)).array().clipUnder(0.0);
        }
    } // for each row
    }

    void Meddis::reset(){
    // Initial values
    LOUDIA_DEBUG("MEDDIS: Resetting...");

    kt = MatrixXR::Constant(1, _channels, g * A / (A + B));
    spont = kt * (M * y);
    spont = spont.array() / ( (kt*l).array() + (y * ( l + r )) );

    c = spont;
    q = (c * ( l + r )).array() / kt.array();
    w = c * r / x;
    }


    int Meddis::channels() const {
    return _channels;
    }

    Real Meddis::sampleRate() const {
    return _sampleRate;
    }
}