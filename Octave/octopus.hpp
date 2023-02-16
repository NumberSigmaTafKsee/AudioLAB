#pragma once
#include <iostream>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
#include <Eigen/Core>

//#include "carlo_mkl.hpp"
//#include "carlo_samples.hpp"

using ArrayXf = Array<float>;
using ArrayXd = Array<double>;
using ArrayXcf = Array<std::complex<float>>;
using ArrayXcd = Array<std::complex<double>>;
using VectorXf = FloatRowVector;
using VectorXd = RowVector;
using VectorXcf= FloatComplexRowVector;
using VectorXcd= ComplexRowVector;
using ColVectorXf = FloatColumnVector;
using ColVectorXd = ColumnVector;
using ColVectorXcf= FloatComplexColumnVector;
using ColVectorXcd= ComplexColumnVector;
using MatrixXf = FloatMatrix;
using MatrixXd = Matrix;
using MatrixXcf= FloatComplexMatrix;
using MatrixXcd= ComplexMatrix;
using Value=octave_value;
using ValueList=octave_value_list;

#include "octopus_rowvector.hpp"
#include "octopus_colvector.hpp"
#include "octopus_matrix.hpp"
//#include "octopus_octavate.hpp"


namespace Octopus
{   
    /*
    struct Application : public octave::application
    {
        Application() {
            forced_interactive(true);
        }
        int execute() {
            return 0;
        }
    };
    */
    struct OctopusValue : public octave_value
    {
        OctopusValue() = default;
        OctopusValue(const octave_value & v) : octave_value(v) {}
        
        OctopusValue(double v) : octave_value(v) {}

        OctopusValue(const ArrayXf& v) : octave_value(v) {}
        OctopusValue(const ArrayXd& v) : octave_value(v) {}
        OctopusValue(const ArrayXcf& v) : octave_value(v) {}
        OctopusValue(const ArrayXcd& v) : octave_value(v) {}

        OctopusValue(const VectorXf& v) : octave_value(v) {}
        OctopusValue(const VectorXd& v) : octave_value(v) {}
        OctopusValue(const VectorXcf& v) : octave_value(v) {}
        OctopusValue(const VectorXcd& v) : octave_value(v) {}

        OctopusValue(const ColVectorXf& v) : octave_value(v) {}
        OctopusValue(const ColVectorXd& v) : octave_value(v) {}
        OctopusValue(const ColVectorXcf& v) : octave_value(v) {}
        OctopusValue(const ColVectorXcd& v) : octave_value(v) {}

        OctopusValue(const MatrixXf& v) : octave_value(v) {}
        OctopusValue(const MatrixXd& v) : octave_value(v) {}
        OctopusValue(const MatrixXcf& v) : octave_value(v) {}
        OctopusValue(const MatrixXcd& v) : octave_value(v) {}
        
        /*        
        OctopusValue(const Casino::sample_vector<float>& v) : octave_value(Octavate(v)) {}
        OctopusValue(const Casino::sample_vector<double>& v) : octave_value(Octavate(v)) {}
        OctopusValue(const Casino::complex_vector<float>& v) : octave_value(Octavate(v)) {}
        OctopusValue(const Casino::complex_vector<double>& v) : octave_value(Octavate(v)) {}

        OctopusValue(const Casino::sample_matrix<float>& v) : octave_value(Octavate(v)) {}
        OctopusValue(const Casino::sample_matrix<double>& v) : octave_value(Octavate(v)) {}
        OctopusValue(const Casino::complex_matrix<float>& v) : octave_value(Octavate(v)) {}
        OctopusValue(const Casino::complex_matrix<double>& v) : octave_value(Octavate(v)) {}
        */

        double getScalarValue() {
            return this->scalar_value();
        }
        
        OctopusRowVectorXf getFloatRowVector() {
            return OctopusRowVectorXf(this->float_row_vector_value());
        }
        OctopusColVectorXf getFloatColVector() {
            return OctopusColVectorXf(this->float_column_vector_value());
        }
        OctopusRowVectorXcf getFloatComplexRowVector() {
            return OctopusRowVectorXcf(this->float_complex_row_vector_value());
        }
        OctopusColVectorXcf getFloatComplexColVector() {
            return OctopusColVectorXcf(this->float_complex_column_vector_value());
        }
        OctopusRowVectorXd getRowVector() {
            return OctopusRowVectorXd(this->row_vector_value());
        }
        OctopusColVectorXd getColVector() {
            return OctopusColVectorXd(this->column_vector_value());
        }
        OctopusRowVectorXcd getComplexRowVector() {
            return OctopusRowVectorXcd(this->complex_row_vector_value());
        }
        OctopusColVectorXcd getComplexColVector() {
            return OctopusColVectorXcd(this->complex_column_vector_value());
        }

        OctopusMatrixXf getFloatMatrix() {
            return OctopusMatrixXf(this->float_matrix_value());
        }
        OctopusMatrixXd getMatrix() {
            return OctopusMatrixXd(this->matrix_value());
        }
        OctopusMatrixXcf getFloatComplexMatrix() {
            return OctopusMatrixXcf(this->float_complex_matrix_value());
        }
        OctopusMatrixXcd getComplexMatrix() {
            return OctopusMatrixXcd(this->complex_matrix_value());
        }

        
    };

    
    
    struct OctopusValueList 
    {
        octave_value_list vlist;

        OctopusValueList() = default;
        OctopusValueList(const octave_value_list& v) : vlist(v) {}        
        OctopusValueList(const OctopusValueList& v) : vlist(v.vlist) {}
        ~OctopusValueList() = default;

        OctopusValue get(size_t i) {
            return OctopusValue((*this)(i));
        }
        void set(size_t i, const OctopusValue & v) {
            vlist(i) = v;
        }
        void set(size_t i, const double val) {
            OctopusValue v(val);
            vlist(i) = v;
        }
        void set(size_t i, const std::string& val) {
            OctopusValue v(val);
            vlist(i) = v;
        }
        
        void clear() {
            vlist.clear();
        }
        Value& operator[](size_t i) {
            return vlist(i);
        }
        Value& operator()(size_t i) {
            return vlist(i);
        }
        
        OctopusValueList& operator = (const OctopusValueList & v) {
            vlist = v.vlist;
            return *this;
        }
        OctopusValueList& operator = (const octave_value_list & v) {
            vlist = v;
            return *this;
        }
        OctopusValue __getitem__(size_t i) { return get(i); }

        void __setitem__(size_t i, const OctopusValue& v) { set(i,v); }

        void __setitem__(size_t i, const double& v) { set(i,Value(v)); }

        void __setitem__(size_t i, const ArrayXf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const ArrayXd& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const ArrayXcf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const ArrayXcd& v) { set(i,Value(v)); }
        
        void __setitem__(size_t i, const VectorXf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const VectorXd& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const VectorXcf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const VectorXcd& v) { set(i,Value(v)); }
        
        void __setitem__(size_t i, const ColVectorXf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const ColVectorXd& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const ColVectorXcf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const ColVectorXcd& v) { set(i,Value(v)); }

        void __setitem__(size_t i, const MatrixXf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const MatrixXd& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const MatrixXcf& v) { set(i,Value(v)); }
        void __setitem__(size_t i, const MatrixXcd& v) { set(i,Value(v)); }
        
        /*
        void __setitem__(size_t i, const Casino::sample_vector<float>& v) { set(i,Value(Octavate(v))); }
        void __setitem__(size_t i, const Casino::sample_vector<double>& v) { set(i,Value(Octavate(v))); }
        void __setitem__(size_t i, const Casino::complex_vector<float>& v) { set(i,Value(Octavate(v))); }
        void __setitem__(size_t i, const Casino::complex_vector<double>& v) { set(i,Value(Octavate(v))); }

        void __setitem__(size_t i, const Casino::sample_matrix<float>& v) { set(i,Value(Octavate(v))); }
        void __setitem__(size_t i, const Casino::sample_matrix<double>& v) { set(i,Value(Octavate(v))); }
        void __setitem__(size_t i, const Casino::complex_matrix<float> & v) { set(i,Value(Octavate(v))); }
        void __setitem__(size_t i, const Casino::complex_matrix<double>& v) { set(i,Value(Octavate(v))); }        
        */
    };

    struct OctaveInterpreter
    {   
        octave::interpreter *interpreter;
        //Application pita;    

        OctaveInterpreter() {                  
            interpreter = new octave::interpreter();
            interpreter->interactive(false);
            interpreter->initialize_history(false);       
            interpreter->initialize();            
            interpreter->execute();
            std::string path = "Matlab";
            octave_value_list p;
            p(0) = path;
            octave_value_list o1 = interpreter->feval("addpath", p, 1);            
            run_script("startup.m");
        }
        ~OctaveInterpreter()
        {
            if(interpreter) delete interpreter;
        }
        
        void run_script(const std::string& s) {
            octave::source_file(s);
        }
        
        OctopusValueList eval_string(const std::string& func, bool silent=false, int noutputs=1)
        {          
            octave_value_list out =interpreter->eval_string(func.c_str(), silent, noutputs);
            return OctopusValueList(out);
        }
        OctopusValueList eval(const std::string& func, const OctopusValueList& inputs, int noutputs=1)
        {          
            octave_value_list out = interpreter->feval(func.c_str(), inputs.vlist, noutputs);
            return OctopusValueList(out);
        }

        OctopusValueList operator()(const std::string& func, const OctopusValueList& inputs, int noutputs=1)
        {
            return eval(func,inputs,noutputs);
        }
        
        void createVar(const std::string& name, const OctopusValue& v, bool global=true)
        {
            interpreter->install_variable(name,v,global);
        }        
        OctopusValue getGlobalVar(const std::string& name)
        {    
            return interpreter->global_varval(name);
        }
        void setGlobalVar(const std::string& name, const OctopusValue& v)
        {
            //interpreter->set_global_value(name,v);
        }
        OctopusValue getVarVal(const std::string& name) {
            return interpreter->varval(name);
        }
        void assign(const std::string& name, const OctopusValue& v) {
            interpreter->assign(name,v);
        }

    };

    struct OctopusVar
    {
        OctaveInterpreter& interp;
        std::string var;
        bool _global;
        OctopusVar(OctaveInterpreter& i, const std::string& name, const Value & v, bool global=true)
        : interp(i),var(name),_global(global)
        {
            interp.createVar(var,v,global);
        }

        Value getValue() { 
            if(_global) return interp.getGlobalVar(var); 
            return interp.getVarVal(var);
        }
        void setValue(const Value& v) { 
            if(_global) interp.setGlobalVar(var,v); 
            else interp.assign(var,v);    
        }
    };
    
    
    struct OctopusFunction
    {
        std::string code;
        std::string name;
        OctaveInterpreter& interp;
        
        OctopusFunction(OctaveInterpreter& i, const std::string& n, const std::string& c)
        : interp(i),code(c),name(n) 
        {                        
            interp.interpreter->eval(code,0);
        }
        OctopusValueList operator()(OctopusValueList & inputs, int numOut)
        {     
            return interp.eval(name.c_str(),inputs,numOut);
        }
    };  

    struct Function
    {
        
        std::string name;
                
        Function(const std::string& f) : name(f) {}
                
        OctopusValueList operator()(const OctopusValueList & input, int num_outputs=1)
        {
            return eval(input,num_outputs);
        }
        
        OctopusMatrixXf operator()(const OctopusMatrixXf & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).float_matrix_value();
        }
        OctopusMatrixXd operator()(const OctopusMatrixXd & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).matrix_value();
        }
        OctopusMatrixXcf operator()(const OctopusMatrixXcf & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).float_complex_matrix_value();
        }
        OctopusMatrixXcd operator()(const OctopusMatrixXcd & m, int num_outputs=1)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).complex_matrix_value();
        }
        OctopusRowVectorXf operator()(const OctopusRowVectorXf & m)  
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).float_row_vector_value();
        }
        OctopusRowVectorXd operator()(const OctopusRowVectorXd & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).row_vector_value();
        }
        OctopusRowVectorXcf operator()(const OctopusRowVectorXcf & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).float_complex_row_vector_value();
        }
        OctopusRowVectorXcd operator()(const OctopusRowVectorXcd & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).complex_row_vector_value();
        }
        OctopusColVectorXf operator()(const OctopusColVectorXf & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).float_column_vector_value();
        }
        OctopusColVectorXd operator()(const OctopusColVectorXd & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).column_vector_value();
        }
        OctopusColVectorXcf operator()(const OctopusColVectorXcf & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).float_complex_column_vector_value();
        }
        OctopusColVectorXcd operator()(const OctopusColVectorXcd & m)
        {
            OctopusValueList input;
            input(0) = m;
            return eval(input)(0).complex_column_vector_value();
        }

        OctopusValueList eval(const OctopusValueList &input, int numOutputs=1);
        
        /*        
        AudioDSP::sample_vector<float> eval(const AudioDSP::sample_vector<float> & input, int numOutputs=1)
        {
            AudioDSP::sample_vector<float> r;
            VectorXf x(input.size());
            for(size_t i = 0; i < input.size(); i++) x(i) = input[i];
            ValueList v;
            v(0) = x;
            octave::feval(name.c_str(),v,numOutputs);
            if(numOutputs > 0)
            {
                x = v(0).float_row_vector_value();
                r.resize(x.size(1));
                for(size_t i = 0; x.size(1); i++) 
                    r[i] = x(i);
            }
            return r;
        }        
        AudioDSP::sample_vector<double> eval(const AudioDSP::sample_vector<double> & input, int numOutputs=1)
        {
            AudioDSP::sample_vector<double> r;
            VectorXd x(input.size());
            for(size_t i = 0; i < input.size(); i++) x(i) = input[i];
            ValueList v;
            v(0) = x;
            octave::feval(name.c_str(),v,numOutputs);
            if(numOutputs > 0)
            {
                x = v(0).row_vector_value();
                r.resize(x.size(1));
                for(size_t i = 0; x.size(1); i++) 
                    r[i] = x(i);
            }
            return r;
        }
        AudioDSP::sample_matrix<float> eval(const AudioDSP::sample_matrix<float> & input, int numOutputs=1)
        {
            AudioDSP::sample_matrix<float> r;
            MatrixXf x(input.rows(),input.cols());
            for(size_t i = 0; i < input.rows(); i++) 
            for(size_t j = 0; j < input.cols(); j++)
                x(i,j) = input(i,j);
            ValueList v;
            v(0) = x;
            octave::feval(name.c_str(),v,numOutputs);
            if(numOutputs > 0)
            {
                x = v(0).float_matrix_value();
                r.resize((size_t)x.rows(),(size_t)x.cols());
                for(size_t i = 0; x.rows(); i++) 
                for(size_t j = 0; x.cols(); j++) 
                    r(i,j) = x(i,j);
            }
            return r;
        }        
        AudioDSP::sample_matrix<double> eval(const AudioDSP::sample_matrix<double> & input, int numOutputs=1)
        {
            AudioDSP::sample_matrix<double> r;
            MatrixXd x(input.rows(),input.cols());
            for(size_t i = 0; i < input.rows(); i++) 
            for(size_t j = 0; j < input.cols(); j++)
                x(i,j) = input(i,j);
            ValueList v;
            v(0) = x;
            octave::feval(name.c_str(),v,numOutputs);
            if(numOutputs > 0)
            {
                x = v(0).matrix_value();
                r.resize((size_t)x.rows(),(size_t)x.cols());
                for(size_t i = 0; x.rows(); i++) 
                for(size_t j = 0; x.cols(); j++) 
                    r(i,j) = x(i,j);
            }
            return r;
        }
        */
        /*
        #ifdef SWIG
        %extend {            
            OctopusValueList __call__(const OctopusValueList & v)
            {
                return octave::feval(v,1);
            }
            AudioDSP::sample_vector<float> __call__(const AudioDSP::sample_vector<float> & input) {
                return octave::feval(input,1);
            }
            AudioDSP::sample_vector<double> __call__(const AudioDSP::sample_vector<double> & input) {
                return octave::feval(input,1);
            }
            AudioDSP::sample_matrix<float> __call__(const AudioDSP::sample_matrix<float> & input) {
                return octave::feval(input,1);
            }
            AudioDSP::sample_matrix<double> __call__(const AudioDSP::sample_matrix<double> & input) {
                return octave::feval(input,1);
            }
        }
        #endif
        */
    };


extern Function octave_fft;
extern Function octave_ifft;
extern Function octave_fft2;
extern Function octave_ifft2;
extern Function octave_fftconv;
extern Function octave_fftfilt;
extern Function octave_fftn;
extern Function octave_fftshift;
extern Function octave_fftw;
extern Function octave_ifftn;
extern Function octave_ifftshift;
extern Function octave_ifht;
extern Function octave_ifourier;
extern Function octave_ifwht;
extern Function octave_ifwt;
extern Function octave_ifwt2;
extern Function octave_buffer;
extern Function octave_chirp;
extern Function octave_cmorwavf;
extern Function octave_gauspuls;
extern Function octave_gmonopuls;
extern Function octave_mexihat;
extern Function octave_meyeraux;
extern Function octave_morlet;
extern Function octave_pulstran;
extern Function octave_rectpuls;
extern Function octave_sawtooth;
extern Function octave_shanwavf;
extern Function octave_shiftdata;
extern Function octave_sigmoid_train;
extern Function octave_specgram;
extern Function octave_square;
extern Function octave_tripuls;
extern Function octave_udecode;
extern Function octave_uencoder;
extern Function octave_unshiftdata;
extern Function octave_findpeaks;
extern Function octave_peak2peak;
extern Function octave_peak2rms;
extern Function octave_rms;
extern Function octave_rssq;
extern Function octave_cconv;
extern Function octave_convmtx;
extern Function octave_wconv;
extern Function octave_xcorr;
extern Function octave_xcorr2;
extern Function octave_xcov;
extern Function octave_filtfilt;
extern Function octave_fltic;
extern Function octave_medfilt1;
extern Function octave_movingrms;
extern Function octave_sgolayfilt;
extern Function octave_sosfilt;
extern Function octave_freqs;
extern Function octave_freqs_plot;
extern Function octave_freqz;
extern Function octave_freqz_plot;
extern Function octave_impz;
extern Function octave_zplane;
extern Function octave_filter;
extern Function octave_filter2;
extern Function octave_fir1;
extern Function octave_fir2;
extern Function octave_firls;
extern Function octave_sinc;
extern Function octave_unwrap;
extern Function octave_bartlett;
extern Function octave_blackman;
extern Function octave_blackmanharris;
extern Function octave_blackmannuttal;
extern Function octave_dftmtx;
extern Function octave_hamming;
extern Function octave_hann;
extern Function octave_hanning;
extern Function octave_pchip;
extern Function octave_periodogram;
extern Function octave_sinetone;
extern Function octave_sinewave;
extern Function octave_spectral_adf;
extern Function octave_spectral_xdf;
extern Function octave_spencer;
extern Function octave_stft;
extern Function octave_synthesis;
extern Function octave_yulewalker;
extern Function octave_polystab;
extern Function octave_residued;
extern Function octave_residuez;
extern Function octave_sos2ss;
extern Function octave_sos2tf;
extern Function octave_sos2zp;
extern Function octave_ss2tf;
extern Function octave_ss2zp;
extern Function octave_tf2sos;
extern Function octave_tf2ss;
extern Function octave_tf2zp;
extern Function octave_zp2sos;
extern Function octave_zp2ss;
extern Function octave_zp2tf;
extern Function octave_besselap;
extern Function octave_besself;
extern Function octave_bilinear;
extern Function octave_buttap;
extern Function octave_butter;
extern Function octave_buttord;
extern Function octave_cheb;
extern Function octave_cheb1ap;
extern Function octave_cheb1ord;
extern Function octave_cheb2ap;
extern Function octave_cheb2ord;
extern Function octave_chebywin;
extern Function octave_cheby1;
extern Function octave_cheby2;
extern Function octave_ellip;
extern Function octave_ellipap;
extern Function octave_ellipord;
extern Function octave_impinvar;
extern Function octave_ncauer;
extern Function octave_pei_tseng_notch;
extern Function octave_sftrans;
extern Function octave_cl2bp;
extern Function octave_kaiserord;
extern Function octave_qp_kaiser;
extern Function octave_remez;
extern Function octave_sgplay;
extern Function octave_bitrevorder;
extern Function octave_cceps;
extern Function octave_cplxreal;
extern Function octave_czt;
extern Function octave_dct;
extern Function octave_dct2;
extern Function octave_dctmtx;
extern Function octave_digitrevorder;
extern Function octave_dst;
extern Function octave_dwt;
extern Function octave_rceps;
extern Function octave_ar_psd;
extern Function octave_cohere;
extern Function octave_cpsd;
extern Function octave_csd;
extern Function octave_db2pow;
extern Function octave_mscohere;
extern Function octave_pburg;
extern Function octave_pow2db;
extern Function octave_pwelch;
extern Function octave_pyulear;
extern Function octave_tfe;
extern Function octave_tfestimate;
extern Function octave___power;
extern Function octave_barthannwin;
extern Function octave_bohmanwin;
extern Function octave_boxcar;
extern Function octave_flattopwin;
extern Function octave_chebwin;
extern Function octave_gaussian;
extern Function octave_gausswin;
extern Function octave_kaiser;
extern Function octave_nuttalwin;
extern Function octave_parzenwin;
extern Function octave_rectwin;
extern Function octave_tukeywin;
extern Function octave_ultrwin;
extern Function octave_welchwin;
extern Function octave_window;
extern Function octave_arburg;
extern Function octave_aryule;
extern Function octave_invfreq;
extern Function octave_invfreqz;
extern Function octave_invfreqs;
extern Function octave_levinson;
extern Function octave_data2fun;
extern Function octave_decimate;
extern Function octave_interp;
extern Function octave_resample;
extern Function octave_upfirdn;
extern Function octave_upsample;
extern Function octave_clustersegment;
extern Function octave_fracshift;
extern Function octave_marcumq;
extern Function octave_primitive;
extern Function octave_sampled2continuous;
extern Function octave_schtrig;
extern Function octave_upsamplefill;
extern Function octave_wkeep;
extern Function octave_wrev;
extern Function octave_zerocrossing;
extern Function octave_fht;
extern Function octave_fwht;
extern Function octave_hilbert;
extern Function octave_idct;
extern Function octave_idct2;
extern Function octave_max;
extern Function octave_mean;
extern Function octave_meansq;
extern Function octave_median;
extern Function octave_min;
extern Function octave_plot;
extern Function octave_pause;
extern Function octave_abs;
extern Function octave_accumarray;
extern Function octave_accumdim;
extern Function octave_acos;
extern Function octave_acosd;
extern Function octave_acosh;
extern Function octave_acot;
extern Function octave_acotd;
extern Function octave_acoth;
extern Function octave_acsc;
extern Function octave_acsch;
extern Function octave_acscd;
extern Function octave_airy;
extern Function octave_adjoint;
extern Function octave_all;
extern Function octave_allow_non_integer_range_as_index;
extern Function octave_amd;
extern Function octave_ancestor;
extern Function octave_and;
extern Function octave_angle;
extern Function octave_annotation;
extern Function octave_anova;
extern Function octave_ans;
extern Function octave_any;
extern Function octave_arch_fit;
extern Function octave_arch_rnd;
extern Function octave_arch_test;
extern Function octave_area;
extern Function octave_arg;
extern Function octave_arrayfun;
extern Function octave_asec;
extern Function octave_asecd;
extern Function octave_asech;
extern Function octave_asin;
extern Function octave_asind;
extern Function octave_asinh;
extern Function octave_assume;
extern Function octave_assumptions;
extern Function octave_atan;
extern Function octave_atand;
extern Function octave_atanh;
extern Function octave_atan2;
extern Function octave_audiodevinfo;
extern Function octave_audioformats;
extern Function octave_audioinfo;
extern Function octave_audioread;
extern Function octave_audiowrite;
extern Function octave_autoreg_matrix;
extern Function octave_autumn;
extern Function octave_axes;
extern Function octave_axis;
extern Function octave_balance;
extern Function octave_bandwidth;
extern Function octave_bar;
extern Function octave_barh;
extern Function octave_bathannwin;
extern Function octave_bartlett_test;
extern Function octave_base2dec;
extern Function octave_base64_decode;
extern Function octave_base64_encode;
extern Function octave_beep;
extern Function octave_beep_on_error;
extern Function octave_bernoulli;
extern Function octave_besseli;
extern Function octave_besseljn;
extern Function octave_besselk;
extern Function octave_bessely;
extern Function octave_beta;
extern Function octave_betacdf;
extern Function octave_betainc;
extern Function octave_betaincinv;
extern Function octave_betainv;
extern Function octave_betain;
extern Function octave_betapdf;
extern Function octave_betarnd;
extern Function octave_bicg;
extern Function octave_bicgstab;
extern Function octave_bin2dec;
extern Function octave_bincoeff;
extern Function octave_binocdf;
extern Function octave_binoinv;
extern Function octave_binopdf;
extern Function octave_binornd;
extern Function octave_bitand;
extern Function octave_bitcmp;
extern Function octave_bitget;
extern Function octave_bitor;
extern Function octave_bitpack;
extern Function octave_bitset;
extern Function octave_bitshift;
extern Function octave_bitunpack;
extern Function octave_bitxor;
extern Function octave_blanks;
extern Function octave_blkdiag;
extern Function octave_blkmm;
extern Function octave_bone;
extern Function octave_box;
extern Function octave_brighten;
extern Function octave_bsxfun;
extern Function octave_builtin;
extern Function octave_bzip2;
extern Function octave_calendar;
extern Function octave_camlight;
extern Function octave_cart2pol;
extern Function octave_cart2sph;
extern Function octave_cast;
extern Function octave_cat;
extern Function octave_catalan;
extern Function octave_cauchy;
extern Function octave_cauchy_cdf;
extern Function octave_cauchy_inv;
extern Function octave_cauchy_pdf;
extern Function octave_cauchy_rnd;
extern Function octave_caxis;
extern Function octave_cbrt;
extern Function octave_ccode;
extern Function octave_ccolamd;
extern Function octave_ceil;
extern Function octave_center;
extern Function octave_centroid;
extern Function octave_cgs;
extern Function octave_chi2cdf;
extern Function octave_chi2inv;
extern Function octave_chi2pdf;
extern Function octave_chi2rnd;
extern Function octave_children;
extern Function octave_chisquare_test_homogeneity;
extern Function octave_chebyshevpoly;
extern Function octave_chebyshevT;
extern Function octave_chebyshevU;
extern Function octave_chol;
extern Function octave_chol2inv;
extern Function octave_choldelete;
extern Function octave_cholinsert;
extern Function octave_colinv;
extern Function octave_cholshift;
extern Function octave_cholupdate;
extern Function octave_chop;
extern Function octave_circshift;
extern Function octave_cla;
extern Function octave_clabel;
extern Function octave_clc;
extern Function octave_clf;
extern Function octave_clock;
extern Function octave_cloglog;
extern Function octave_cmpermute;
extern Function octave_cmunique;
extern Function octave_coeffs;
extern Function octave_colamd;
extern Function octave_colloc;
extern Function octave_colon;
extern Function octave_colorbar;
extern Function octave_colorcube;
extern Function octave_colormap;
extern Function octave_colperm;
extern Function octave_columns;
extern Function octave_comet;
extern Function octave_compan;
extern Function octave_compass;
extern Function octave_complex;
extern Function octave_computer;
extern Function octave_cond;
extern Function octave_condeig;
extern Function octave_condest;
extern Function octave_conj;
extern Function octave_contour;
extern Function octave_contour3;
extern Function octave_contourc;
extern Function octave_contourf;
extern Function octave_contrast;
extern Function octave_conv;
extern Function octave_conv2;
extern Function octave_convhull;
extern Function octave_convhulln;
extern Function octave_cool;
extern Function octave_copper;
extern Function octave_copyfile;
extern Function octave_copyobj;
extern Function octave_cor_test;
extern Function octave_cos;
extern Function octave_cosd;
extern Function octave_cosh;
extern Function octave_coshint;
extern Function octave_cosint;
extern Function octave_cot;
extern Function octave_cotd;
extern Function octave_coth;
extern Function octave_cov;
extern Function octave_cplxpair;
extern Function octave_cputime;
extern Function octave_cross;
extern Function octave_csc;
extern Function octave_cscd;
extern Function octave_csch;
extern Function octave_cstrcat;
extern Function octave_cstrcmp;
extern Function octave_csvread;
extern Function octave_csvwrite;
extern Function octave_csymamd;
extern Function octave_ctime;
extern Function octave_ctranspose;
extern Function octave_cubehelix;
extern Function octave_cummax;
extern Function octave_cummin;
extern Function octave_cumprod;
extern Function octave_cumsum;
extern Function octave_cumtrapz;
extern Function octave_cylinder;
extern Function octave_daspect;
extern Function octave_daspk;
extern Function octave_dasrt_options;
extern Function octave_dassl;
extern Function octave_dassl_options;
extern Function octave_date;
extern Function octave_datenum;
extern Function octave_datestr;
extern Function octave_datetick;
extern Function octave_dawson;
extern Function octave_dbclear;
extern Function octave_dbcont;
extern Function octave_dbdown;
extern Function octave_dblist;
extern Function octave_dblquad;
extern Function octave_dbquit;
extern Function octave_dbstack;
extern Function octave_dbstatus;
extern Function octave_dbstep;
extern Function octave_dbstop;
extern Function octave_dbtype;
extern Function octave_dbup;
extern Function octave_dbwhere;
extern Function octave_deal;
extern Function octave_deblank;
extern Function octave_dec2base;
extern Function octave_dec2hex;
extern Function octave_deconv;
extern Function octave_deg2rad;
extern Function octave_del2;
extern Function octave_delaunay;
extern Function octave_delaunayn;
extern Function octave_det;
extern Function octave_detrend;
extern Function octave_diag;
extern Function octave_diff;
extern Function octave_diffpara;
extern Function octave_diffuse;
extern Function octave_digits;
extern Function octave_dilog;
extern Function octave_dir;
extern Function octave_dirac;
extern Function octave_discrete_cdf;
extern Function octave_discrete_inv;
extern Function octave_discrete_pdf;
extern Function octave_discrete_rnd;
extern Function octave_disp;
extern Function octave_display;
extern Function octave_divergence;
extern Function octave_dimread;
extern Function octave_dimwrite;
extern Function octave_dmperm;
extern Function octave_do_string_escapes;
extern Function octave_doc;
extern Function octave_dot;
extern Function octave_double;
extern Function octave_downsample;
extern Function octave_dsearch;
extern Function octave_dsearchn;
extern Function octave_dsolve;
extern Function octave_dup2;
extern Function octave_duplication_matrix;
extern Function octave_durblevinson;
extern Function octave_e;
extern Function octave_ei;
extern Function octave_eig;
extern Function octave_ellipke;
extern Function octave_ellipsoid;
extern Function octave_ellipticCE;
extern Function octave_ellipticCK;
extern Function octave_ellipticCPi;
extern Function octave_ellipticE;
extern Function octave_ellipticF;
extern Function octave_ellipticK;
extern Function octave_ellipticPi;
extern Function octave_empirical_cdf;
extern Function octave_empirical_inv;
extern Function octave_empirical_pdf;
extern Function octave_empirical_rnd;
extern Function octave_end;
extern Function octave_endgrent;
extern Function octave_endpwent;
extern Function octave_eomday;
extern Function octave_eps;
extern Function octave_eq;
extern Function octave_equationsToMatrix;
extern Function octave_erf;
extern Function octave_erfc;
extern Function octave_erfinv;
extern Function octave_erfi;
extern Function octave_errno;
extern Function octave_error;
extern Function octave_error_ids;
extern Function octave_errorbar;
extern Function octave_etime;
extern Function octave_etree;
extern Function octave_etreeplot;
extern Function octave_eulier;
extern Function octave_eulergamma;
extern Function octave_evalin;
extern Function octave_exp;
extern Function octave_expand;
extern Function octave_expcdf;
extern Function octave_expint;
extern Function octave_expinv;
extern Function octave_expm;
extern Function octave_expm1;
extern Function octave_exppdf;
extern Function octave_exprnd;
extern Function octave_eye;
extern Function octave_ezcontour;
extern Function octave_ezcontourf;
extern Function octave_ezmesh;
extern Function octave_explot;
extern Function octave_ezplot3;
extern Function octave_ezsurf;
extern Function octave_ezpolar;
extern Function octave_ezsurfc;
extern Function octave_f_test_regression;
extern Function octave_factor;
extern Function octave_factorial;
extern Function octave_false;
extern Function octave_fcdf;
extern Function octave_fclear;
extern Function octave_fcntl;
extern Function octave_fdisp;
extern Function octave_feather;
extern Function octave_ff2n;
extern Function octave_fibonacci;
extern Function octave_find;
extern Function octave_findsym;
extern Function octave_finiteset;
extern Function octave_finv;
extern Function octave_fix;
extern Function octave_flintmax;
extern Function octave_flip;
extern Function octave_flipir;
extern Function octave_flipud;
extern Function octave_floor;
extern Function octave_fminbnd;
extern Function octave_fminunc;
extern Function octave_formula;
extern Function octave_fortran;
extern Function octave_fourier;
extern Function octave_fpdf;
extern Function octave_fplot;
extern Function octave_frac;
extern Function octave_fractdiff;
extern Function octave_frame2im;
extern Function octave_freport;
extern Function octave_fresneic;
extern Function octave_frnd;
extern Function octave_fskipl;
extern Function octave_fsolve;
extern Function octave_full;
extern Function octave_fwhm;
extern Function octave_fzero;
extern Function octave_gallery;
extern Function octave_gamcdf;
extern Function octave_gaminv;
extern Function octave_gamma;
extern Function octave_gammainc;
extern Function octave_gammaln;
extern Function octave_gca;
extern Function octave_gcbf;
extern Function octave_gcbo;
extern Function octave_gcd;
extern Function octave_ge;
extern Function octave_geocdf;
extern Function octave_geoinv;
extern Function octave_geopdf;
extern Function octave_geornd;
extern Function octave_givens;
extern Function octave_glpk;
extern Function octave_gmres;
extern Function octave_gmtime;
extern Function octave_gnplot_binary;
extern Function octave_gplot;
extern Function octave_gradient;
extern Function octave_gray;
extern Function octave_gray2ind;
extern Function octave_gt;
extern Function octave_gunzip;
extern Function octave_gzip;
extern Function octave_hadamard;
extern Function octave_hankel;
extern Function octave_harmonic;
extern Function octave_has;
extern Function octave_hash;
extern Function octave_heaviside;
extern Function octave_help;
extern Function octave_hess;
extern Function octave_hex2dec;
extern Function octave_hex2num;
extern Function octave_hilb;
extern Function octave_hilbert_curve;
extern Function octave_hist;
extern Function octave_horner;
extern Function octave_horzcat;
extern Function octave_hot;
extern Function octave_housh;
extern Function octave_hsv2rgb;
extern Function octave_hurst;
extern Function octave_hygecdf;
extern Function octave_hygeinv;
extern Function octave_hygepdf;
extern Function octave_hygernd;
extern Function octave_hypergeom;
extern Function octave_hypot;
extern Function octave_I;
extern Function octave_ichol;
extern Function octave_idist;
extern Function octave_idivide;
extern Function octave_igamma;
extern Function octave_ilaplace;
extern Function octave_ilu;
extern Function octave_im2double;
extern Function octave_im2frame;
extern Function octave_im2int16;
extern Function octave_im2single;
extern Function octave_im2uint16;
extern Function octave_im2uint8;
extern Function octave_imag;
extern Function octave_image;
extern Function octave_imagesc;
extern Function octave_imfinfo;
extern Function octave_imformats;
extern Function octave_importdata;
extern Function octave_imread;
extern Function octave_imshow;
extern Function octave_imwrite;
extern Function octave_ind2gray;
extern Function octave_ind2rgb;
extern Function octave_int2sub;
extern Function octave_index;
extern Function octave_inf;
extern Function octave_inpolygon;
extern Function octave_input;
extern Function octave_interp1;
extern Function octave_interp2;
extern Function octave_interp3;
extern Function octave_intersect;
extern Function octave_intmin;
extern Function octave_inv;
extern Function octave_invhilb;
extern Function octave_inimpinvar;
extern Function octave_ipermute;
extern Function octave_iqr;
extern Function octave_isa;
extern Function octave_isequal;
extern Function octave_ishermitian;
extern Function octave_isprime;
extern Function octave_jit_enable;
extern Function octave_kbhit;
extern Function octave_kendall;
extern Function octave_kron;
extern Function octave_kurtosis;
extern Function octave_laplace;
extern Function octave_laplace_cdf;
extern Function octave_laplace_inv;
extern Function octave_laplace_pdf;
extern Function octave_laplace_rnd;
extern Function octave_laplacian;
extern Function octave_lcm;
extern Function octave_ldivide;
extern Function octave_le;
extern Function octave_legendre;
extern Function octave_length;
extern Function octave_lgamma;
extern Function octave_limit;
extern Function octave_line;
extern Function octave_linprog;
extern Function octave_linsolve;
extern Function octave_linspace;
extern Function octave_load;
extern Function octave_log;
extern Function octave_log10;
extern Function octave_log1p;
extern Function octave_log2;
extern Function octave_logical;
extern Function octave_logistic_cdf;
extern Function octave_logistic_inv;
extern Function octave_logistic_pdf;
extern Function octave_logistic_regression;
extern Function octave_logit;
extern Function octave_loglog;
extern Function octave_loglogerr;
extern Function octave_logm;
extern Function octave_logncdf;
extern Function octave_logninv;
extern Function octave_lognpdf;
extern Function octave_lognrnd;
extern Function octave_lognspace;
extern Function octave_lookup;
extern Function octave_lscov;
extern Function octave_lsode;
extern Function octave_lsqnonneg;
extern Function octave_lt;
extern Function octave_magic;
extern Function octave_manova;
extern Function octave_minus;
extern Function octave_mkpp;
extern Function octave_mldivide;
extern Function octave_mod;
extern Function octave_moment;
extern Function octave_mpoles;
extern Function octave_mpower;
extern Function octave_mrdivide;
extern Function octave_mu2lin;
extern Function octave_na;
extern Function octave_nan;
extern Function octave_nextpow2;
extern Function octave_nnz;
extern Function octave_nonzeros;
extern Function octave_norm;
extern Function octave_normcdf;
extern Function octave_normest;
extern Function octave_normest1;
extern Function octave_norminv;
extern Function octave_normpdf;
extern Function octave_normrnd;
extern Function octave_nth_element;
extern Function octave_nth_root;
extern Function octave_null;
extern Function octave_numel;
extern Function octave_ode23;
extern Function octave_ode45;
extern Function octave_ols;
extern Function octave_ones;
extern Function octave_prod;
extern Function octave_power;
extern Function octave_sin;
extern Function octave_sqrt;
extern Function octave_sum;
extern Function octave_sumsq;
extern Function octave_tan;
extern Function octave_tanh;
extern Function octave_sinh;
extern Function octave_bin_values;
extern Function octave_catmullrom;
extern Function octave_csape;
extern Function octave_csapi;
extern Function octave_csaps;
extern Function octave_csaps_sel;
extern Function octave_dedup;
extern Function octave_fnder;
extern Function octave_fnplt;
extern Function octave_fnval;
extern Function octave_regularization;
extern Function octave_regularization2D;
extern Function octave_tpaps;
extern Function octave_tps_val;
extern Function octave_tps_val_der;
//extern Function octave_rms;
extern Function octave_normalize;
extern Function octave_gaindb;
extern Function octave_crestfactor;
extern Function octave_uquant;
extern Function octave_firwin;
extern Function octave_firkaiser;
extern Function octave_fir2long;
extern Function octave_long2fir;
extern Function octave_freqwin;
extern Function octave_firfilter;
extern Function octave_blfilter;
extern Function octave_warpedblfilter;
extern Function octave_freqfilter;
extern Function octave_pfilt;
extern Function octave_magresp;
extern Function octave_transferfunction;
extern Function octave_pgrdelay;
extern Function octave_rampup;
extern Function octave_rampdown;
extern Function octave_thresh;
extern Function octave_largestr;
extern Function octave_largestn;
extern Function octave_dynlimit;
extern Function octave_groupthresh;
extern Function octave_rgb2jpeg;
extern Function octave_jpeg2rgb;
extern Function octave_qam4;
extern Function octave_iqam4;
extern Function octave_semiaudplot;
extern Function octave_audtofreq;
extern Function octave_freqtoaud;
extern Function octave_audspace;
extern Function octave_audspacebw;
extern Function octave_erbtofreq;
extern Function octave_freqtoerb;
extern Function octave_erbspace;

}

