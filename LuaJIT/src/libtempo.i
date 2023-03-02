%module libtempo
%{

#include "curve_utils.h"
#include "defines.h"
#include "fourier_utils.h"
#include "math_utils.hpp"
#include "mat_utils.h"
#include "signal_utils.h"
#include "tempogram_utils.h"
#include "resample.h"

#include <armadillo>

#include "tempogram_processing.h"

using namespace libtempo;
using namespace libtempo::curve_utils;
using namespace libtempo::fourier_utils;
using namespace libtempo::mat_utils;
using namespace libtempo::utils::math;
using namespace libtempo::signal_utils;
using namespace libtempo::tempogram_utils;
using namespace libtempo::tempogram_processing;
#include <fftw/fftw.h>
using namespace fftw;
%}

namespace libtempo::curve_utils {

    std::vector<uvec> split_curve(const vec &curve);
    std::vector<uvec> join_adjacent_segments(const std::vector<uvec> &segments);
    vec correct_curve_by_length(const vec &measurements, int min_length);
    vec correct_curve_by_confidence(const vec &measurements, const vec &confidence, float threshold = 0.85);

    struct Section 
    {
        Section();
        Section(double start, double end, double bpm, double offset = 0);

        double start, end, bpm, offset;
    };

    //std::ostream &operator<<(std::ostream &os, const Section &section);
    std::vector<Section> tempo_segments_to_sections(const std::vector<uvec> &segments, const vec &curve, const vec &t,
                                                    double bpm_reference = DEFAULT_REF_TEMPO);
    void split_section(const Section &section, std::vector<Section> &sections, double max_section_len = 60);
    std::vector<Section> merge_sections(const std::vector<Section> &sections, float threshold = 0.5f);
    Section average_sections(const std::vector<Section> &sections);
    void extract_offset(const vec &novelty_curve, Section &section, const std::vector<int> &tempo_multiples,
                        int feature_rate, float bpm_doubt_window = 2, double bpm_doubt_step = 0.1);
    void correct_offset(Section &section, int smallest_fraction_shift = 4);
};

namespace libtempo::fourier_utils {

    cx_mat
    compute_fourier_coefficients(vec &t, const vec &s, const vec &window, int n_overlap, const vec &f, double sr);

    template<typename T>
    mat stft(float &feature_rate, vec &t, vec &f, const Col<T> &signal, int sr, const vec &window,
             std::tuple<int, int> coefficient_range, int n_fft = -1, int hop_length = -1);
    
    template<typename T>
    mat stft(float &feature_rate, vec &t, vec &f, const Col<T> &s, int sr, const vec &window, int n_fft = -1,
             int hop_length = -1);
}


namespace libtempo::mat_utils {

    uword find_nearest(const vec &a, double value);
    mat subtract_mean(const mat &tempogram);
    cx_mat colwise_normalize_p1(const cx_mat &feature, unsigned int p, double threshold);
    mat rowwise_normalize(const mat &tempogram);
    uvec argmax(const mat &tempogram);
    vec max_bucket(const mat &tempogram, const vec &axis_lut);
    vec pad_vec(const vec &x, unsigned int n_before, unsigned int n_after);
    mat pad_mat(const mat &data, unsigned int n_before, unsigned int n_after, bool repeat = true);
    vec smooth_filter_subtract(const vec &novelty_curve, int sr, int hop_length);


    template<typename T>
    void mat_interp1_nearest(const Col<T>& XG, const Mat<T>& YG, const Col<T>& XI, Mat<T>& YI);
};

namespace libtempo::utils::math {

    int fix(double x);
    vec sinc_fac(const vec &x);

    arma::vec my_hanning( const arma::uword N );
    int calc_gcd(int n1, int n2);

    int quotient_ceil(int n1, int n2);
};


namespace libtempo {
    
    class Resampler {
    public:
        Resampler(int up_rate, int down_rate, const vec &coefs);
        virtual ~Resampler();
        int out_count(int in_count);
        vec apply(vec &in);
        int coefs_per_phase();
    };

    vec upfirdn(const vec &signal, int up_rate, int down_rate, const vec &filter);
    vec firls(int length, vec freq, const vec &amplitude);
    vec resample(const vec &signal, int upfactor, int downfactor);
}

namespace libtempo::signal_utils {
    std::tuple<vec, vec> generate_pulse(double bpm, int window_length, int feature_rate, int shift_by=0);
}

namespace libtempo::tempogram_utils {
    vec extract_confidence(const mat &tempogram, const bool &rollover = true);
    mat tempogram_apply_window(const mat &tempogram, const vec &window);
    mat tempogram_include_triplets(const mat &tempogram, const vec &axis_lut, float weight = 2.f);
    mat smoothen_tempogram(const mat &tempogram, const vec &axis_lut, int temporal_unit_size = 100,
                           float triplet_weight = 0.8f);
    vec extract_tempo_curve(const mat &tempogram, const vec &axis_lut);
    double bpm_to_cyclic(double bpm, double ref_tempo = 60);
}



namespace libtempo::tempogram_processing
{
    vec audio_to_novelty_curve(int &feature_rate_ret, const fvec &signal, int sr, int window_length = -1, int hop_length = -1,
                                                double compression_c = 1000, bool log_compression = true,
                                                int resample_feature_rate = 200);


    cx_mat novelty_curve_to_tempogram_dft(vec &t, const vec &novelty_curve, const vec &bpms, double feature_rate,
                                                                int tempo_window, int hop_length = -1);

    mat tempogram_to_cyclic_tempogram(vec &y_axis, const cx_mat &tempogram, const vec &bpm, int octave_divider = 120,
                                                       int ref_tempo = DEFAULT_REF_TEMPO);
}

%template(stft_float) libtempo::fourier_utils::stft<float>;
%template(mat_interp1_nearest_float) libtempo::mat_utils::mat_interp1_nearest<float>;
