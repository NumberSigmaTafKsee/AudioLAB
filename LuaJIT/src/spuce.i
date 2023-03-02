%module spuce
%{

#include "filters/biquad.h"
#include "filters/fir.h"
#include "filters/iir.h"

//#include "filters/fir_adapt.h"
#include "filters/fir_coeff.h"
#include "filters/fir_decim.h"
#include "filters/fir_interp.h"
#include "filters/fir_inv_dft.h"

#include "filters/design_fir.h"
#include "filters/design_iir.h"
#include "filters/design_window.h"


#include "filters/iir_coeff.h"
#include "filters/iir_df.h"
#include "filters/iir_comb.h"
#include "filters/iir_shelf.h"
#include "filters/iir_1st.h"

//#include "filters/iir_allpass1_sections.h"
//#include "filters/iir_allpass1_sections_variable_delay.h"
//#include "filters/iir_allpass_variable_cascade.h"


#include "filters/allpass_1.h"
#include "filters/allpass_2nd.h"
#include "filters/allpass.h"

#include "filters/audio_equalizer.h"

#include "filters/butterworth_allpass.h"
#include "filters/butterworth_iir.h"
#include "filters/butterworth_fir.h"

#include "filters/calculate_decimator_taps.h"
#include "filters/cascaded_cic.h"
#include "filters/create_remez_lpfir.h"
#include "filters/cutboost.h"
#include "filters/decimator.h"

#include "filters/elliptic_allpass.h"
#include "filters/elliptic_iir.h"
#include "filters/farrow.h"
#include "filters/farrow_upsampler.h"
#include "filters/find_roots.h"

#include "filters/gaussian_fir.h"

#include "filters/lagrange.h"

#include "filters/notch_allpass.h"
#include "filters/notch_comb.h"
#include "filters/notch_iir.h"

#include "filters/raised_cosine.h"
#include "filters/raised_cosine_imp.h"

#include "filters/remez_estimate.h"
#include "filters/remez_fir.h"

#include "filters/root_raised_cosine.h"
#include "filters/root_raised_cosine_imp.h"

#include "filters/running_average.h"
#include "filters/running_sum.h"

#include "filters/scic.h"

#include "filters/shelf_allpass1.h"
#include "filters/sinc_fir.h"
#include "filters/sinc_helper.h"
#include "filters/transform_fir.h"
#include "filters/window.h"

#include "filters/cic.h"
#include "filters/cascaded_cic.h"

#include "filters/farrow_upsampler.h"
#include "filters/lagrange.h"

using namespace spuce;

%}

%include "filters/biquad.h"
%include "filters/fir.h"
%include "filters/iir.h"

%include "filters/cic.h"
%include "filters/cascaded_cic.h"

%include "filters/iir_coeff.h"
%include "filters/iir_df.h"
%include "filters/iir_comb.h"
%include "filters/iir_shelf.h"
%include "filters/iir_1st.h"

//%include "filters/iir_allpass1_sections.h"
//%include "filters/iir_allpass1_sections_variable_delay.h"
//%include "filters/iir_allpass_variable_cascade.h"

//%include "filters/fir_adapt.h"
%include "filters/fir_coeff.h"
%include "filters/fir_decim.h"
%include "filters/fir_interp.h"
%include "filters/fir_inv_dft.h"


//%include "filters/allpass_1.h"
//%include "filters/allpass_2nd.h"
//%include "filters/allpass.h"

//%include "filters/audio_equalizer.h"

%include "filters/butterworth_allpass.h"
%include "filters/butterworth_iir.h"
%include "filters/butterworth_fir.h"

%include "filters/calculate_decimator_taps.h"
%include "filters/cascaded_cic.h"
%include "filters/create_remez_lpfir.h"
%include "filters/cutboost.h"

%include "filters/decimator.h"

%include "filters/design_fir.h"
%include "filters/design_iir.h"
%include "filters/design_window.h"

%include "filters/elliptic_allpass.h"
%include "filters/elliptic_iir.h"

%include "filters/farrow.h"
%include "filters/farrow_upsampler.h"
%include "filters/find_roots.h"
%include "filters/fir_adapt.h"
%include "filters/fir_coeff.h"
%include "filters/fir_decim.h"

%include "filters/gaussian_fir.h"

%include "filters/notch_allpass.h"
%include "filters/notch_comb.h"
%include "filters/notch_iir.h"

%include "filters/raised_cosine.h"
%include "filters/raised_cosine_imp.h"

%include "filters/remez_estimate.h"
%include "filters/remez_fir.h"

%include "filters/root_raised_cosine.h"
%include "filters/root_raised_cosine_imp.h"

%include "filters/running_average.h"
%include "filters/running_sum.h"

%include "filters/scic.h"
%include "filters/cic.h"

%include "filters/shelf_allpass1.h"
%include "filters/sinc_fir.h"
%include "filters/sinc_helper.h"
%include "filters/transform_fir.h"
%include "filters/window.h"
%include "filters/farrow_upsampler.h"
%include "filters/lagrange.h"

//%template(float_allpass) spuce::allpass<float>;
//%template(double_allpass) spuce::allpass<double>;

//%template(float_allpass_1) spuce::allpass_1<float>;
//%template(double_allpass_1) spuce::allpass_1<double>;

//%template(float_allpass_2) spuce::allpass_2nd<float>;
//%template(double_allpass_2) spuce::allpass_2nd<double>;

%template(float_biquad) spuce::biquad<float>;
%template(double_biquad) spuce::biquad<double>;

%template(float_cascaded_cic) spuce::cascaded_cic<float>;
%template(double_cascaded_cic) spuce::cascaded_cic<double>;

%template(float_cic) spuce::cic<float>;
%template(double_cic) spuce::cic<double>;

%template(float_cutboost) spuce::cutboost<float>;
%template(double_cutboost) spuce::cutboost<double>;

%template(float_decimator) spuce::decimator<float>;
%template(double_decimator) spuce::decimator<double>;

%template(float_farrow) spuce::farrow<float>;
%template(double_farrow) spuce::farrow<double>;

%template(float_fir) spuce::fir<float>;
%template(double_fir) spuce::fir<double>;

//%template(float_fir_adapt) spuce::fir_adapt<float>;
//%template(double_fir_adapt) spuce::fir_adapt<double>;

%template(float_fir_coeff) spuce::fir_coeff<float>;
%template(double_fir_coeff) spuce::fir_coeff<double>;

%template(float_fir_decim) spuce::fir_decim<float>;
%template(double_fir_decim) spuce::fir_decim<double>;

%template(float_fir_interp) spuce::fir_interp<float>;
%template(double_fir_interp) spuce::fir_interp<double>;

%template(float_iir) spuce::iir<float>;
%template(double_iir) spuce::iir<double>;

%template(float_iir_1) spuce::iir_1st<float>;
%template(double_iir_1) spuce::iir_1st<double>;

//%template(float_iir_allpass1_sections) spuce::iir_allpass1_sections<float>;
//%template(double_iir_allpass1_sections) spuce::iir_allpass1_sections<double>;

//%template(float_iir_allpass1_sections_variable_delay) spuce::iir_allpass1_sections_variable_delay<float>;
//%template(double_iir_allpass1_sections_variable_delay) spuce::iir_allpass1_sections_variable_delay<double>;

//%template(float_iir_allpass_variable_cascade) spuce::iir_allpass_variable_cascade<float>;
//%template(double_iir_allpass_variable_cascade) spuce::iir_allpass_variable_cascade<double>;


%template(float_iir_comb) spuce::iir_comb<float>;
%template(double_iir_comb) spuce::iir_comb<double>;

%template(float_iir_df) spuce::iir_df<float>;
%template(double_iir_df) spuce::iir_df<double>;

%template(float_iir_shelf) spuce::iir_shelf<float>;
%template(double_iir_shelf) spuce::iir_shelf<double>;

%template(float_lagrange) spuce::lagrange<float>;
%template(double_lagrange) spuce::lagrange<double>;

%template(float_notch_allpass) spuce::notch_allpass<float>;
%template(double_notch_allpass) spuce::notch_allpass<double>;

%template(float_notch_comb) spuce::notch_comb<float>;
%template(double_notch_comb) spuce::notch_comb<double>;

%template(float_notch_iir) spuce::notch_iir<float>;
%template(double_notch_iir) spuce::notch_iir<double>;

%template(float_raised_cosine) spuce::raised_cosine<float>;
%template(double_raised_cosine) spuce::raised_cosine<double>;

%template(float_root_raised_cosine_quantized) spuce::root_raised_cosine_quantized<float>;
%template(double_root_raised_cosine_quantized) spuce::root_raised_cosine_quantized<double>;

%template(float_running_average) spuce::running_average<float>;
%template(double_running_average) spuce::running_average<double>;

%template(float_running_sum) spuce::running_sum<float>;
%template(double_running_sum) spuce::running_sum<double>;

%template(float_scic) spuce::scic<float>;
%template(double_scic) spuce::scic<double>;

//%template(float_shelf_allpass1) spuce::shelf_allpass1<float>;
//%template(double_shelf_allpass1) spuce::shelf_allpass1<double>;





