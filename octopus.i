%module octopus
%{
#include "octopus.hpp"
//#include "octopus_functions.hpp"
//#include "audiodsp_std_samples.hpp"
//#include "octopus_octavate.hpp"
using namespace Octopus;
%}


%include "std_math.i"
%include "std_string.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

// Get the STL typemaps
//%include "stl.i"

// Handle standard exceptions
%include "std_exception.i"

//%include "octave.i"
//%include "octopus_rowvector.hpp"
%include "octopus_rowvectorxf.hpp"
%include "octopus_rowvectorxd.hpp"
%include "octopus_rowvectorxcf.hpp"
%include "octopus_rowvectorxcd.hpp"
//%include "octopus_colvector.hpp"
%include "octopus_colvectorxf.hpp"
%include "octopus_colvectorxd.hpp"
%include "octopus_colvectorxcf.hpp"
%include "octopus_colvectorxcd.hpp"
//%include "octopus_matrix.hpp"
%include "octopus_matrixxf.hpp"
%include "octopus_matrixxd.hpp"
%include "octopus_matrixxcf.hpp"
%include "octopus_matrixxcd.hpp"
//%include "octopus_octavate.hpp"
//%include "octopus_functions.hpp"
%include "octopus.hpp"

%inline
%{
namespace Octopus
{
  OctopusValueList evalFunc(const std::string& func, const OctopusValueList& inputs, int noutputs=1)
  {          
      octave_value_list out = octave::feval(func.c_str(), inputs.vlist, noutputs);
      return OctopusValueList(out);
  }

  OctopusValueList Function::eval(const OctopusValueList &input, int numOutputs)
  {
      return evalFunc(name,input,numOutputs);
  }
    
  Function octave_fft= Function("fft");
  Function octave_ifft= Function("ifft");
  Function octave_fft2= Function("fft2");
  Function octave_ifft2= Function("ifft2");
  Function octave_fftconv= Function("fftconv");
  Function octave_fftfilt= Function("fftfilt");
  Function octave_fftn= Function("fftn");
  Function octave_fftshift= Function("fftshift");
  Function octave_fftw= Function("fftw");
  Function octave_ifftn= Function("ifftn");
  Function octave_ifftshift= Function("ifftshift");
  Function octave_ifht= Function("ifht");
  Function octave_ifourier= Function("ifourier");
  Function octave_ifwht= Function("ifwht");
  Function octave_ifwt= Function("ifwt");
  Function octave_ifwt2= Function("ifwt2");
  Function octave_buffer= Function("buffer");
  Function octave_chirp= Function("chirp");
  Function octave_cmorwavf= Function("cmorwavf");
  Function octave_gauspuls= Function("gauspuls");
  Function octave_gmonopuls= Function("gmonopuls");
  Function octave_mexihat= Function("mexihat");
  Function octave_meyeraux= Function("meyeraux");
  Function octave_morlet= Function("morlet");
  Function octave_pulstran= Function("pulstran");
  Function octave_rectpuls= Function("rectpuls");
  Function octave_sawtooth= Function("sawtooth");
  Function octave_shanwavf= Function("shanwavf");
  Function octave_shiftdata= Function("shiftdata");
  Function octave_sigmoid_train= Function("sigmoid_train");
  Function octave_specgram= Function("specgram");
  Function octave_square= Function("square");
  Function octave_tripuls= Function("tripuls");
  Function octave_udecode= Function("udecode");
  Function octave_uencoder= Function("uencoder");
  Function octave_unshiftdata= Function("unshiftdata");
  Function octave_findpeaks= Function("findpeaks");
  Function octave_peak2peak= Function("peak2peak");
  Function octave_peak2rms= Function("peak2rms");
  Function octave_rms= Function("rms");
  Function octave_rssq= Function("rssq");
  Function octave_cconv= Function("cconv");
  Function octave_convmtx= Function("convmtx");
  Function octave_wconv= Function("wconv");
  Function octave_xcorr= Function("xcorr");
  Function octave_xcorr2= Function("xcorr2");
  Function octave_xcov= Function("xcov");
  Function octave_filtfilt= Function("filtfilt");
  Function octave_fltic= Function("fltic");
  Function octave_medfilt1= Function("medfilt1");
  Function octave_movingrms= Function("movingrms");
  Function octave_sgolayfilt= Function("sgolayfilt");
  Function octave_sosfilt= Function("sosfilt");
  Function octave_freqs= Function("freqs");
  Function octave_freqs_plot= Function("freqs_plot");
  Function octave_freqz= Function("freqz");
  Function octave_freqz_plot= Function("freqz_plot");
  Function octave_impz= Function("impz");
  Function octave_zplane= Function("zplane");
  Function octave_filter= Function("filter");
  Function octave_filter2= Function("filter2");
  Function octave_fir1= Function("fir1");
  Function octave_fir2= Function("fir2");
  Function octave_firls= Function("firls");
  Function octave_sinc= Function("sinc");
  Function octave_unwrap= Function("unwrap");
  Function octave_bartlett= Function("bartlett");
  Function octave_blackman= Function("blackman");
  Function octave_blackmanharris= Function("blackmanharris");
  Function octave_blackmannuttal= Function("blackmannuttal");
  Function octave_dftmtx= Function("dftmtx");
  Function octave_hamming= Function("hamming");
  Function octave_hann= Function("hann");
  Function octave_hanning= Function("hanning");
  Function octave_pchip= Function("pchip");
  Function octave_periodogram= Function("periodogram");
  Function octave_sinetone= Function("sinetone");
  Function octave_sinewave= Function("sinewave");
  Function octave_spectral_adf= Function("spectral_adf");
  Function octave_spectral_xdf= Function("spectral_xdf");
  Function octave_spencer= Function("spencer");
  Function octave_stft= Function("stft");
  Function octave_synthesis= Function("synthesis");
  Function octave_yulewalker= Function("yulewalker");
  Function octave_polystab= Function("polystab");
  Function octave_residued= Function("residued");
  Function octave_residuez= Function("residuez");
  Function octave_sos2ss= Function("sos2ss");
  Function octave_sos2tf= Function("sos2tf");
  Function octave_sos2zp= Function("sos2zp");
  Function octave_ss2tf= Function("ss2tf");
  Function octave_ss2zp= Function("ss2zp");
  Function octave_tf2sos= Function("tf2sos");
  Function octave_tf2ss= Function("tf2ss");
  Function octave_tf2zp= Function("tf2zp");
  Function octave_zp2sos= Function("zp2sos");
  Function octave_zp2ss= Function("zp2ss");
  Function octave_zp2tf= Function("zp2tf");
  Function octave_besselap= Function("besselap");
  Function octave_besself= Function("besself");
  Function octave_bilinear= Function("bilinear");
  Function octave_buttap= Function("buttap");
  Function octave_butter= Function("butter");
  Function octave_buttord= Function("buttord");
  Function octave_cheb= Function("cheb");
  Function octave_cheb1ap= Function("cheb1ap");
  Function octave_cheb1ord= Function("cheb1ord");
  Function octave_cheb2ap= Function("cheb2ap");
  Function octave_cheb2ord= Function("cheb2ord");
  Function octave_chebywin= Function("chebywin");
  Function octave_cheby1= Function("cheby1");
  Function octave_cheby2= Function("cheby2");
  Function octave_ellip= Function("ellip");
  Function octave_ellipap= Function("ellipap");
  Function octave_ellipord= Function("ellipord");
  Function octave_impinvar= Function("impinvar");
  Function octave_ncauer= Function("ncauer");
  Function octave_pei_tseng_notch= Function("pei_tseng_notch");
  Function octave_sftrans= Function("sftrans");
  Function octave_cl2bp= Function("cl2bp");
  Function octave_kaiserord= Function("kaiserord");
  Function octave_qp_kaiser= Function("qp_kaiser");
  Function octave_remez= Function("remez");
  Function octave_sgplay= Function("sgplay");
  Function octave_bitrevorder= Function("bitrevorder");
  Function octave_cceps= Function("cceps");
  Function octave_cplxreal= Function("cplxreal");
  Function octave_czt= Function("czt");
  Function octave_dct= Function("dct");
  Function octave_dct2= Function("dct2");
  Function octave_dctmtx= Function("dctmtx");
  Function octave_digitrevorder= Function("digitrevorder");
  Function octave_dst= Function("dst");
  Function octave_dwt= Function("dwt");
  Function octave_rceps= Function("rceps");
  Function octave_ar_psd= Function("ar_psd");
  Function octave_cohere= Function("cohere");
  Function octave_cpsd= Function("cpsd");
  Function octave_csd= Function("csd");
  Function octave_db2pow= Function("db2pow");
  Function octave_mscohere= Function("mscohere");
  Function octave_pburg= Function("pburg");
  Function octave_pow2db= Function("pow2db");
  Function octave_pwelch= Function("pwelch");
  Function octave_pyulear= Function("pyulear");
  Function octave_tfe= Function("tfe");
  Function octave_tfestimate= Function("tfestimate");
  Function octave___power= Function("__power");
  Function octave_barthannwin= Function("barthannwin");
  Function octave_bohmanwin= Function("bohmanwin");
  Function octave_boxcar= Function("boxcar");
  Function octave_flattopwin= Function("flattopwin");
  Function octave_chebwin= Function("chebwin");
  Function octave_gaussian= Function("gaussian");
  Function octave_gausswin= Function("gausswin");
  Function octave_kaiser= Function("kaiser");
  Function octave_nuttalwin= Function("nuttalwin");
  Function octave_parzenwin= Function("parzenwin");
  Function octave_rectwin= Function("rectwin");
  Function octave_tukeywin= Function("tukeywin");
  Function octave_ultrwin= Function("ultrwin");
  Function octave_welchwin= Function("welchwin");
  Function octave_window= Function("window");
  Function octave_arburg= Function("arburg");
  Function octave_aryule= Function("aryule");
  Function octave_invfreq= Function("invfreq");
  Function octave_invfreqz= Function("invfreqz");
  Function octave_invfreqs= Function("invfreqs");
  Function octave_levinson= Function("levinson");
  Function octave_data2fun= Function("data2fun");
  Function octave_decimate= Function("decimate");
  Function octave_interp= Function("interp");
  Function octave_resample= Function("resample");
  Function octave_upfirdn= Function("upfirdn");
  Function octave_upsample= Function("upsample");
  Function octave_clustersegment= Function("clustersegment");
  Function octave_fracshift= Function("fracshift");
  Function octave_marcumq= Function("marcumq");
  Function octave_primitive= Function("primitive");
  Function octave_sampled2continuous= Function("sampled2continuous");
  Function octave_schtrig= Function("schtrig");
  Function octave_upsamplefill= Function("upsamplefill");
  Function octave_wkeep= Function("wkeep");
  Function octave_wrev= Function("wrev");
  Function octave_zerocrossing= Function("zerocrossing");
  Function octave_fht= Function("fht");
  Function octave_fwht= Function("fwht");
  Function octave_hilbert= Function("hilbert");
  Function octave_idct= Function("idct");
  Function octave_idct2= Function("idct2");
  Function octave_max= Function("max");
  Function octave_mean= Function("mean");
  Function octave_meansq= Function("meansq");
  Function octave_median= Function("median");
  Function octave_min= Function("min");
  Function octave_plot= Function("plot");
  Function octave_pause= Function("pause");
  Function octave_abs= Function("abs");
  Function octave_accumarray= Function("accumarray");
  Function octave_accumdim= Function("accumdim");
  Function octave_acos= Function("acos");
  Function octave_acosd= Function("acosd");
  Function octave_acosh= Function("acosh");
  Function octave_acot= Function("acot");
  Function octave_acotd= Function("acotd");
  Function octave_acoth= Function("acoth");
  Function octave_acsc= Function("acsc");
  Function octave_acsch= Function("acsch");
  Function octave_acscd= Function("acscd");
  Function octave_airy= Function("airy");
  Function octave_adjoint= Function("adjoint");
  Function octave_all= Function("all");
  Function octave_allow_non_integer_range_as_index= Function("allow_non_integer_range_as_index");
  Function octave_amd= Function("amd");
  Function octave_ancestor= Function("ancestor");
  Function octave_and= Function("and");
  Function octave_angle= Function("angle");
  Function octave_annotation= Function("annotation");
  Function octave_anova= Function("anova");
  Function octave_ans= Function("ans");
  Function octave_any= Function("any");
  Function octave_arch_fit= Function("arch_fit");
  Function octave_arch_rnd= Function("arch_rnd");
  Function octave_arch_test= Function("arch_test");
  Function octave_area= Function("area");
  Function octave_arg= Function("arg");
  Function octave_arrayfun= Function("arrayfun");
  Function octave_asec= Function("asec");
  Function octave_asecd= Function("asecd");
  Function octave_asech= Function("asech");
  Function octave_asin= Function("asin");
  Function octave_asind= Function("asind");
  Function octave_asinh= Function("asinh");
  Function octave_assume= Function("assume");
  Function octave_assumptions= Function("assumptions");
  Function octave_atan= Function("atan");
  Function octave_atand= Function("atand");
  Function octave_atanh= Function("atanh");
  Function octave_atan2= Function("atan2");
  Function octave_audiodevinfo= Function("audiodevinfo");
  Function octave_audioformats= Function("audioformats");
  Function octave_audioinfo= Function("audioinfo");
  Function octave_audioread= Function("audioread");
  Function octave_audiowrite= Function("audiowrite");
  Function octave_autoreg_matrix= Function("autoreg_matrix");
  Function octave_autumn= Function("autumn");
  Function octave_axes= Function("axes");
  Function octave_axis= Function("axis");
  Function octave_balance= Function("balance");
  Function octave_bandwidth= Function("bandwidth");
  Function octave_bar= Function("bar");
  Function octave_barh= Function("barh");
  Function octave_bathannwin= Function("bathannwin");
  Function octave_bartlett_test= Function("bartlett_test");
  Function octave_base2dec= Function("base2dec");
  Function octave_base64_decode= Function("base64_decode");
  Function octave_base64_encode= Function("base64_encode");
  Function octave_beep= Function("beep");
  Function octave_beep_on_error= Function("beep_on_error");
  Function octave_bernoulli= Function("bernoulli");
  Function octave_besseli= Function("besseli");
  Function octave_besseljn= Function("besseljn");
  Function octave_besselk= Function("besselk");
  Function octave_bessely= Function("bessely");
  Function octave_beta= Function("beta");
  Function octave_betacdf= Function("betacdf");
  Function octave_betainc= Function("betainc");
  Function octave_betaincinv= Function("betaincinv");
  Function octave_betainv= Function("betainv");
  Function octave_betain= Function("betain");
  Function octave_betapdf= Function("betapdf");
  Function octave_betarnd= Function("betarnd");
  Function octave_bicg= Function("bicg");
  Function octave_bicgstab= Function("bicgstab");
  Function octave_bin2dec= Function("bin2dec");
  Function octave_bincoeff= Function("bincoeff");
  Function octave_binocdf= Function("binocdf");
  Function octave_binoinv= Function("binoinv");
  Function octave_binopdf= Function("binopdf");
  Function octave_binornd= Function("binornd");
  Function octave_bitand= Function("bitand");
  Function octave_bitcmp= Function("bitcmp");
  Function octave_bitget= Function("bitget");
  Function octave_bitor= Function("bitor");
  Function octave_bitpack= Function("bitpack");
  Function octave_bitset= Function("bitset");
  Function octave_bitshift= Function("bitshift");
  Function octave_bitunpack= Function("bitunpack");
  Function octave_bitxor= Function("bitxor");
  Function octave_blanks= Function("blanks");
  Function octave_blkdiag= Function("blkdiag");
  Function octave_blkmm= Function("blkmm");
  Function octave_bone= Function("bone");
  Function octave_box= Function("box");
  Function octave_brighten= Function("brighten");
  Function octave_bsxfun= Function("bsxfun");
  Function octave_builtin= Function("builtin");
  Function octave_bzip2= Function("bzip2");
  Function octave_calendar= Function("calendar");
  Function octave_camlight= Function("camlight");
  Function octave_cart2pol= Function("cart2pol");
  Function octave_cart2sph= Function("cart2sph");
  Function octave_cast= Function("cast");
  Function octave_cat= Function("cat");
  Function octave_catalan= Function("catalan");
  Function octave_cauchy= Function("cauchy");
  Function octave_cauchy_cdf= Function("cauchy_cdf");
  Function octave_cauchy_inv= Function("cauchy_inv");
  Function octave_cauchy_pdf= Function("cauchy_pdf");
  Function octave_cauchy_rnd= Function("cauchy_rnd");
  Function octave_caxis= Function("caxis");
  Function octave_cbrt= Function("cbrt");
  Function octave_ccode= Function("ccode");
  Function octave_ccolamd= Function("ccolamd");
  Function octave_ceil= Function("ceil");
  Function octave_center= Function("center");
  Function octave_centroid= Function("centroid");
  Function octave_cgs= Function("cgs");
  Function octave_chi2cdf= Function("chi2cdf");
  Function octave_chi2inv= Function("chi2inv");
  Function octave_chi2pdf= Function("chi2pdf");
  Function octave_chi2rnd= Function("chi2rnd");
  Function octave_children= Function("children");
  Function octave_chisquare_test_homogeneity= Function("chisquare_test_homogeneity");
  Function octave_chebyshevpoly= Function("chebyshevpoly");
  Function octave_chebyshevT= Function("chebyshevT");
  Function octave_chebyshevU= Function("chebyshevU");
  Function octave_chol= Function("chol");
  Function octave_chol2inv= Function("chol2inv");
  Function octave_choldelete= Function("choldelete");
  Function octave_cholinsert= Function("cholinsert");
  Function octave_colinv= Function("colinv");
  Function octave_cholshift= Function("cholshift");
  Function octave_cholupdate= Function("cholupdate");
  Function octave_chop= Function("chop");
  Function octave_circshift= Function("circshift");
  Function octave_cla= Function("cla");
  Function octave_clabel= Function("clabel");
  Function octave_clc= Function("clc");
  Function octave_clf= Function("clf");
  Function octave_clock= Function("clock");
  Function octave_cloglog= Function("cloglog");
  Function octave_cmpermute= Function("cmpermute");
  Function octave_cmunique= Function("cmunique");
  Function octave_coeffs= Function("coeffs");
  Function octave_colamd= Function("colamd");
  Function octave_colloc= Function("colloc");
  Function octave_colon= Function("colon");
  Function octave_colorbar= Function("colorbar");
  Function octave_colorcube= Function("colorcube");
  Function octave_colormap= Function("colormap");
  Function octave_colperm= Function("colperm");
  Function octave_columns= Function("columns");
  Function octave_comet= Function("comet");
  Function octave_compan= Function("compan");
  Function octave_compass= Function("compass");
  Function octave_complex= Function("complex");
  Function octave_computer= Function("computer");
  Function octave_cond= Function("cond");
  Function octave_condeig= Function("condeig");
  Function octave_condest= Function("condest");
  Function octave_conj= Function("conj");
  Function octave_contour= Function("contour");
  Function octave_contour3= Function("contour3");
  Function octave_contourc= Function("contourc");
  Function octave_contourf= Function("contourf");
  Function octave_contrast= Function("contrast");
  Function octave_conv= Function("conv");
  Function octave_conv2= Function("conv2");
  Function octave_convhull= Function("convhull");
  Function octave_convhulln= Function("convhulln");
  Function octave_cool= Function("cool");
  Function octave_copper= Function("copper");
  Function octave_copyfile= Function("copyfile");
  Function octave_copyobj= Function("copyobj");
  Function octave_cor_test= Function("cor_test");
  Function octave_cos= Function("cos");
  Function octave_cosd= Function("cosd");
  Function octave_cosh= Function("cosh");
  Function octave_coshint= Function("coshint");
  Function octave_cosint= Function("cosint");
  Function octave_cot= Function("cot");
  Function octave_cotd= Function("cotd");
  Function octave_coth= Function("coth");
  Function octave_cov= Function("cov");
  Function octave_cplxpair= Function("cplxpair");
  Function octave_cputime= Function("cputime");
  Function octave_cross= Function("cross");
  Function octave_csc= Function("csc");
  Function octave_cscd= Function("cscd");
  Function octave_csch= Function("csch");
  Function octave_cstrcat= Function("cstrcat");
  Function octave_cstrcmp= Function("cstrcmp");
  Function octave_csvread= Function("csvread");
  Function octave_csvwrite= Function("csvwrite");
  Function octave_csymamd= Function("csymamd");
  Function octave_ctime= Function("ctime");
  Function octave_ctranspose= Function("ctranspose");
  Function octave_cubehelix= Function("cubehelix");
  Function octave_cummax= Function("cummax");
  Function octave_cummin= Function("cummin");
  Function octave_cumprod= Function("cumprod");
  Function octave_cumsum= Function("cumsum");
  Function octave_cumtrapz= Function("cumtrapz");
  Function octave_cylinder= Function("cylinder");
  Function octave_daspect= Function("daspect");
  Function octave_daspk= Function("daspk");
  Function octave_dasrt_options= Function("dasrt_options");
  Function octave_dassl= Function("dassl");
  Function octave_dassl_options= Function("dassl_options");
  Function octave_date= Function("date");
  Function octave_datenum= Function("datenum");
  Function octave_datestr= Function("datestr");
  Function octave_datetick= Function("datetick");
  Function octave_dawson= Function("dawson");
  Function octave_dbclear= Function("dbclear");
  Function octave_dbcont= Function("dbcont");
  Function octave_dbdown= Function("dbdown");
  Function octave_dblist= Function("dblist");
  Function octave_dblquad= Function("dblquad");
  Function octave_dbquit= Function("dbquit");
  Function octave_dbstack= Function("dbstack");
  Function octave_dbstatus= Function("dbstatus");
  Function octave_dbstep= Function("dbstep");
  Function octave_dbstop= Function("dbstop");
  Function octave_dbtype= Function("dbtype");
  Function octave_dbup= Function("dbup");
  Function octave_dbwhere= Function("dbwhere");
  Function octave_deal= Function("deal");
  Function octave_deblank= Function("deblank");
  Function octave_dec2base= Function("dec2base");
  Function octave_dec2hex= Function("dec2hex");
  Function octave_deconv= Function("deconv");
  Function octave_deg2rad= Function("deg2rad");
  Function octave_del2= Function("del2");
  Function octave_delaunay= Function("delaunay");
  Function octave_delaunayn= Function("delaunayn");
  Function octave_det= Function("det");
  Function octave_detrend= Function("detrend");
  Function octave_diag= Function("diag");
  Function octave_diff= Function("diff");
  Function octave_diffpara= Function("diffpara");
  Function octave_diffuse= Function("diffuse");
  Function octave_digits= Function("digits");
  Function octave_dilog= Function("dilog");
  Function octave_dir= Function("dir");
  Function octave_dirac= Function("dirac");
  Function octave_discrete_cdf= Function("discrete_cdf");
  Function octave_discrete_inv= Function("discrete_inv");
  Function octave_discrete_pdf= Function("discrete_pdf");
  Function octave_discrete_rnd= Function("discrete_rnd");
  Function octave_disp= Function("disp");
  Function octave_display= Function("display");
  Function octave_divergence= Function("divergence");
  Function octave_dimread= Function("dimread");
  Function octave_dimwrite= Function("dimwrite");
  Function octave_dmperm= Function("dmperm");
  Function octave_do_string_escapes= Function("do_string_escapes");
  Function octave_doc= Function("doc");
  Function octave_dot= Function("dot");
  Function octave_double= Function("double");
  Function octave_downsample= Function("downsample");
  Function octave_dsearch= Function("dsearch");
  Function octave_dsearchn= Function("dsearchn");
  Function octave_dsolve= Function("dsolve");
  Function octave_dup2= Function("dup2");
  Function octave_duplication_matrix= Function("duplication_matrix");
  Function octave_durblevinson= Function("durblevinson");
  Function octave_e= Function("e");
  Function octave_ei= Function("ei");
  Function octave_eig= Function("eig");
  Function octave_ellipke= Function("ellipke");
  Function octave_ellipsoid= Function("ellipsoid");
  Function octave_ellipticCE= Function("ellipticCE");
  Function octave_ellipticCK= Function("ellipticCK");
  Function octave_ellipticCPi= Function("ellipticCPi");
  Function octave_ellipticE= Function("ellipticE");
  Function octave_ellipticF= Function("ellipticF");
  Function octave_ellipticK= Function("ellipticK");
  Function octave_ellipticPi= Function("ellipticPi");
  Function octave_empirical_cdf= Function("empirical_cdf");
  Function octave_empirical_inv= Function("empirical_inv");
  Function octave_empirical_pdf= Function("empirical_pdf");
  Function octave_empirical_rnd= Function("empirical_rnd");
  Function octave_end= Function("end");
  Function octave_endgrent= Function("endgrent");
  Function octave_endpwent= Function("endpwent");
  Function octave_eomday= Function("eomday");
  Function octave_eps= Function("eps");
  Function octave_eq= Function("eq");
  Function octave_equationsToMatrix= Function("equationsToMatrix");
  Function octave_erf= Function("erf");
  Function octave_erfc= Function("erfc");
  Function octave_erfinv= Function("erfinv");
  Function octave_erfi= Function("erfi");
  Function octave_errno= Function("errno");
  Function octave_error= Function("error");
  Function octave_error_ids= Function("error_ids");
  Function octave_errorbar= Function("errorbar");
  Function octave_etime= Function("etime");
  Function octave_etree= Function("etree");
  Function octave_etreeplot= Function("etreeplot");
  Function octave_eulier= Function("eulier");
  Function octave_eulergamma= Function("eulergamma");
  Function octave_evalin= Function("evalin");
  Function octave_exp= Function("exp");
  Function octave_expand= Function("expand");
  Function octave_expcdf= Function("expcdf");
  Function octave_expint= Function("expint");
  Function octave_expinv= Function("expinv");
  Function octave_expm= Function("expm");
  Function octave_expm1= Function("expm1");
  Function octave_exppdf= Function("exppdf");
  Function octave_exprnd= Function("exprnd");
  Function octave_eye= Function("eye");
  Function octave_ezcontour= Function("ezcontour");
  Function octave_ezcontourf= Function("ezcontourf");
  Function octave_ezmesh= Function("ezmesh");
  Function octave_explot= Function("explot");
  Function octave_ezplot3= Function("ezplot3");
  Function octave_ezsurf= Function("ezsurf");
  Function octave_ezpolar= Function("ezpolar");
  Function octave_ezsurfc= Function("ezsurfc");
  Function octave_f_test_regression= Function("f_test_regression");
  Function octave_factor= Function("factor");
  Function octave_factorial= Function("factorial");
  Function octave_false= Function("false");
  Function octave_fcdf= Function("fcdf");
  Function octave_fclear= Function("fclear");
  Function octave_fcntl= Function("fcntl");
  Function octave_fdisp= Function("fdisp");
  Function octave_feather= Function("feather");
  Function octave_ff2n= Function("ff2n");
  Function octave_fibonacci= Function("fibonacci");
  Function octave_find= Function("find");
  Function octave_findsym= Function("findsym");
  Function octave_finiteset= Function("finiteset");
  Function octave_finv= Function("finv");
  Function octave_fix= Function("fix");
  Function octave_flintmax= Function("flintmax");
  Function octave_flip= Function("flip");
  Function octave_flipir= Function("flipir");
  Function octave_flipud= Function("flipud");
  Function octave_floor= Function("floor");
  Function octave_fminbnd= Function("fminbnd");
  Function octave_fminunc= Function("fminunc");
  Function octave_formula= Function("formula");
  Function octave_fortran= Function("fortran");
  Function octave_fourier= Function("fourier");
  Function octave_fpdf= Function("fpdf");
  Function octave_fplot= Function("fplot");
  Function octave_frac= Function("frac");
  Function octave_fractdiff= Function("fractdiff");
  Function octave_frame2im= Function("frame2im");
  Function octave_freport= Function("freport");
  Function octave_fresneic= Function("fresneic");
  Function octave_frnd= Function("frnd");
  Function octave_fskipl= Function("fskipl");
  Function octave_fsolve= Function("fsolve");
  Function octave_full= Function("full");
  Function octave_fwhm= Function("fwhm");
  Function octave_fzero= Function("fzero");
  Function octave_gallery= Function("gallery");
  Function octave_gamcdf= Function("gamcdf");
  Function octave_gaminv= Function("gaminv");
  Function octave_gamma= Function("gamma");
  Function octave_gammainc= Function("gammainc");
  Function octave_gammaln= Function("gammaln");
  Function octave_gca= Function("gca");
  Function octave_gcbf= Function("gcbf");
  Function octave_gcbo= Function("gcbo");
  Function octave_gcd= Function("gcd");
  Function octave_ge= Function("ge");
  Function octave_geocdf= Function("geocdf");
  Function octave_geoinv= Function("geoinv");
  Function octave_geopdf= Function("geopdf");
  Function octave_geornd= Function("geornd");
  Function octave_givens= Function("givens");
  Function octave_glpk= Function("glpk");
  Function octave_gmres= Function("gmres");
  Function octave_gmtime= Function("gmtime");
  Function octave_gnplot_binary= Function("gnplot_binary");
  Function octave_gplot= Function("gplot");
  Function octave_gradient= Function("gradient");
  Function octave_gray= Function("gray");
  Function octave_gray2ind= Function("gray2ind");
  Function octave_gt= Function("gt");
  Function octave_gunzip= Function("gunzip");
  Function octave_gzip= Function("gzip");
  Function octave_hadamard= Function("hadamard");
  Function octave_hankel= Function("hankel");
  Function octave_harmonic= Function("harmonic");
  Function octave_has= Function("has");
  Function octave_hash= Function("hash");
  Function octave_heaviside= Function("heaviside");
  Function octave_help= Function("help");
  Function octave_hess= Function("hess");
  Function octave_hex2dec= Function("hex2dec");
  Function octave_hex2num= Function("hex2num");
  Function octave_hilb= Function("hilb");
  Function octave_hilbert_curve= Function("hilbert_curve");
  Function octave_hist= Function("hist");
  Function octave_horner= Function("horner");
  Function octave_horzcat= Function("horzcat");
  Function octave_hot= Function("hot");
  Function octave_housh= Function("housh");
  Function octave_hsv2rgb= Function("hsv2rgb");
  Function octave_hurst= Function("hurst");
  Function octave_hygecdf= Function("hygecdf");
  Function octave_hygeinv= Function("hygeinv");
  Function octave_hygepdf= Function("hygepdf");
  Function octave_hygernd= Function("hygernd");
  Function octave_hypergeom= Function("hypergeom");
  Function octave_hypot= Function("hypot");
  Function octave_I= Function("I");
  Function octave_ichol= Function("ichol");
  Function octave_idist= Function("idist");
  Function octave_idivide= Function("idivide");
  Function octave_igamma= Function("igamma");
  Function octave_ilaplace= Function("ilaplace");
  Function octave_ilu= Function("ilu");
  Function octave_im2double= Function("im2double");
  Function octave_im2frame= Function("im2frame");
  Function octave_im2int16= Function("im2int16");
  Function octave_im2single= Function("im2single");
  Function octave_im2uint16= Function("im2uint16");
  Function octave_im2uint8= Function("im2uint8");
  Function octave_imag= Function("imag");
  Function octave_image= Function("image");
  Function octave_imagesc= Function("imagesc");
  Function octave_imfinfo= Function("imfinfo");
  Function octave_imformats= Function("imformats");
  Function octave_importdata= Function("importdata");
  Function octave_imread= Function("imread");
  Function octave_imshow= Function("imshow");
  Function octave_imwrite= Function("imwrite");
  Function octave_ind2gray= Function("ind2gray");
  Function octave_ind2rgb= Function("ind2rgb");
  Function octave_int2sub= Function("int2sub");
  Function octave_index= Function("index");
  Function octave_inf= Function("Inf");
  Function octave_inpolygon= Function("inpolygon");
  Function octave_input= Function("input");
  Function octave_interp1= Function("interp1");
  Function octave_interp2= Function("interp2");
  Function octave_interp3= Function("interp3");
  Function octave_intersect= Function("intersect");
  Function octave_intmin= Function("intmin");
  Function octave_inv= Function("inv");
  Function octave_invhilb= Function("invhilb");
  Function octave_inimpinvar= Function("inimpinvar");
  Function octave_ipermute= Function("ipermute");
  Function octave_iqr= Function("iqr");
  Function octave_isa= Function("isa");
  Function octave_isequal= Function("isequal");
  Function octave_ishermitian= Function("ishermitian");
  Function octave_isprime= Function("isprime");
  Function octave_jit_enable= Function("jit_enable");
  Function octave_kbhit= Function("kbhit");
  Function octave_kendall= Function("kendall");
  Function octave_kron= Function("kron");
  Function octave_kurtosis= Function("kurtosis");
  Function octave_laplace= Function("laplace");
  Function octave_laplace_cdf= Function("laplace_cdf");
  Function octave_laplace_inv= Function("laplace_inv");
  Function octave_laplace_pdf= Function("laplace_pdf");
  Function octave_laplace_rnd= Function("laplace_rnd");
  Function octave_laplacian= Function("laplacian");
  Function octave_lcm= Function("lcm");
  Function octave_ldivide= Function("ldivide");
  Function octave_le= Function("le");
  Function octave_legendre= Function("legendre");
  Function octave_length= Function("length");
  Function octave_lgamma= Function("lgamma");
  Function octave_limit= Function("limit");
  Function octave_line= Function("line");
  Function octave_linprog= Function("linprog");
  Function octave_linsolve= Function("linsolve");
  Function octave_linspace= Function("linspace");
  Function octave_load= Function("load");
  Function octave_log= Function("log");
  Function octave_log10= Function("log10");
  Function octave_log1p= Function("log1p");
  Function octave_log2= Function("log2");
  Function octave_logical= Function("logical");
  Function octave_logistic_cdf= Function("logistic_cdf");
  Function octave_logistic_inv= Function("logistic_inv");
  Function octave_logistic_pdf= Function("logistic_pdf");
  Function octave_logistic_regression= Function("logistic_regression");
  Function octave_logit= Function("logit");
  Function octave_loglog= Function("loglog");
  Function octave_loglogerr= Function("loglogerr");
  Function octave_logm= Function("logm");
  Function octave_logncdf= Function("logncdf");
  Function octave_logninv= Function("logninv");
  Function octave_lognpdf= Function("lognpdf");
  Function octave_lognrnd= Function("lognrnd");
  Function octave_lognspace= Function("lognspace");
  Function octave_lookup= Function("lookup");
  Function octave_lscov= Function("lscov");
  Function octave_lsode= Function("lsode");
  Function octave_lsqnonneg= Function("lsqnonneg");
  Function octave_lt= Function("lt");
  Function octave_magic= Function("magic");
  Function octave_manova= Function("manova");
  Function octave_minus= Function("minus");
  Function octave_mkpp= Function("mkpp");
  Function octave_mldivide= Function("mldivide");
  Function octave_mod= Function("mod");
  Function octave_moment= Function("moment");
  Function octave_mpoles= Function("mpoles");
  Function octave_mpower= Function("mpower");
  Function octave_mrdivide= Function("mrdivide");
  Function octave_mu2lin= Function("mu2lin");
  Function octave_na= Function("NA");
  Function octave_nan= Function("NaN");
  Function octave_nextpow2= Function("nextpow2");
  Function octave_nnz= Function("nnz");
  Function octave_nonzeros= Function("nonzeros");
  Function octave_norm= Function("norm");
  Function octave_normcdf= Function("normcdf");
  Function octave_normest= Function("normest");
  Function octave_normest1= Function("normest1");
  Function octave_norminv= Function("norminv");
  Function octave_normpdf= Function("normpdf");
  Function octave_normrnd= Function("normrnd");
  Function octave_nth_element= Function("nth_element");
  Function octave_nth_root= Function("nth_root");
  Function octave_null= Function("null");
  Function octave_numel= Function("numel");
  Function octave_ode23= Function("ode23");
  Function octave_ode45= Function("ode45");
  Function octave_ols= Function("ols");
  Function octave_ones= Function("ones");
  Function octave_prod= Function("prod");
  Function octave_power= Function("power");
  Function octave_sin= Function("sin");
  Function octave_sqrt= Function("sqrt");
  Function octave_sum= Function("sum");
  Function octave_sumsq= Function("sumsq");
  Function octave_tan= Function("tan");
  Function octave_tanh= Function("tanh");
  Function octave_sinh= Function("sinh");
  Function octave_bin_values= Function("bin_values");
  Function octave_catmullrom= Function("catmullrom");
  Function octave_csape= Function("csape");
  Function octave_csapi= Function("csapi");
  Function octave_csaps= Function("csaps");
  Function octave_csaps_sel= Function("csaps_sel");
  Function octave_dedup= Function("dedup");
  Function octave_fnder= Function("fnder");
  Function octave_fnplt= Function("fnplt");
  Function octave_fnval= Function("fnval");
  Function octave_regularization= Function("regularization");
  Function octave_regularization2D= Function("regularization2D");
  Function octave_tpaps= Function("tpaps");
  Function octave_tps_val= Function("tps_val");
  Function octave_tps_val_der= Function("tps_val_der");
  //Function octave_rms= Function("rms");
  Function octave_normalize= Function("normalize");
  Function octave_gaindb= Function("gaindb");
  Function octave_crestfactor= Function("crestfactor");
  Function octave_uquant= Function("uquant");
  Function octave_firwin= Function("firwin");
  Function octave_firkaiser= Function("firkaiser");
  Function octave_fir2long= Function("fir2long");
  Function octave_long2fir= Function("long2fir");
  Function octave_freqwin= Function("freqwin");
  Function octave_firfilter= Function("firfilter");
  Function octave_blfilter= Function("blfilter");
  Function octave_warpedblfilter= Function("warpedblfilter");
  Function octave_freqfilter= Function("freqfilter");
  Function octave_pfilt= Function("pfilt");
  Function octave_magresp= Function("magresp");
  Function octave_transferfunction= Function("transferfunction");
  Function octave_pgrdelay= Function("pgrdelay");
  Function octave_rampup= Function("rampup");
  Function octave_rampdown= Function("rampdown");
  Function octave_thresh= Function("thresh");
  Function octave_largestr= Function("largestr");
  Function octave_largestn= Function("largestn");
  Function octave_dynlimit= Function("dynlimit");
  Function octave_groupthresh= Function("groupthresh");
  Function octave_rgb2jpeg= Function("rgb2jpeg");
  Function octave_jpeg2rgb= Function("jpeg2rgb");
  Function octave_qam4= Function("qam4");
  Function octave_iqam4= Function("iqam4");
  Function octave_semiaudplot= Function("semiaudplot");
  Function octave_audtofreq= Function("audtofreq");
  Function octave_freqtoaud= Function("freqtoaud");
  Function octave_audspace= Function("audspace");
  Function octave_audspacebw= Function("audspacebw");
  Function octave_erbtofreq= Function("erbtofreq");
  Function octave_freqtoerb= Function("freqtoerb");
  Function octave_erbspace= Function("erbspace");
}

Octopus::OctopusValueList* convert_octave_value_list(lua_State * L)
{
  Octopus::OctopusValueList* ptr;
  SWIG_ConvertPtr(L,-1,(void**)&ptr,SWIGTYPE_p_Octopus__OctopusValueList,0);
  return ptr;
}

%}