
#include <iostream>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
#include <Eigen/Core>

using ArrayXf = Array<float>;
using ArrayXd = Array<double>;
using ArrayXcf = Array<std::complex<float>>;
using ArrayXcd = Array<std::complex<double>>;
using VectorXf = FloatRowVector;
using VectorXd = RowVector;
using VectorXcf= FloatComplexRowVector;
using VectorXcd= ComplexRowVector;
//using ColVectorXf = FloatColVector;
//using ColVectorXd = ColVector;
//using ColVectorXcf= FloatComplexColVector;
//using ColVectorXcd= ComplexColVector;
using MatrixXf = FloatMatrix;
using MatrixXd = Matrix;
using MatrixXcf= FloatComplexMatrix;
using MatrixXcd= ComplexMatrix;
using Value=octave_value;
using ValueList=octave_value_list;

namespace Octave
{
  struct Application : public octave::application
  {
    Application() {
      forced_interactive(true);
    }
    int execute() {
      return 0;
    }
  };

  struct Octopus
  {   
      octave::interpreter *interpreter;
      Application pita;    

      Octopus() {                  
        interpreter = new octave::interpreter();
        interpreter->interactive(false);
        interpreter->initialize_history(false);       
        interpreter->initialize();            
        interpreter->execute();
        std::string path = ".";
        octave_value_list p;
        p(0) = path;
        octave_value_list o1 = interpreter->feval("addpath", p, 1);            
        run_script("startup.m");
      }
      ~Octopus()
      {
          if(interpreter) delete interpreter;
      }
      
      void run_script(const std::string& s) {
        octave::source_file(s);
      }
      
      ValueList eval_string(std::string func, bool silent=false, int noutputs=1)
      {          
        octave_value_list out =interpreter->eval_string(func.c_str(), silent, noutputs);
        return out;
      }
      ValueList eval(std::string func, ValueList inputs, int noutputs=1)
      {          
        octave_value_list out =interpreter->feval(func.c_str(), inputs, noutputs);
        return out;
      }

      ValueList operator()(std::string func, ValueList inputs, int noutputs=1)
      {
        return eval(func,inputs,noutputs);
      }
      
  };

  Eigen::VectorXf Eigenize(const VectorXf & m)
  {
      Eigen::VectorXf r(m.cols());
      for(size_t i = 0; i < m.cols(); i++)        
              r(i) = m(i);
      return r;
  }
  Eigen::VectorXd Eigenize(const VectorXd & m)
  {
      Eigen::VectorXd r(m.cols());
      for(size_t i = 0; i < m.cols(); i++)        
              r(i) = m(i);
      return r;
  }
  Eigen::VectorXcf Eigenize(const VectorXcf & m)
  {
      Eigen::VectorXcf r(m.cols());
      for(size_t i = 0; i < m.cols(); i++)        
              r(i) = m(i);
      return r;
  }
  Eigen::VectorXcd Eigenize(const VectorXcd & m)
  {
      Eigen::VectorXcd r(m.cols());
      for(size_t i = 0; i < m.rows(); i++)        
              r(i) = m(i);
      return r;
  }

  Eigen::MatrixXf Eigenize(const MatrixXf & m)
  {
      Eigen::MatrixXf r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }
  Eigen::MatrixXd Eigenize(const MatrixXd & m)
  {
      Eigen::MatrixXd r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }
  Eigen::MatrixXcf Eigenize(const MatrixXcf & m)
  {
      Eigen::MatrixXcf r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }
  Eigen::MatrixXcd Eigenize(const MatrixXcd & m)
  {
      Eigen::MatrixXcd r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }

  VectorXf Octavate(const Eigen::VectorXf & m)
  {
      VectorXf r(m.rows());
      for(size_t i = 0; i < m.rows(); i++)        
              r(i) = m(i);
      return r;
  }
  VectorXd Octavate(const Eigen::VectorXd & m)
  {
      VectorXd r(m.rows());
      for(size_t i = 0; i < m.rows(); i++)        
              r(i) = m(i);
      return r;
  }
  VectorXcf Octavate(const Eigen::VectorXcf & m)
  {
      VectorXcf r(m.rows());
      for(size_t i = 0; i < m.rows(); i++)        
              r(i) = m(i);
      return r;
  }
  VectorXcd Octavate(const Eigen::VectorXcd & m)
  {
      VectorXcd r(m.rows());
      for(size_t i = 0; i < m.rows(); i++)        
              r(i) = m(i);
      return r;
  }

  MatrixXf Octavate(const Eigen::MatrixXf & m)
  {
      MatrixXf r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }
  MatrixXd Octavate(const Eigen::MatrixXd & m)
  {
      MatrixXd r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }
  MatrixXcf Octavate(const Eigen::MatrixXcf & m)
  {
      MatrixXcf r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }
  MatrixXcd Octavate(const Eigen::MatrixXcd & m)
  {
      MatrixXcd r(m.rows(),m.cols());
      for(size_t i = 0; i < m.rows(); i++)
          for(size_t j = 0; j < m.cols(); j++)
              r(i,j) = m(i,j);
      return r;
  }


  struct Function
  {
    std::string name;
    

    Function(const std::string& f) : name(f) {}

    ValueList operator()()
    {
        ValueList input;
        int num_outputs=0;
        return octave::feval(name.c_str(),input,num_outputs);
    }

    ValueList operator()(ValueList & input, int num_outputs=1)
    {
        return octave::feval(name.c_str(),input,num_outputs);
    }
  };

    // why no?
    //#define def(x) Function octave_##x(#x)  
    #define def(x) Function x(#x)  
  
    def(fft);
    def(ifft);
    def(fft2);
    def(ifft2);
    def(fftconv);
    def(fftfilt);
    def(fftn);
    def(fftshift);
    def(fftw);
    def(ifftn);
    def(ifftshift);
    def(ifht);
    def(ifourier);
    def(ifwht);
    def(ifwt);
    def(ifwt2);
    def(buffer);
    def(chirp);
    def(cmorwavf);  
    def(gauspuls);
    def(gmonopuls);
    def(mexihat);
    def(meyeraux);  
    def(morlet);
    def(pulstran);
    def(rectpuls);
    def(sawtooth);
    def(shanwavf);
    def(shiftdata);
    def(sigmoid_train);
    def(specgram);
    def(square);
    def(tripuls);
    def(udecode);
    def(uencoder);
    def(unshiftdata);
    def(findpeaks);
    def(peak2peak);
    def(peak2rms);
    def(rms);
    def(rssq);
    def(cconv);
    def(convmtx);  
    def(wconv);
    def(xcorr);
    def(xcorr2);
    def(xcov);
    def(filtfilt);
    def(fltic);
    def(medfilt1);
    def(movingrms);
    def(sgolayfilt);
    def(sosfilt);
    def(freqs);
    def(freqs_plot);
    def(freqz);
    def(freqz_plot);
    def(impz);
    def(zplane);
    def(filter);
    def(filter2);  
    def(fir1);
    def(fir2);
    def(firls);
    def(sinc);
    def(unwrap);

    def(bartlett);
    def(blackman);  
    def(blackmanharris);
    def(blackmannuttal);
    def(dftmtx);
    def(hamming);
    def(hann);
    def(hanning);
    def(pchip);
    def(periodogram);
    def(sinetone);
    def(sinewave);
    def(spectral_adf);
    def(spectral_xdf);
    def(spencer);
    def(stft);
    def(synthesis);
    def(yulewalker);
    def(polystab);
    def(residued);
    def(residuez);
    def(sos2ss);
    def(sos2tf);
    def(sos2zp);
    def(ss2tf);
    def(ss2zp);
    def(tf2sos);
    def(tf2ss);
    def(tf2zp);
    def(zp2sos);
    def(zp2ss);
    def(zp2tf);
    def(besselap);
    def(besself);
    def(bilinear);
    def(buttap);
    def(butter);
    def(buttord);
    def(cheb);
    def(cheb1ap);
    def(cheb1ord);
    def(cheb2ap);
    def(cheb2ord);
    def(chebywin);
    def(cheby1);
    def(cheby2);
    def(ellip);
    def(ellipap);  
    def(ellipord);
    def(impinvar);
    def(ncauer);
    def(pei_tseng_notch);
    def(sftrans);
    def(cl2bp);
    def(kaiserord);
    def(qp_kaiser);
    def(remez);
    def(sgplay);
    def(bitrevorder);
    def(cceps);
    def(cplxreal);
    def(czt);
    def(dct);
    def(dct2);  
    def(dctmtx);
    def(digitrevorder);
    def(dst);
    def(dwt);
    def(rceps);
    def(ar_psd);
    def(cohere);
    def(cpsd);
    def(csd);
    def(db2pow);
    def(mscohere);
    def(pburg);
    def(pow2db);
    def(pwelch);
    def(pyulear);
    def(tfe);
    def(tfestimate);
    def(__power);
    def(barthannwin);
    def(bohmanwin);
    def(boxcar);
    def(flattopwin);
    def(chebwin);
    def(gaussian);
    def(gausswin);
    def(kaiser);  
    def(nuttalwin);
    def(parzenwin);
    def(rectwin);
    def(tukeywin);
    def(ultrwin);
    def(welchwin);
    def(window);
    def(arburg);
    def(aryule);
    def(invfreq);
    def(invfreqz);
    def(invfreqs);
    def(levinson);
    def(data2fun);
    def(decimate);
    //def(interp);
    def(resample);
    def(upfirdn);
    def(upsample);
    def(clustersegment);
    def(fracshift);
    def(marcumq);
    def(primitive);
    def(sampled2continuous);
    def(schtrig);
    def(upsamplefill);
    def(wkeep);
    def(wrev);
    def(zerocrossing);


    def(fht);
    def(fwht);  
    def(hilbert);
    def(idct);
    def(idct2);

    def(max);
    def(mean);
    def(meansq);
    def(median);
    def(min);

    def(plot);
    def(pause);

    def(abs);
    def(accumarray);
    def(accumdim);
    def(acos);
    def(acosd);
    def(acosh);
    def(acot);
    def(acotd);
    def(acoth);
    def(acsc);
    def(acsch);
    def(acscd);
    def(airy);
    def(adjoint);
    def(all);
    def(allow_non_integer_range_as_index);
    def(amd);
    def(ancestor);
    //def(and);
    def(angle);
    def(annotation);
    def(anova);
    def(ans);
    def(any);    
    def(arch_fit);
    def(arch_rnd);
    def(arch_test);
    def(area);
    def(arg);
    def(arrayfun);  
    def(asec);
    def(asecd);
    def(asech);
    def(asin);
    def(asind);
    def(asinh);
    def(assume);
    def(assumptions);
    def(atan);
    def(atand);
    def(atanh);
    def(atan2);
    def(audiodevinfo);
    def(audioformats);
    def(audioinfo);
    def(audioread);
    def(audiowrite);
    def(autoreg_matrix);
    def(autumn);
    def(axes);
    def(axis);
    def(balance);
    def(bandwidth);


    def(bar);
    def(barh);
    def(bathannwin);
    def(bartlett_test);
    def(base2dec);
    def(base64_decode);
    def(base64_encode);
    def(beep);
    def(beep_on_error);
    def(bernoulli);  
    def(besseli);
    def(besseljn);
    def(besselk);
    def(bessely);
    def(beta);
    def(betacdf);
    def(betainc);
    def(betaincinv);
    def(betainv);
    def(betain);
    def(betapdf);
    def(betarnd);
    def(bicg);
    def(bicgstab);  
    def(bin2dec);
    def(bincoeff);
    def(binocdf);
    def(binoinv);
    def(binopdf);
    def(binornd);
    //def(bitand);
    def(bitcmp);
    def(bitget);
    //def(bitor);
    def(bitpack);  
    def(bitset);
    def(bitshift);
    def(bitunpack);
    def(bitxor);
    def(blanks);
    def(blkdiag);
    def(blkmm);
    def(bone);
    def(box);  
    def(brighten);
    def(bsxfun);
    def(builtin);
    def(bzip2);

    def(calendar);
    def(camlight);
    def(cart2pol);
    def(cart2sph);
    def(cast);
    def(cat);
    def(catalan);
    def(cauchy);
    def(cauchy_cdf);
    def(cauchy_inv);
    def(cauchy_pdf);
    def(cauchy_rnd);
    def(caxis);
    def(cbrt);  
    def(ccode);
    def(ccolamd);  
    def(ceil);
    def(center);
    def(centroid);
    def(cgs);  
    def(chi2cdf);
    def(chi2inv);
    def(chi2pdf);
    def(chi2rnd);
    def(children);  
    def(chisquare_test_homogeneity);  
    def(chebyshevpoly);
    def(chebyshevT);
    def(chebyshevU);
    def(chol);
    def(chol2inv);
    def(choldelete);
    def(cholinsert);
    def(colinv);
    def(cholshift);
    def(cholupdate);
    def(chop);
    def(circshift);  
    def(cla);
    def(clabel);
    def(clc);
    def(clf);
    def(clock);
    def(cloglog);  
    def(cmpermute);
    def(cmunique);
    def(coeffs);  
    def(colamd);
    def(colloc);
    def(colon);
    def(colorbar);
    def(colorcube);
    def(colormap);
    def(colperm);
    def(columns);
    def(comet);
    def(compan);
    def(compass);
    def(complex);
    def(computer);
    def(cond);
    def(condeig);
    def(condest);
    def(conj);
    def(contour);
    def(contour3);
    def(contourc);
    def(contourf);
    def(contrast);
    def(conv);
    def(conv2);
    def(convhull);
    def(convhulln);  
    def(cool);
    def(copper);
    def(copyfile);
    def(copyobj);
    def(cor_test);
    def(cos);
    def(cosd);
    def(cosh);
    def(coshint);
    def(cosint);
    def(cot);
    def(cotd);
    def(coth);
    def(cov);
    def(cplxpair);    
    def(cputime);
    def(cross);
    def(csc);
    def(cscd);
    def(csch);  
    def(cstrcat);
    def(cstrcmp);
    def(csvread);
    def(csvwrite);
    def(csymamd);
    def(ctime);
    def(ctranspose);
    def(cubehelix);
    def(cummax);
    def(cummin);
    def(cumprod);
    def(cumsum);
    def(cumtrapz);
    def(cylinder);

    def(daspect);
    def(daspk);
    def(dasrt_options);
    def(dassl);
    def(dassl_options);  
    def(date);
    def(datenum);
    def(datestr);
    def(datetick);
    def(dawson);  
    def(dbclear);
    def(dbcont);
    def(dbdown);
    def(dblist);
    def(dblquad);
    def(dbquit);
    def(dbstack);
    def(dbstatus);
    def(dbstep);
    def(dbstop);
    def(dbtype);
    def(dbup);
    def(dbwhere);  
    def(deal);
    def(deblank);
    def(dec2base);
    def(dec2hex);  
    def(deconv);
    def(deg2rad);
    def(del2);
    def(delaunay);
    def(delaunayn);  
    def(det);
    def(detrend);  
    def(diag);
    def(diff);
    def(diffpara);
    def(diffuse);  
    def(digits);
    def(dilog);
    def(dir);
    def(dirac);  
    def(discrete_cdf);
    def(discrete_inv);
    def(discrete_pdf);
    def(discrete_rnd);
    def(disp);
    def(display);
    def(divergence);
    def(dimread);
    def(dimwrite);
    def(dmperm);
    def(do_string_escapes);
    def(doc);
    def(dot);
    //def(double);
    def(downsample);
    def(dsearch);
    def(dsearchn);
    def(dsolve);  
    def(dup2);
    def(duplication_matrix);
    def(durblevinson);


    def(e);
    def(ei);
    def(eig);
    def(ellipke);  
    def(ellipsoid);
    def(ellipticCE);
    def(ellipticCK);
    def(ellipticCPi);
    def(ellipticE);
    def(ellipticF);
    def(ellipticK);
    def(ellipticPi);
    def(empirical_cdf);
    def(empirical_inv);
    def(empirical_pdf);
    def(empirical_rnd);
    def(end);
    def(endgrent);
    def(endpwent);
    def(eomday);
    def(eps);
    def(eq);
    def(equationsToMatrix);
    def(erf);
    def(erfc);
    def(erfinv);
    def(erfi);
    //def(errno);
    def(error);
    def(error_ids);
    def(errorbar);
    def(etime);
    def(etree);
    def(etreeplot);
    def(eulier);
    def(eulergamma);
    def(evalin);
    def(exp);
    def(expand);
    def(expcdf);
    def(expint);
    def(expinv);
    def(expm);
    def(expm1);
    def(exppdf);
    def(exprnd);
    def(eye);
    def(ezcontour);
    def(ezcontourf);
    def(ezmesh);
    def(explot);
    def(ezplot3);
    def(ezsurf);
    def(ezpolar);
    def(ezsurfc);

    def(f_test_regression);
    def(factor);
    def(factorial);
    //def(false);
    def(fcdf);
    def(fclear);
    def(fcntl);
    def(fdisp);
    def(feather);
    def(ff2n);  
    def(fibonacci);  
    def(find);  
    def(findsym);
    def(finiteset);
    def(finv);
    def(fix);  
    def(flintmax);
    def(flip);
    def(flipir);
    def(flipud);
    def(floor);
    def(fminbnd);
    def(fminunc);
    def(formula);
    def(fortran);  
    def(fourier);
    def(fpdf);
    def(fplot);
    def(frac);
    def(fractdiff);
    def(frame2im);
    def(freport);  
    def(fresneic);
    def(frnd);
    def(fskipl);
    def(fsolve);
    def(full);
    def(fwhm);  
    def(fzero);

    def(gallery);
    def(gamcdf);
    def(gaminv);
    def(gamma);
    def(gammainc);
    def(gammaln);    
    def(gca);
    def(gcbf);
    def(gcbo);
    def(gcd);
    def(ge);
    def(geocdf);
    def(geoinv);
    def(geopdf);
    def(geornd);
    def(givens);
    def(glpk);  
    def(gmres);
    def(gmtime);
    def(gnplot_binary);
    def(gplot);
    def(gradient);
    def(gray);
    def(gray2ind);
    def(gt);
    def(gunzip);
    def(gzip);

    def(hadamard);  
    def(hankel);  
    def(harmonic);
    def(has);
    def(hash);
    def(heaviside);
    def(help);
    def(hess);
    def(hex2dec);
    def(hex2num);
    def(hilb);  
    def(hilbert_curve);
    def(hist);
    def(horner);
    def(horzcat);
    def(hot);
    def(housh);
    def(hsv2rgb);
    def(hurst);
    def(hygecdf);
    def(hygeinv);
    def(hygepdf);
    def(hygernd);
    def(hypergeom);
    def(hypot);

    def(I);
    def(ichol);  
    def(idist);
    def(idivide);  
    def(igamma);  
    def(ilaplace);
    def(ilu);
    def(im2double);
    def(im2frame);
    def(im2int16);
    def(im2single);
    def(im2uint16);
    def(im2uint8);
    def(imag);
    def(image);
    def(imagesc);
    def(imfinfo);
    def(imformats);
    def(importdata);  
    def(imread);
    def(imshow);
    def(imwrite);
    def(ind2gray);
    def(ind2rgb);
    def(int2sub);
    def(index);
    def(Inf);
    def(inpolygon);
    def(input);  
    def(interp1);
    def(interp2);
    def(interp3);
    def(intersect);
    def(intmin);
    def(inv);
    def(invhilb);
    def(inimpinvar);
    def(ipermute);
    def(iqr);
    def(isa);
    def(isequal);
    def(ishermitian);
    def(isprime);

    def(jit_enable);

    def(kbhit);
    def(kendall);
    def(kron);
    def(kurtosis);

    def(laplace);
    def(laplace_cdf);
    def(laplace_inv);
    def(laplace_pdf);
    def(laplace_rnd);
    def(laplacian);
    def(lcm);
    def(ldivide);
    def(le);
    def(legendre);
    def(length);
    def(lgamma);
    def(limit);
    def(line);
    def(linprog);
    def(linsolve);
    def(linspace);
    def(load);
    def(log);
    def(log10);
    def(log1p);
    def(log2);
    def(logical);
    def(logistic_cdf);
    def(logistic_inv);
    def(logistic_pdf);
    def(logistic_regression);
    def(logit);
    def(loglog);
    def(loglogerr);
    def(logm);
    def(logncdf);
    def(logninv);
    def(lognpdf);
    def(lognrnd);
    def(lognspace);
    def(lookup);
    def(lscov);
    def(lsode);
    def(lsqnonneg);
    def(lt);

    def(magic);
    def(manova);
    def(minus);
    def(mkpp);
    def(mldivide);
    def(mod);
    def(moment);    
    def(mpoles);
    def(mpower);
    def(mrdivide);
    def(mu2lin);

    def(NA);
    def(NaN);
    def(nextpow2);
    def(nnz);
    def(nonzeros);
    def(norm);
    def(normcdf);
    def(normest);
    def(normest1);
    def(norminv);
    def(normpdf);
    def(normrnd);
    def(nth_element);
    def(nth_root);
    def(null);
    def(numel);

    def(ode23);
    def(ode45);
    def(ols);
    def(ones);

    def(prod);
    def(sin);
    def(sqrt);
    def(sum);
    def(sumsq);

    def(tan);
    def(tanh);
    def(sinh);




// image
// fuzzy-logic-toolkit

/*  
// splines  
namespace splines
{  
  def(bin_values);
  def(catmullrom);  
  def(csape);
  def(csapi);
  def(csaps);
  def(csaps_sel);
  def(dedup);
  def(fnder);
  def(fnplt);
  def(fnval);
  def(regularization);
  def(regularization2D);
  def(tpaps);
  def(tps_val);
  def(tps_val_der);
}
*/
  /* ltfat = not installed
namspace ltfat {  
  def(rms);
  def(normalize);
  def(gaindb);
  def(crestfactor);
  def(uquant);
  def(firwin);
  def(firkaiser);
  def(fir2long);
  def(long2fir);
  def(freqwin);
  def(firfilter);
  def(blfilter);
  def(warpedblfilter);
  def(freqfilter);
  def(pfilt);
  def(magresp);
  def(transferfunction);
  def(pgrdelay);
  def(rampup);
  def(rampdown);
  def(thresh);
  def(largestr);
  def(largestn);
  def(dynlimit);
  def(groupthresh);
  def(rgb2jpeg);
  def(jpeg2rgb);
  def(qam4);
  def(iqam4);
  def(semiaudplot);
  def(audtofreq);
  def(freqtoaud);
  def(audspace);
  def(audspacebw);
  def(erbtofreq);
  def(freqtoerb);
  def(erbspace);
  */
};