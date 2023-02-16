%module arma
%{
#include "SimpleArmadillo.h"
using namespace Matlab;
%}
%include "std_complex.i"

%include "SimpleArmadillo.h"

%template(float_complex) std::complex<float>;
%template(double_domplex) std::complex<double>;

%template(flt_mat) Matlab::Matrix<float>;
%template(dbl_mat) Matlab::Matrix<double>;
%template(fltcomplex_mat) Matlab::Matrix<std::complex<float>>;
%template(dblcomplex_mat) Matlab::Matrix<std::complex<double>>;
%template(int_mat) Matlab::Matrix<int>;
%template(uint_mat) Matlab::Matrix<unsigned>;
%template(long_mat) Matlab::Matrix<long>;
%template(ulong_mat) Matlab::Matrix<unsigned long>;
%template(llong_mat) Matlab::Matrix<long long>;
%template(ullong_mat) Matlab::Matrix<unsigned long long>;

%template(flt_vec) Matlab::Vector<float>;
%template(dbl_vec) Matlab::Vector<double>;
%template(fltcomplex_vec) Matlab::Vector<std::complex<float>>;
%template(dblcomplex_vec) Matlab::Vector<std::complex<double>>;


%typemap(in) arma::vec {
    $1 = Matlab::Vector<double>;
}
%typemap(out) arma::vec {
    $result = Matlab::Vector<double>($1);
}
%typemap(in) arma::mat {
    $1 = Matlab::Matrix<double>;
}
%typemap(out) arma::mat {
    $result = Matlab::Matrix<double>($1);
}
%typemap(in) arma::dmat {
    $1 = Matlab::Matrix<double>;
}
%typemap(out) arma::dmat {
    $result = Matlab::Matrix<double>($1);
}
%typemap(in) arma::fmat {
    $1 = Matlab::Matrix<float>;
}
%typemap(out) arma::fmat {
    $result = Matlab::Matrix<float>($1);
}
%typemap(in) arma::cx_mat {
    $1 = Matlab::Matrix<std::complex<double>>;
}
%typemap(out) arma::cx_mat {
    $result = Matlab::Matrix<std::complex<double>>($1);
}
%typemap(in) arma::cx_dmat {
    $1 = Matlab::Matrix<std::complex<double>>;
}
%typemap(out) arma::cx_dmat {
    $result = Matlab::Matrix<std::complex<double>>($1);
}
%typemap(in) arma::cx_fmat {
    $1 = Matlab::Matrix<std::complex<float>>;
}
%typemap(out) arma::cx_fmat {
    $result = Matlab::Matrix<std::complex<float>>($1);
}
%typemap(in) arma::umat {
    $1 = Matlab::Matrix<unsigned>;
}
%typemap(out) arma::umat {
    $result = Matlab::Matrix<unsigned>($1);
}
%typemap(in) arma::imat {
    $1 = Matlab::Matrix<int>;
}
%typemap(out) arma::imat {
    $result = Matlab::Matrix<int>($1);
}