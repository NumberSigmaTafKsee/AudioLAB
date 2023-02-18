
%{
#include "Eigen.h"
using namespace SimpleEigen;
%}

%include "eigen.i"

%template(ArrayXf) SimpleEigen::Array<float>;
%template(ArrayXd) SimpleEigen::Array<double>;
%template(ArrayXi) SimpleEigen::Array<int>;

%template(ArrayXXf) SimpleEigen::Array2D<float>;
%template(ArrayXXd) SimpleEigen::Array2D<double>;
%template(ArrayXXi) SimpleEigen::Array2D<int>;

%template(VectorXf) SimpleEigen::ColVector<float>;
%template(VectorXd) SimpleEigen::ColVector<double>;
%template(VectorXi) SimpleEigen::ColVector<int>;

%template(RowVectorXf) SimpleEigen::RowVector<float>;
%template(RowVectorXd) SimpleEigen::RowVector<double>;
%template(RowVectorXi) SimpleEigen::RowVector<int>;

%template(MatrixXf) SimpleEigen::Matrix<float>;
%template(MatrixXd) SimpleEigen::Matrix<double>;
%template(MatrixXi) SimpleEigen::Matrix<int>;
