%module eigen
%{
#define EIGEN_USE_BLAS 1
#define EIGEN_USE_LAPACKE 1
#include <Eigen/Core>
#include <complex>
#include "eigen-array.hpp"
#include "eigen-vector.hpp"
//#include "eigen-vectorwiseop.hpp"
#include <cfloat>
%}
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_string.i"


%include "eigen-array.i"
%include "eigen-matrix.i"
%include "eigen-array.hpp"
%include "eigen-vector.hpp"
//%include "eigen-vectorwiseop.hpp"

%template(Vectorf)  Eigen::Matrix<float,1,Eigen::Dynamic,Eigen::RowMajor>;
%template(Vectorcf) Eigen::Matrix<std::complex<float>,1,Eigen::Dynamic,Eigen::RowMajor>;
%template(VectorXf) Vector<float>;
%template(VectorXcf) Vector<std::complex<float>>;

%template(Arrayf)   Eigen::Array<float,1,Eigen::Dynamic,Eigen::RowMajor>;
%template(Arraycf)  Eigen::Array<std::complex<float>,1,Eigen::Dynamic,Eigen::RowMajor>;
%template(ArrayXf)  Array<float>;
%template(ArrayXcf) Array<std::complex<float>>;


%template(MatrixXf) Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
%template(MatrixXd) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
%template(MatrixXcf) Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
%template(MatrixXcd) Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

%template(VectorwiseRowOpFloat) Eigen::VectorwiseOp<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>, Eigen::Horizontal>;
//%template(VectorWiseRowOpFloat) VectorwiseRowOp<float>;

%inline %{
    // all the cool kids call it the hadamard product.
    template<typename T>
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> 
    hadamard(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & a, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &b)
    {
        return Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(a.cwiseProduct(b));
    }
    

    /// Neural Network functions

    template<typename T>
    void sigmoid(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {       
        
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.array();
        t = -t;
        m = (1 / (1 + t.exp()));
    }

    
    template<typename T>
    void sigmoidd(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    { 
        
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.array();
        m = (t * ( 1 - t ));
    }

    
    template<typename T>
    void tanH(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {   
        
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.array();  
        m = (t.tanh());
    }

    
    template<typename T>
    void tanHd(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {   
        
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.array();
        m = (1 - (t*t));
    }

    template<typename T>
    void relu(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {                  
        m = m.cwiseMax(0).eval();    
    }

    template<typename T>
    void relud(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {          
        for(size_t i = 0; i < m.rows(); i++)
            for(size_t j = 0; j < m.cols(); j++)
            {
                T x = m(i,j);
                if(x > FLT_MIN) m(i,j) = 1;
                else m(i,j) = 0;
            }
    }    

    
    template<typename T>
    void softmax(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m)
    {        
        T summed = m.sum();
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> t = m.array();        
        m = (t.exp() / summed);
    }
    
    
    // Identity or Linear
    template<typename T>
    void noalias(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & m, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & rhs)  
    {
        m.noalias() = rhs;
    }

    template<typename T>
    void mish_activate(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & Z, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> & A)
    {
        // h(x) = tanh(softplus(x)) = (1 + exp(x))^2 - 1
        //                            ------------------
        //                            (1 + exp(x))^2 + 1
        // Let s = exp(-abs(x)), t = 1 + s
        // If x >= 0, then h(x) = (t^2 - s^2) / (t^2 + s^2)
        // If x <= 0, then h(x) = (t^2 - 1) / (t^2 + 1)
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> S;
        S = ((-Z.array().abs()).exp()());
        A.array() = (S.array() + (T)1).square();
        S.noalias() = (Z.array() >= (T)0).select(S.cwiseAbs2(),(T)1);
        A.array() = (A.array() - S.array()) /
                                (A.array() + S.array());
        A.array() *= Z.array();                                
    }

    template<typename T>
    void mish_apply_jacobian(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& Z, const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& A,
                                            const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& F, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& G)
    {
        // Let h(x) = tanh(softplus(x))
        // Mish'(x) = h(x) + x * h'(x)
        // h'(x) = tanh'(softplus(x)) * softplus'(x)
        //       = [1 - h(x)^2] * exp(x) / (1 + exp(x))
        //       = [1 - h(x)^2] / (1 + exp(-x))
        // Mish'(x) = h(x) + [x - Mish(x) * h(x)] / (1 + exp(-x))
        // A = Mish(Z) = Z .* h(Z) => h(Z) = A ./ Z, h(0) = 0.6
        G.noalias() = (Z.array() == (T)0).select((T)0.6, A.cwiseQuotient(Z));
        G.array() += (Z.array() - A.array() * G.array()) / ((T)1 + (-Z).array().exp());
        G.array() *= F.array();
    }    
%}

%template(hadamard_float) hadamard<float>;
%template(sigmoid_float) sigmoid<float>;
%template(sigmoid_deriv_float) sigmoidd<float>;
%template(tanh_float) tanH<float>;
%template(tanh_deriv_float) tanHd<float>;
%template(relu_float) relu<float>;
%template(relu_deriv_float) relud<float>;
%template(softmax_float) softmax<float>;
