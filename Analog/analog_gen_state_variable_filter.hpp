// http://dafx.de/paper-archive/2020/proceedings/papers/DAFx2020_paper_70.pdf
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "GenericSoundObject.hpp"

/**
 * Generalized State-Variable Filter structure,
 * as defined by Werner and McClellan (https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_70.pdf)
 * 
 * I've added some nonlinear elements of my own design.
 */
template<typename DSP>
class GeneralSVF : public GSSoundProcessor<DSP> {
public:
    GeneralSVF();

    void reset();
    void calcCoefs(DSP r, DSP k, DSP wc);
    void setDrive (DSP newDrive);

    inline DSP process(DSP x) {
        Eigen::Matrix<DSP, 4, 1> v = A_tilde * v_n1 + B_tilde * x;
        DSP y = (C_tilde * v_n1)(0, 0) + D_tilde(0, 0) * x;
        nonlinearity(v);
        v_n1 = v;
        return y * makeup;
    }
	void ProcessSIMD(size_t n, DSP * in, DSP * out)
	{
		#pragma omp simd aligned(in,out)
		for(size_t i = 0; i < n; i++)
		{
			const DSP x = in[i];
	        Eigen::Matrix<DSP, 4, 1> v = A_tilde * v_n1 + B_tilde * x;
			DSP y = (C_tilde * v_n1)(0, 0) + D_tilde(0, 0) * x;
			nonlinearity(v);
			v_n1 = v;
			out[i] = y * makeup;
		}
	}
	void ProcessBlock(size_t n, DSP * in, DSP * out) {
		ProcessSIMD(n,in,out);
	}
	void ProcessInplace(size_t n, DSP * out) {
		ProcessSIMD(n,out);
	}
	
    inline void nonlinearity(Eigen::Matrix<DSP, 4, 1>& v) {
        v(0,0) = std::tanh(v(0,0) * drive) * invDrive;
        v(2,0) = std::tanh(v(2,0) * drive) * invDrive;
        v(3,0) = std::tanh(v(3,0) * drive) * invDrive;
    }

private:
    Eigen::Matrix<DSP, 4, 4> A;
    Eigen::Matrix<DSP, 4, 1> B;
    Eigen::Matrix<DSP, 1, 4> C;

    Eigen::Matrix<DSP, 4, 4> A_tilde;
    Eigen::Matrix<DSP, 4, 1> B_tilde;
    Eigen::Matrix<DSP, 1, 4> C_tilde;
    Eigen::Matrix<DSP, 1, 1> D_tilde;

    DSP g;
    Eigen::Matrix<DSP, 4, 1> v_n1;

    DSP drive = 1.0f;
    DSP invDrive = 1.0f;
    DSP makeup = 1.0f;
};


template<typename DSP>
GeneralSVF<DSP>::GeneralSVF() {
    A << 0.0f, 1.0f, 0.0f, 0.0f,
        -1.0f, 0.0f, 0.0f, 0.0f,
         0.0f,-1.0f, 0.0f, 1.0f,
         0.0f, 0.0f,-1.0f, 0.0f;

    B << 1.0f, 0.0f, 0.0f, 0.0f;
    C << 0.0f, 0.0f, 0.0f,-1.0f;
}

template<typename DSP>
void GeneralSVF<DSP>::reset() {
    v_n1 = Eigen::Matrix<DSP, 4, 1>::Zero();
}

template<typename DSP>
void GeneralSVF<DSP>::setDrive (DSP newDrive) {
    drive = newDrive;
    invDrive = 1.0f / drive;
    makeup = std::max(1.0f, (DSP) std::pow(drive, 0.75f));
}

template<typename DSP>
void GeneralSVF<DSP>::calcCoefs(DSP r, DSP k, DSP wc) {
    // cook parameters
    g = std::tan(wc);
    A(0, 0) = -2.0f * r;
    A(0, 3) = 4.0f * k * r * r;
    A(2, 2) = A(0, 0);

    // cook discrete state-space matrices
    Eigen::Matrix<DSP, 4, 4> I = Eigen::Matrix<DSP, 4, 4>::Identity();
    Eigen::Matrix<DSP, 4, 4> H = g * (I - g * A).inverse();
    A_tilde = 2 * H * A + I;
    B_tilde = 2 * H * B;
    C_tilde = C * (H * A + I);
    D_tilde = C * B;
}
