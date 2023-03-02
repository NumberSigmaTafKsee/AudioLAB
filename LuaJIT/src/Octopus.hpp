#pragma once

#include <iostream>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>




// this is not for calculation
// it is meant to marshal data from Octave
// the calculations are done by Octave
// just need a convenient way to return vectors and matrix 
// a scalar is just one column (1,1)
template<typename T>
struct OctopusMatrix : public std::vector<T>
{
    size_t rows,cols;
    OctopusMatrix() = default;
    OctopusMatrix(size_t m, size_t n) : std::vector<T>(m*n)
    {
        rows = m;
        cols = n;
    }
    T operator()(int i, int j) const { return (*this)[i*cols + j]; }
    T& operator()(int i, int j) { return (*this)[i*cols + j]; }

    std::vector<T> get_row(int i) {
        std::vector<T> r(cols);
        memcpy(r.data(),this->data()+i*cols,cols*sizeof(T));
        return r;
    }        
    std::vector<T> get_col(int j) {
        std::vector<T> r(rows);
        for(int i = 0; i < rows; i++)
            r[i] = (*this)(i,j);
        return r;
    }
    void resize(size_t i, size_t j) {
      std::vector<T>::resize(i*j);
      rows = i;
      cols = j;
    }
    OctopusMatrix<T>& operator = (const OctopusMatrix<T>& m) {      
      rows = m.rows;
      cols = m.cols;
      if(m.size() > 0)
      {
        resize(rows,cols);
        std::copy(m.begin(),m.end(),std::vector<T>::begin());      
      }
      return *this;
    }
};

struct OctopusMatrixType {
  OctopusMatrix<double> md;
  OctopusMatrix<std::complex<double>> mc;
  enum {
    DOUBLE,
    COMPLEX,
  };
  int type = DOUBLE;

  
  OctopusMatrixType& operator = (const OctopusMatrixType & t) {
    type = t.type;
    md   = t.md;
    mc   = t.mc;
    return *this;
  }
};

octave::interpreter interpreter;

struct Octopus
{   
    Octopus() {
      interpreter.initialize_history(false);
      interpreter.initialize();
      interpreter.execute();
      std::string path = ".";
      octave_value_list p;
      p(0) = path;
      octave_value_list o1 = interpreter.feval("addpath", p, 1);            
      run_script("startup.m");
    }
    ~Octopus()
    {
        
    }
    
    void run_script(const std::string& s) {
      octave::source_file(s);
    }
    
    octave_value_list eval(std::string func, octave_value_list inputs, int noutputs=1)
    {          
      octave_value_list out =interpreter.feval(func.c_str(), inputs, noutputs);
      return out;
    }

    octave_value_list operator()(std::string func, octave_value_list inputs, int noutputs=1)
    {
      return eval(func,inputs,noutputs);
    }
    
    template<typename T>
    std::list<OctopusMatrixType> feval(const std::string & func, const OctopusMatrix<T>& m, int noutputs=0) {      
      octave_value_list l;
      for(int i = 0; i < m.rows; i++)
      {
        if(m.cols == 1) {
          l(i) = m(i,0);
        }
        else {
          RowVector rv(m.cols);        
          for(int j = 0; j < m.cols; j++)
          {
            rv(j) = m(i,j);
          }
          l(i) = rv;
        }
      }
      octave_value_list out =interpreter.feval(func.c_str(), l, noutputs);
      std::list<OctopusMatrixType> olist;
      for(size_t i = 0; i < out.length(); i++)
      {
        OctopusMatrix<T>  q;
        OctopusMatrix<std::complex<T>>  qc;
        OctopusMatrixType r;
        
        auto x = out(i);
        if(x.is_scalar_type()) {          
          q.resize(1,1);
          q(0,0) = x.double_value();
          r.type = OctopusMatrixType::DOUBLE;
          r.md   = q;
          olist.push_back(r);
        }        
        else if(x.is_complex_scalar()) {
          qc.resize(1,1);
          Complex c = x.complex_value();
          qc(0,0) = c.real();
          qc(0,1) = c.imag();
          r.type = OctopusMatrixType::COMPLEX;
          r.mc   = qc;
          olist.push_back(r);
        }        
        else if(x.is_matrix_type()) {          
          auto matrix = x.matrix_value();          
          q.resize(matrix.rows(),matrix.cols());
          for(size_t j = 0; j < matrix.rows(); j++)
            for(size_t k=0; k < matrix.cols(); k++)
              q(j,k) = matrix(j,k);
          r.type = OctopusMatrixType::DOUBLE;
          r.md   = q;          
          olist.push_back(r);          
        }        
        else if(x.is_complex_matrix()) {
          auto matrix = x.complex_matrix_value();          
          qc.resize(matrix.rows(),matrix.cols());          
          for(size_t j = 0; j < matrix.rows(); j++)
            for(size_t k=0; k < matrix.cols(); k++)
            {
                qc(j,k*2) = matrix(j,k).real();
                qc(j,k*2+1) = matrix(j,k).imag();
            }
            r.type = OctopusMatrixType::COMPLEX;
            r.mc   = qc;
            olist.push_back(r);
        }               
      }
      return olist;
    }
};

Octopus oct;

std::vector<std::complex<double>> fft(std::vector<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i).real(input[i]);     
    c(i).imag(0);
  }
  olist(0) = c;
  octave_value_list output = oct.eval("fft",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  std::vector<std::complex<double>> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<std::complex<double>> cfft(std::vector<std::complex<double>> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i) = input[i];
  }
  olist(0) = c;
  octave_value_list output = oct.eval("fft",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  std::vector<std::complex<double>> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> ifft(std::vector<std::complex<double>> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i) = input[i];    
  }
  olist(0) = c;
  octave_value_list output = oct.eval("ifft",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<std::complex<double>> cifft(std::vector<std::complex<double>> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i) = input[i];    
  }
  olist(0) = c;
  octave_value_list output = oct.eval("ifft",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  std::vector<std::complex<double>> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

void butterworth()
{
  octave_idx_type n = 2;
  octave_value_list in,out;
  in(0) = 4;
  in(1) = 0.5;  
  out = oct.eval("butterworth_filter",in,2);  
  std::cout << out.length() << std::endl;
  std::cout << out(0).length() << std::endl;  
  auto x = out(0).matrix_value();
  std::cout << x.rows() << "," << x.cols() << std::endl;
  std::cout << x << std::endl;
  std::cout << out(1).double_value() << std::endl;
}

std::vector<double> fftconv(std::vector<double> & x, std::vector<double> & y)
{
  octave_value_list olist;
  Matrix m1 = Matrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    m1(i) = (x[i]);    
  }
  Matrix m2 = Matrix(1,y.size());      
  for(size_t i = 0; i < y.size(); i++) {
    m2(i) = (y[i]);    
  }
  olist(0) = m1;
  olist(1) = m2;
  octave_value_list output = oct.eval("fftconv",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}
std::vector<std::complex<double>> cfftconv(std::vector<std::complex<double>> & x, std::vector<std::complex<double>> & y)
{
  octave_value_list olist;
  ComplexMatrix c1 = ComplexMatrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    c1(i) = (x[i]);    
  }
  ComplexMatrix c2 = ComplexMatrix(1,y.size());      
  for(size_t i = 0; i < y.size(); i++) {
    c2(i) = (y[i]);    
  }
  olist(0) = c1;
  olist(1) = c2;
  octave_value_list output = oct.eval("fftconv",olist,1);    
  ComplexMatrix x1 = output(0).matrix_value();    
  std::vector<std::complex<double>> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}


std::vector<double> fftfilt(std::vector<double> & b, std::vector<double> & x)
{
  octave_value_list olist;
  Matrix m1 = Matrix(1,b.size());      
  for(size_t i = 0; i < b.size(); i++) {
    m1(i) = (b[i]);    
  }
  Matrix m2 = Matrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    m2(i) = (x[i]);    
  }
  olist(0) = m1;
  olist(1) = m2;
  octave_value_list output = oct.eval("fftfilt",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}


std::vector<double> filter(std::vector<double> & b, std::vector<double> & a, std::vector<double> & x)
{
  octave_value_list olist;
  Matrix m1 = Matrix(1,b.size());      
  for(size_t i = 0; i < b.size(); i++) {
    m1(i) = (b[i]);    
  }
  Matrix m2 = Matrix(1,a.size());      
  for(size_t i = 0; i < a.size(); i++) {
    m2(i) = (a[i]);    
  }
  Matrix m3 = Matrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    m3(i) = (x[i]);    
  }
  olist(0) = m1;
  olist(1) = m2;
  olist(2) = m3;

  octave_value_list output = oct.eval("filter",olist,1);    
  
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}


std::vector<std::complex<double>> freqz(std::vector<std::complex<double>> & b, std::vector<std::complex<double>> & a)
{
  octave_value_list olist;
  ComplexMatrix m1 = Matrix(1,b.size());      
  for(size_t i = 0; i < b.size(); i++) {
    m1(i) = (b[i]);    
  }
  ComplexMatrix m2 = Matrix(1,a.size());      
  for(size_t i = 0; i < a.size(); i++) {
    m2(i) = (a[i]);    
  }
  olist(0) = m1;
  olist(1) = m2;
  octave_value_list output = oct.eval("freqz",olist,1);    
  ComplexMatrix x1 = output(0).matrix_value();    
  std::vector<std::complex<double>> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> sinc(std::vector<double> & x)
{
  octave_value_list olist;
  Matrix m1 = Matrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    m1(i) = (x[i]);    
  }  
  olist(0) = m1;  
  octave_value_list output = oct.eval("sinc",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> unwrap(std::vector<double> & x)
{
  octave_value_list olist;
  Matrix m1 = Matrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    m1(i) = (x[i]);    
  }  
  olist(0) = m1;  
  octave_value_list output = oct.eval("unwrap",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> bartlett(int m)
{
  octave_value_list olist;  
  olist(0) = m;
  octave_value_list output = oct.eval("bartlett",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> hamming(int m)
{
  octave_value_list olist;  
  olist(0) = m;
  octave_value_list output = oct.eval("hamming",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> sinetone(double freq, double rate, double sec, double ampl)
{
  octave_value_list olist;  
  olist(0) = freq;
  olist(1) = rate;
  olist(2) = sec;
  olist(3) = ampl;
  octave_value_list output = oct.eval("sinetone",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<double> sinewave(int m, double n, double d)
{
  octave_value_list olist;  
  olist(0) = m;
  olist(1) = n;
  olist(2) = d;
  octave_value_list output = oct.eval("sinewave",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

std::vector<std::complex<double>> stft(std::vector<double> & x, int win_size)
{
  octave_value_list olist;  
  Matrix m1 = Matrix(1,x.size());      
  for(size_t i = 0; i < x.size(); i++) {
    m1(i) = (x[i]);    
  }  
  olist(0) = m1;  
  olist(1) = win_size;
  octave_value_list output = oct.eval("stft",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  std::vector<std::complex<double>> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

/*
std::vector<double> synthesis(std::vector<std::complex<double>> & stft, int win_size)
{
  octave_value_list olist;  
  ComplexMatrix m1 = ComplexMatrix(1,stft.size());      
  for(size_t i = 0; i < stft.size(); i++) {
    m1(i) = (stft[i]);    
  }  
  olist(0) = m1;  
  olist(1) = win_size;
  octave_value_list output = oct.eval("synthesis",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  std::vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}
*/