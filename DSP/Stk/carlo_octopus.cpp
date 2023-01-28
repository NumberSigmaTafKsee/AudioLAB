#include "Octopus.hpp"

octave::interpreter interpreter;
Octopus oct;

// fft
// stft
//

Casino::complex_vector<double> rfft(Casino::sample_vector<double> & input)
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
  Casino::complex_vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

Casino::complex_vector<double> fft(Casino::complex_vector<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i)  = (input[i]);         
  }
  olist(0) = c;
  octave_value_list output = oct.eval("fft",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  Casino::complex_vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

Casino::sample_vector<double> rifft(Casino::complex_vector<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i) = input[i];    
  }
  olist(0) = c;
  octave_value_list output = oct.eval("ifft",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  Casino::sample_vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}

Casino::complex_vector<double> ifft(Casino::complex_vector<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.size(); i++) {
    c(i) = input[i];    
  }
  olist(0) = c;
  octave_value_list output = oct.eval("ifft",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  Casino::sample_vector<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.cols(); i++)
  {    
    r[i] = x1(i);        
  }
  return r;
}


Casino::complex_matrix<double> fft2(Casino::complex_matrix<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.rows(); i++) {
    for(size_t j = 0; j < input.cols(); j++)
    {
        c(i,j) = (input(i,j));             
    }
  }
  olist(0) = c;
  octave_value_list output = oct.eval("fft2",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  Casino::complex_matrix<double> r(x1.rows(),x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.rows(); i++)
  for(size_t j = 0; j < x1.cols(); j++)
  {    
    r(i,j) = x1(i,j);        
  }
  return r;
}

Casino::complex_matrix<double> rfft2(Casino::sample_matrix<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.rows(); i++) {
    for(size_t j = 0; j < input.cols(); j++)
    {
        c(i,j).real(input(i,j));     
        c(i,j).imag(0);
    }
  }
  olist(0) = c;
  octave_value_list output = oct.eval("fft2",olist,1);    
  ComplexMatrix x1 = output(0).complex_matrix_value();    
  Casino::complex_matrix<double> r(x1.rows(),x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  
  for(size_t i = 0; i < x1.rows(); i++)
  for(size_t j = 0; j < x1.cols(); j++)
  {    
    r(i,j) = x1(i,j);        
  }
  return r;
}

Casino::sample_matrix<double> rifft2(Casino::complex_matrix<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.rows(); i++) {
    for(size_t j = 0; j < input.cols(); j++) {
        c(i,j) = input(i,j);    
  }
  olist(0) = c;
  octave_value_list output = oct.eval("ifft2",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  Casino::sample_matrix<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  

  for(size_t i = 0; i < x1.rows(); i++)
  for(size_t j = 0; j < x1.cols(); j++)
  {    
    r(i,j) = x1(i,j);        
  }
  return r;
}

Casino::complex_matrix<double> cifft2(Casino::complex_matrix<double> & input)
{
  octave_value_list olist;
  ComplexMatrix c = ComplexMatrix(1,input.size());      
  for(size_t i = 0; i < input.rows(); i++) {
    for(size_t j = 0; j < input.cols(); j++) {
        c(i,j) = input(i,j);    
  }
  olist(0) = c;
  octave_value_list output = oct.eval("ifft2",olist,1);    
  Matrix x1 = output(0).matrix_value();    
  Casino::complex_matrix<double> r(x1.cols());
  //std::cout << x1.rows() << "," << x1.cols() << std::endl;  

  for(size_t i = 0; i < x1.rows(); i++)
  for(size_t j = 0; j < x1.cols(); j++)
  {    
    r(i,j) = x1(i,j);        
  }
  return r;
}
