#pragma once

template<typename T>
struct CVec
{
  uint64_t length;  /**< length of buffer = (requested length)/2 + 1 */
  std::vector<T> norm;          /**< norm array of size ::cvec_t.length */
  std::vector<T> phas;          /**< phase array of size ::cvec_t.length */

  CVec(uint64_t len) {      
    length = len;
    norm.resize(len);
    phas.resize(len);
  }
  ~CVec()
  {
    
  }

  CVec& operator = (const CVec<T> & c) {
    norm = c.norm;
    phas = c.phase;
    length = c.length;
    return *this;
  }

  std::pair<T&,T&> operator()(size_t i) {
    return std::pair<T&,T&>(norm[i],phas[i]);
  }
  std::pair<T,T> operator()(size_t i) const {
    return std::pair<T,T>(norm[i],phas[i]);
  }

  void set_norm(T data, size_t pos) {
    norm[pos] = data;
  }
  void set_phase(T data, size_t pos) {
    phas[pos] = data;
  }
  T& get_norm(size_t pos) {
    return norm[pos];
  }
  T& get_phase(size_t pos) {
    return phas[pos];
  }
  T* get_norm() { return norm; }
  T* get_phase() { return phas; }

  void fill_norm(T v) {
    #pragma omp simd
    for(size_t i = 0; i < length; i++) norm[i] = v;
  }
  void fill_phase(T v) {
    #pragma omp simd
    for(size_t i = 0; i < length; i++) phas[i] = v;
  }
  void zero_norm() { fill_norm((T)0); }
  void zero_phase() { fill_phase((T)0); }
  void ones_norm() { fill_norm((T)1); }
  void ones_phase() { fill_phase((T)1); }
  void zeros() {
    zero_norm();
    zero_phase();
  }

  void logmag(T lambda)
  {
    #if defined(HAVE_INTEL_IPP)
      aubio_ippsMulC(norm, lambda, norm, (int)length);
      aubio_ippsAddC(norm, 1.0, norm, (int)length);
      aubio_ippsLn(norm, norm, (int)length);
    #else      
      #pragma omp simd
      for (size_t j=0; j< s->length; j++) {
        s->norm[j] = std::log(lambda * norm[j] + 1);
      }
    #endif
  }

  void print() {
    std::cout << "norm: ";
    for(size_t i = 0; i < length; i++) std::cout << norm[i] << ",";
    std::cout << std::endl;
    std::cout << "phas: ";
    for(size_t i = 0; i < length; i++) std::cout << phas[i] << ",";
    std::cout << std::endl;
  }


};

std::ostream& operator << (std::ostream &o, const CVec<T> & v)
{
    o << "norm: ";
    for(size_t i = 0; i < length; i++) o << norm[i] << ",";
    o << std::endl;
    o << "phas: ";
    for(size_t i = 0; i < length; i++) o << phas[i] << ",";
    o << std::endl;
    return o;
}