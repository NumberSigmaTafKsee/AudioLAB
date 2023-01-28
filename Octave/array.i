////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1996-2022 The Octave Project Developers
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

%{
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
%}

//%include "octave/Array-fwd.h"

// not all this stuff is going to work right now or is even needed
template <typename T, typename Alloc>
class Array
{
public:

  void make_unique (void);
  typedef T element_type;
  typedef T value_type;
  //! Used for operator(), and returned by numel() and size()
  //! (beware: signed integer)
  typedef octave_idx_type size_type;
  typedef typename ref_param<T>::type crefT;

  typedef bool (*compare_fcn_type) (typename ref_param<T>::type,
                                    typename ref_param<T>::type);


  //! Empty ctor (0 by 0).
  Array (void);
  explicit Array (const dim_vector& dv);
  explicit Array (const dim_vector& dv, const T& val);
  explicit Array (T *ptr, const dim_vector& dv,
                  const Alloc& xallocator = Alloc ());
  Array (const Array<T, Alloc>& a, const dim_vector& dv);

  //! Constructor from standard library sequence containers.
  template<template <typename...> class Container>
  Array (const Container<T>& a, const dim_vector& dv);

  //! Type conversion case.
  template <typename U, typename A = Alloc>
  Array (const Array<U, A>& a);
  Array (const Array<T, Alloc>& a);
  Array (Array<T, Alloc>&& a);

  virtual ~Array (void);
  Array<T, Alloc>& operator = (const Array<T, Alloc>& a);
  Array<T, Alloc>& operator = (Array<T, Alloc>&& a);
  
  void fill (const T& val);
  void clear (void);
  void clear (const dim_vector& dv);

  void clear (octave_idx_type r, octave_idx_type c);
  octave_idx_type numel (void) const;
  Array<T, Alloc> as_column (void) const;
  Array<T, Alloc> as_row (void) const;
  Array<T, Alloc> as_matrix (void) const;
  octave_idx_type dim1 (void) const;
  octave_idx_type rows (void) const;
  octave_idx_type dim2 (void) const;
  octave_idx_type cols (void) const;
  octave_idx_type columns (void) const;  
  octave_idx_type dim3 (void) const;
  octave_idx_type pages (void) const;
  size_type size (const size_type d) const;
  std::size_t byte_size (void) const;
  const dim_vector& dims (void) const;
  Array<T, Alloc> squeeze (void) const;
  octave_idx_type compute_index (octave_idx_type i, octave_idx_type j) const;
  octave_idx_type compute_index (octave_idx_type i, octave_idx_type j,
                                 octave_idx_type k) const;
  octave_idx_type compute_index (const Array<octave_idx_type>& ra_idx) const;

  octave_idx_type compute_index_unchecked (const Array<octave_idx_type>& ra_idx) const;
  T& xelem (octave_idx_type n);
  crefT xelem (octave_idx_type n) const;

  T& xelem (octave_idx_type i, octave_idx_type j);
  crefT xelem (octave_idx_type i, octave_idx_type j) const;
  T& xelem (octave_idx_type i, octave_idx_type j, octave_idx_type k);
  crefT xelem (octave_idx_type i, octave_idx_type j, octave_idx_type k) const;
  T& xelem (const Array<octave_idx_type>& ra_idx);
  crefT xelem (const Array<octave_idx_type>& ra_idx) const;
  T& checkelem (octave_idx_type n);
  T& checkelem (octave_idx_type i, octave_idx_type j);
  T& checkelem (octave_idx_type i, octave_idx_type j, octave_idx_type k);
  T& checkelem (const Array<octave_idx_type>& ra_idx);
  T& elem (octave_idx_type n);
  T& elem (octave_idx_type i, octave_idx_type j);
  T& elem (octave_idx_type i, octave_idx_type j, octave_idx_type k);
  T& elem (const Array<octave_idx_type>& ra_idx);
  
  //T& operator () (octave_idx_type n);
  //T& operator () (octave_idx_type i, octave_idx_type j);
  //T& operator () (octave_idx_type i, octave_idx_type j, octave_idx_type k);
  //T& operator () (const Array<octave_idx_type>& ra_idx);

  %extend {
    T __getitem__(octave_idx_type n) {
        return (*$self)(n);
    }
    void __setitem__(octave_idx_type n, const T& value) {
        (*$self)(n) = value;
    }
    T& get2(octave_idx_type i, octave_idx_type j) {
        return (*$self)(i,j);
    }
    T& get3(octave_idx_type i, octave_idx_type j, octave_idx_type k) {
        return (*$self)(i,j,k);
    }
    void set2(octave_idx_type i, octave_idx_type j, const T& value) {
        (*$self)(i,j) = value;
    }
    void set3(octave_idx_type i, octave_idx_type j, octave_idx_type k, const T& value) {
        (*$self)(i,j,k) = value;
    }
    T& getN(const Array<octave_idx_type> & ra) {
        return (*$self)(ra);
    }
    void setN(const Array<octave_idx_type> & ra, const T& value) {
        (*$self)(ra) = value;
    }
  }

  crefT checkelem (octave_idx_type n) const;
  crefT checkelem (octave_idx_type i, octave_idx_type j) const;
  crefT checkelem (octave_idx_type i, octave_idx_type j,
                                octave_idx_type k) const;
  crefT checkelem (const Array<octave_idx_type>& ra_idx) const;
  crefT elem (octave_idx_type n) const;

  crefT elem (octave_idx_type i, octave_idx_type j) const;
  crefT elem (octave_idx_type i, octave_idx_type j, octave_idx_type k) const;
  crefT elem (const Array<octave_idx_type>& ra_idx) const;
  crefT operator () (octave_idx_type n) const;
  crefT operator () (octave_idx_type i, octave_idx_type j) const;
  crefT operator () (octave_idx_type i, octave_idx_type j,
                     octave_idx_type k) const;
  crefT operator () (const Array<octave_idx_type>& ra_idx) const;
  Array<T, Alloc> column (octave_idx_type k) const;  
  Array<T, Alloc> page (octave_idx_type k) const;

  Array<T, Alloc> linear_slice (octave_idx_type lo, octave_idx_type up) const;
  Array<T, Alloc> reshape (octave_idx_type nr, octave_idx_type nc) const;
  Array<T, Alloc> reshape (const dim_vector& new_dims) const;
  Array<T, Alloc> permute (const Array<octave_idx_type>& vec, bool inv = false) const;
  Array<T, Alloc> ipermute (const Array<octave_idx_type>& vec) const;
  bool issquare (void) const;
  bool isempty (void) const;
  bool isvector (void) const;
  bool is_nd_vector (void) const;
  
  
  Array<T, Alloc> transpose (void) const;
  Array<T, Alloc> hermitian (T (*fcn) (const T&) = nullptr) const;
  const T * data (void) const;

  T * fortran_vec (void);
  bool is_shared (void);
  int ndims (void) const;
  
  Array<T, Alloc> index (const octave::idx_vector& i) const;
  Array<T, Alloc> index (const octave::idx_vector& i, const octave::idx_vector& j) const;
  Array<T, Alloc> index (const Array<octave::idx_vector>& ia) const;
  
  //T resize_fill_value (void) const;
  
  void resize2 (octave_idx_type nr, octave_idx_type nc, const T& rfv);
  void resize2 (octave_idx_type nr, octave_idx_type nc);
  void resize1 (octave_idx_type n, const T& rfv);
  void resize1 (octave_idx_type n);
  void resize (const dim_vector& dv, const T& rfv);
  void resize (const dim_vector& dv);
  
  
  Array<T, Alloc> index (const octave::idx_vector& i, bool resize_ok, const T& rfv) const;
  Array<T, Alloc> index (const octave::idx_vector& i, bool resize_ok) const;
  Array<T, Alloc> index (const octave::idx_vector& i, const octave::idx_vector& j,
                               bool resize_ok,
                               const T& rfv) const;
  Array<T, Alloc> index (const octave::idx_vector& i, const octave::idx_vector& j,
                  bool resize_ok) const;
  Array<T, Alloc> index (const Array<octave::idx_vector>& ia, bool resize_ok,
                               const T& rfv) const;
  Array<T, Alloc> index (const Array<octave::idx_vector>& ia, bool resize_ok) const;
  
  void assign (const octave::idx_vector& i, const Array<T, Alloc>& rhs, const T& rfv);
  void assign (const octave::idx_vector& i, const Array<T, Alloc>& rhs);
  void assign (const octave::idx_vector& i, const octave::idx_vector& j,
                            const Array<T, Alloc>& rhs,
                            const T& rfv);
  void assign (const octave::idx_vector& i, const octave::idx_vector& j, const Array<T, Alloc>& rhs);
  void assign (const Array<octave::idx_vector>& ia, const Array<T, Alloc>& rhs, const T& rfv);
  void assign (const Array<octave::idx_vector>& ia, const Array<T, Alloc>& rhs);
  
  void delete_elements (const octave::idx_vector& i);
  void delete_elements (int dim, const octave::idx_vector& i);
  void delete_elements (const Array<octave::idx_vector>& ia);
  
  Array<T, Alloc>& insert (const Array<T, Alloc>& a, const Array<octave_idx_type>& idx);
  Array<T, Alloc>& insert (const Array<T, Alloc>& a, octave_idx_type r, octave_idx_type c);
  void maybe_economize (void);
  void print_info (std::ostream& os, const std::string& prefix) const;
  
  Array<T, Alloc> sort (int dim = 0, sortmode mode = ASCENDING) const;
  Array<T, Alloc> sort (Array<octave_idx_type>& sidx, int dim = 0,
                                     sortmode mode = ASCENDING) const;

  //! Ordering is auto-detected or can be specified.
  sortmode issorted (sortmode mode = UNSORTED) const;

  //! Sort by rows returns only indices.
  Array<octave_idx_type> sort_rows_idx (sortmode mode = ASCENDING) const;

  //! Ordering is auto-detected or can be specified.
  sortmode is_sorted_rows (sortmode mode = UNSORTED) const;

  //! Do a binary lookup in a sorted array.  Must not contain NaNs.
  //! Mode can be specified or is auto-detected by comparing 1st and last element.
  octave_idx_type lookup (const T& value, sortmode mode = UNSORTED) const;

  //! Ditto, but for an array of values, specializing on the case when values
  //! are sorted.  NaNs get the value N.
  Array<octave_idx_type> lookup (const Array<T, Alloc>& values,
                                                sortmode mode = UNSORTED) const;

  //! Count nonzero elements.
  octave_idx_type nnz (void) const;

  //! Find indices of (at most n) nonzero elements.  If n is specified,
  //! backward specifies search from backward.
  Array<octave_idx_type> find (octave_idx_type n = -1,
                                              bool backward = false) const;

  //! Returns the n-th element in increasing order, using the same
  //! ordering as used for sort.  n can either be a scalar index or a
  //! contiguous range.
  Array<T, Alloc> nth_element (const octave::idx_vector& n, int dim = 0) const;

  //! Get the kth super or subdiagonal.  The zeroth diagonal is the
  //! ordinary diagonal.
  Array<T, Alloc> diag (octave_idx_type k = 0) const;

  Array<T, Alloc> diag (octave_idx_type m, octave_idx_type n) const;

  //! Concatenation along a specified (0-based) dimension, equivalent
  //! to cat().  dim = -1 corresponds to dim = 0 and dim = -2
  //! corresponds to dim = 1, but apply the looser matching rules of
  //! vertcat/horzcat.
  static Array<T, Alloc>
  cat (int dim, octave_idx_type n, const Array<T, Alloc> *array_list);
    /*
  //! Apply function fcn to each element of the Array<T, Alloc>.  This function
  //! is optimized with a manually unrolled loop.
#if defined (OCTAVE_HAVE_STD_PMR_POLYMORPHIC_ALLOCATOR)
  template <typename U, typename F,
            typename A = std::pmr::polymorphic_allocator<U>>
#else
  template <typename U, typename F, typename A = std::allocator<U>>
#endif
  Array<U, A>
  map (F fcn) const
  {
    octave_idx_type len = numel ();

    const T *m = data ();

    Array<U, A> result (dims ());
    U *p = result.fortran_vec ();

    octave_idx_type i;
    for (i = 0; i < len - 3; i += 4)
      {
        octave_quit ();

        p[i] = fcn (m[i]);
        p[i+1] = fcn (m[i+1]);
        p[i+2] = fcn (m[i+2]);
        p[i+3] = fcn (m[i+3]);
      }

    octave_quit ();

    for (; i < len; i++)
      p[i] = fcn (m[i]);

    return result;
  }

  //@{
  //! Overloads for function references.
#if defined (OCTAVE_HAVE_STD_PMR_POLYMORPHIC_ALLOCATOR)
  template <typename U, typename A = std::pmr::polymorphic_allocator<U>>
#else
  template <typename U, typename A = std::allocator<U>>
#endif
  Array<U, A>
  map (U (&fcn) (T)) const
  { return map<U, U (&) (T), A> (fcn); }

#if defined (OCTAVE_HAVE_STD_PMR_POLYMORPHIC_ALLOCATOR)
  template <typename U, typename A = std::pmr::polymorphic_allocator<U>>
#else
  template <typename U, typename A = std::allocator<U>>
#endif
  Array<U, A>
  map (U (&fcn) (const T&)) const
  { return map<U, U (&) (const T&), A> (fcn); }
  //@}

  //! Generic any/all test functionality with arbitrary predicate.
  template <typename F, bool zero>
  bool test (F fcn) const
  {
    return octave::any_all_test<F, T, zero> (fcn, data (), numel ());
  }

  //@{
  //! Simpler calls.
  template <typename F>
  bool test_any (F fcn) const
  { return test<F, false> (fcn); }

  template <typename F>
  bool test_all (F fcn) const
  { return test<F, true> (fcn); }
  //@}
  */
  //@{
  //! Overloads for function references.
  bool test_any (bool (&fcn) (T)) const;
  bool test_any (bool (&fcn) (const T&)) const;
  bool test_all (bool (&fcn) (T)) const;
  bool test_all (bool (&fcn) (const T&)) const;
  template <typename U, typename A> friend class Array;
  bool optimize_dimensions (const dim_vector& dv);
};

%ignore rec_permute_helper;
%ignore rec_index_helper;
%ignore rec_resize_helper;
%include <octave/Array.cc>
