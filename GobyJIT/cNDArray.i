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


class ComplexNDArray : public MArray<Complex>
{
public:

  ComplexNDArray (void);
  ComplexNDArray (const dim_vector& dv);
  ComplexNDArray (const dim_vector& dv, const Complex& val);
  ComplexNDArray (const ComplexNDArray& a) ;
  /*
  template <typename U>
  ComplexNDArray (const MArray<U>& a) : MArray<Complex> (a) { }

  template <typename U>
  ComplexNDArray (const Array<U>& a) : MArray<Complex> (a) { }
  */

  ComplexNDArray (const charNDArray&);

  ComplexNDArray& operator = (const ComplexNDArray& a);

  // unary operations
  boolNDArray operator ! (void) const;

  // FIXME: this is not quite the right thing.

  bool any_element_is_nan (void) const;
  bool any_element_is_inf_or_nan (void) const;
  bool all_elements_are_real (void) const;
  bool all_integers (double& max_val, double& min_val) const;
  bool too_large_for_float (void) const;

  boolNDArray all (int dim = -1) const;
  boolNDArray any (int dim = -1) const;

  ComplexNDArray cumprod (int dim = -1) const;
  ComplexNDArray cumsum (int dim = -1) const;
  ComplexNDArray prod (int dim = -1) const;
  ComplexNDArray sum (int dim = -1) const;
  ComplexNDArray xsum (int dim = -1) const;
  ComplexNDArray sumsq (int dim = -1) const;
  ComplexNDArray
  concat (const ComplexNDArray& rb, const Array<octave_idx_type>& ra_idx);
  ComplexNDArray
  concat (const NDArray& rb, const Array<octave_idx_type>& ra_idx);

  ComplexNDArray max (int dim = -1) const;
  ComplexNDArray
  max (Array<octave_idx_type>& index, int dim = -1) const;
  ComplexNDArray min (int dim = -1) const;
  ComplexNDArray
  min (Array<octave_idx_type>& index, int dim = -1) const;

  ComplexNDArray cummax (int dim = -1) const;
  ComplexNDArray
  cummax (Array<octave_idx_type>& index, int dim = -1) const;
  ComplexNDArray cummin (int dim = -1) const;
  ComplexNDArray
  cummin (Array<octave_idx_type>& index, int dim = -1) const;

  ComplexNDArray
  diff (octave_idx_type order = 1, int dim = -1) const;

  ComplexNDArray&
  insert (const NDArray& a, octave_idx_type r, octave_idx_type c);
  ComplexNDArray&
  insert (const ComplexNDArray& a, octave_idx_type r, octave_idx_type c);
  ComplexNDArray&
  insert (const ComplexNDArray& a, const Array<octave_idx_type>& ra_idx);

  NDArray abs (void) const;
  boolNDArray isnan (void) const;
  boolNDArray isinf (void) const;
  boolNDArray isfinite (void) const;

  friend ComplexNDArray conj (const ComplexNDArray& a);

  ComplexNDArray fourier (int dim = 1) const;
  ComplexNDArray ifourier (int dim = 1) const;

  ComplexNDArray fourier2d (void) const;
  ComplexNDArray ifourier2d (void) const;

  ComplexNDArray fourierNd (void) const;
  ComplexNDArray ifourierNd (void) const;

  ComplexNDArray squeeze (void) const;

  static void
  increment_index (Array<octave_idx_type>& ra_idx,
                   const dim_vector& dimensions,
                   int start_dimension = 0);

  static octave_idx_type
  compute_index (Array<octave_idx_type>& ra_idx,
                 const dim_vector& dimensions);

  // i/o
    /*
  friend std::ostream& operator << (std::ostream& os,
                                               const ComplexNDArray& a);
  friend std::istream& operator >> (std::istream& is,
                                               ComplexNDArray& a);
    */
  //  bool all_elements_are_real (void) const;
  //  bool all_integers (double& max_val, double& min_val) const;

  ComplexNDArray diag (octave_idx_type k = 0) const;

  ComplexNDArray diag (octave_idx_type m, octave_idx_type n) const;

  ComplexNDArray& changesign (void);

};
