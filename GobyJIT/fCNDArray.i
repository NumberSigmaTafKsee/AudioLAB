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

class FloatComplexNDArray : public MArray<FloatComplex>
{
public:

  FloatComplexNDArray (void) : MArray<FloatComplex> ();
  FloatComplexNDArray (const dim_vector& dv);
  FloatComplexNDArray (const dim_vector& dv, const FloatComplex& val);
  FloatComplexNDArray (const FloatComplexNDArray& a);
  /*
  template <typename U>
  FloatComplexNDArray (const MArray<U>& a) : MArray<FloatComplex> (a) { }

  template <typename U>
  FloatComplexNDArray (const Array<U>& a) : MArray<FloatComplex> (a) { }
    */
  FloatComplexNDArray (const charNDArray&);

  FloatComplexNDArray& operator = (const FloatComplexNDArray& a);

  // unary operations
  boolNDArray operator ! (void) const;

  // FIXME: this is not quite the right thing.

  bool any_element_is_nan (void) const;
  bool any_element_is_inf_or_nan (void) const;
  bool all_elements_are_real (void) const;
  bool all_integers (float& max_val, float& min_val) const;
  bool too_large_for_float (void) const;

  boolNDArray all (int dim = -1) const;
  boolNDArray any (int dim = -1) const;

  FloatComplexNDArray cumprod (int dim = -1) const;
  FloatComplexNDArray cumsum (int dim = -1) const;
  FloatComplexNDArray prod (int dim = -1) const;
  ComplexNDArray dprod (int dim = -1) const;
  FloatComplexNDArray sum (int dim = -1) const;
  ComplexNDArray dsum (int dim = -1) const;
  FloatComplexNDArray sumsq (int dim = -1) const;
  FloatComplexNDArray
  concat (const FloatComplexNDArray& rb, const Array<octave_idx_type>& ra_idx);
  FloatComplexNDArray
  concat (const FloatNDArray& rb, const Array<octave_idx_type>& ra_idx);

  FloatComplexNDArray max (int dim = -1) const;
  FloatComplexNDArray
  max (Array<octave_idx_type>& index, int dim = -1) const;
  FloatComplexNDArray min (int dim = -1) const;
  FloatComplexNDArray
  min (Array<octave_idx_type>& index, int dim = -1) const;

  FloatComplexNDArray cummax (int dim = -1) const;
  FloatComplexNDArray
  cummax (Array<octave_idx_type>& index, int dim = -1) const;
  FloatComplexNDArray cummin (int dim = -1) const;
  FloatComplexNDArray
  cummin (Array<octave_idx_type>& index, int dim = -1) const;

  FloatComplexNDArray
  diff (octave_idx_type order = 1, int dim = -1) const;

  FloatComplexNDArray&
  insert (const NDArray& a, octave_idx_type r, octave_idx_type c);
  FloatComplexNDArray&
  insert (const FloatComplexNDArray& a, octave_idx_type r, octave_idx_type c);
  FloatComplexNDArray&
  insert (const FloatComplexNDArray& a, const Array<octave_idx_type>& ra_idx);

  FloatNDArray abs (void) const;
  boolNDArray isnan (void) const;
  boolNDArray isinf (void) const;
  boolNDArray isfinite (void) const;

  friend FloatComplexNDArray conj (const FloatComplexNDArray& a);

  FloatComplexNDArray fourier (int dim = 1) const;
  FloatComplexNDArray ifourier (int dim = 1) const;

  FloatComplexNDArray fourier2d (void) const;
  FloatComplexNDArray ifourier2d (void) const;

  FloatComplexNDArray fourierNd (void) const;
  FloatComplexNDArray ifourierNd (void) const;

  FloatComplexNDArray squeeze (void) const;
  
  static void
  increment_index (Array<octave_idx_type>& ra_idx,
                   const dim_vector& dimensions, int start_dimension = 0);

  static octave_idx_type
  compute_index (Array<octave_idx_type>& ra_idx, const dim_vector& dimensions);

  // i/o
 /*
  friend std::ostream& operator << (std::ostream& os,
                                               const FloatComplexNDArray& a);
  friend std::istream& operator >> (std::istream& is,
                                               FloatComplexNDArray& a);
*/
  //  bool all_elements_are_real (void) const;
  //  bool all_integers (float& max_val, float& min_val) const;

  FloatComplexNDArray diag (octave_idx_type k = 0) const;

  FloatComplexNDArray
  diag (octave_idx_type m, octave_idx_type n) const;

  FloatComplexNDArray& changesign (void);

};