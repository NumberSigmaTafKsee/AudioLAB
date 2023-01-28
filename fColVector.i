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

class FloatColumnVector : public MArray<float>
{
public:

  FloatColumnVector (void) : MArray<float> (dim_vector (0, 1));

  explicit FloatColumnVector (octave_idx_type n);
  explicit FloatColumnVector (const dim_vector& dv);
  FloatColumnVector (octave_idx_type n, float val);
  FloatColumnVector (const FloatColumnVector& a);
  FloatColumnVector (const MArray<float>& a);
  FloatColumnVector (const Array<float>& a);

  FloatColumnVector& operator = (const FloatColumnVector& a);

  bool operator == (const FloatColumnVector& a) const;
  bool operator != (const FloatColumnVector& a) const;

  // destructive insert/delete/reorder operations

  FloatColumnVector&
  insert (const FloatColumnVector& a, octave_idx_type r);

  FloatColumnVector& fill (float val);
  FloatColumnVector&
  fill (float val, octave_idx_type r1, octave_idx_type r2);

  FloatColumnVector stack (const FloatColumnVector& a) const;

  FloatRowVector transpose (void) const;

  //friend FloatColumnVector real (const FloatComplexColumnVector& a);
  //friend FloatColumnVector imag (const FloatComplexColumnVector& a);

  // resize is the destructive equivalent for this one

  FloatColumnVector
  extract (octave_idx_type r1, octave_idx_type r2) const;

  FloatColumnVector
  extract_n (octave_idx_type r1, octave_idx_type n) const;

  // matrix by column vector -> column vector operations

  friend FloatColumnVector
  operator * (const FloatMatrix& a, const FloatColumnVector& b);

  // diagonal matrix by column vector -> column vector operations

  friend FloatColumnVector
  operator * (const FloatDiagMatrix& a, const FloatColumnVector& b);

  // other operations

  float min (void) const;
  float max (void) const;

  FloatColumnVector abs (void) const;

  // i/o
 /*
  friend std::ostream&
  operator << (std::ostream& os, const FloatColumnVector& a);
  friend std::istream&
  operator >> (std::istream& is, FloatColumnVector& a);
*/
  void resize (octave_idx_type n, const float& rfv = 0);
  void clear (octave_idx_type n);
};