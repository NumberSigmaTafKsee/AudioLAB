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

class FloatRowVector : public MArray<float>
{
public:

  FloatRowVector (void);
  explicit FloatRowVector (octave_idx_type n);
  explicit FloatRowVector (const dim_vector& dv);
  FloatRowVector (octave_idx_type n, float val);
  FloatRowVector (const FloatRowVector& a);
  FloatRowVector (const MArray<float>& a);
  FloatRowVector (const Array<float>& a);

  FloatRowVector& operator = (const FloatRowVector& a);

  bool operator == (const FloatRowVector& a) const;
  bool operator != (const FloatRowVector& a) const;

  // destructive insert/delete/reorder operations

  FloatRowVector&
  insert (const FloatRowVector& a, octave_idx_type c);

  FloatRowVector& fill (float val);
  FloatRowVector&
  fill (float val, octave_idx_type c1, octave_idx_type c2);

  FloatRowVector append (const FloatRowVector& a) const;

  FloatColumnVector transpose (void) const;

  //friend FloatRowVector real (const FloatComplexRowVector& a);
  //friend FloatRowVector imag (const FloatComplexRowVector& a);

  // resize is the destructive equivalent for this one

  FloatRowVector
  extract (octave_idx_type c1, octave_idx_type c2) const;

  FloatRowVector
  extract_n (octave_idx_type c1, octave_idx_type n) const;

  // row vector by matrix -> row vector

  friend FloatRowVector
  operator * (const FloatRowVector& a, const FloatMatrix& b);

  // other operations

  float min (void) const;
  float max (void) const;

  // i/o
  /*
  friend std::ostream&
  operator << (std::ostream& os, const FloatRowVector& a);
  friend std::istream&
  operator >> (std::istream& is, FloatRowVector& a);
  */
  void resize (octave_idx_type n, const float& rfv = 0);
  void clear (octave_idx_type n);
};
