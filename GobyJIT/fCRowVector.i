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

class FloatComplexRowVector : public MArray<FloatComplex>
{
  friend class FloatComplexColumnVector;

public:

  FloatComplexRowVector (void);
  explicit FloatComplexRowVector (octave_idx_type n);
  explicit FloatComplexRowVector (const dim_vector& dv);
  FloatComplexRowVector (octave_idx_type n, const FloatComplex& val);
  FloatComplexRowVector (const FloatComplexRowVector& a);
  FloatComplexRowVector (const MArray<FloatComplex>& a);
  FloatComplexRowVector (const Array<FloatComplex>& a);
  explicit FloatComplexRowVector (const FloatRowVector& a);
  FloatComplexRowVector& operator = (const FloatComplexRowVector& a);

  bool operator == (const FloatComplexRowVector& a) const;
  bool operator != (const FloatComplexRowVector& a) const;

  // destructive insert/delete/reorder operations

  FloatComplexRowVector&
  insert (const FloatRowVector& a, octave_idx_type c);
  FloatComplexRowVector&
  insert (const FloatComplexRowVector& a, octave_idx_type c);

  FloatComplexRowVector& fill (float val);
  FloatComplexRowVector& fill (const FloatComplex& val);
  FloatComplexRowVector&
  fill (float val, octave_idx_type c1, octave_idx_type c2);
  FloatComplexRowVector&
  fill (const FloatComplex& val, octave_idx_type c1, octave_idx_type c2);

  FloatComplexRowVector append (const FloatRowVector& a) const;
  FloatComplexRowVector
  append (const FloatComplexRowVector& a) const;

  FloatComplexColumnVector hermitian (void) const;
  FloatComplexColumnVector transpose (void) const;

  friend FloatComplexRowVector
  conj (const FloatComplexRowVector& a);

  // resize is the destructive equivalent for this one

  FloatComplexRowVector
  extract (octave_idx_type c1, octave_idx_type c2) const;

  FloatComplexRowVector
  extract_n (octave_idx_type c1, octave_idx_type n) const;

  // row vector by row vector -> row vector operations

  FloatComplexRowVector& operator += (const FloatRowVector& a);
  FloatComplexRowVector& operator -= (const FloatRowVector& a);

  // row vector by matrix -> row vector

  friend FloatComplexRowVector
  operator * (const FloatComplexRowVector& a, const FloatComplexMatrix& b);

  friend FloatComplexRowVector
  operator * (const FloatRowVector& a, const FloatComplexMatrix& b);

  // other operations

  FloatComplex min (void) const;
  FloatComplex max (void) const;

  // i/o
    /*
  friend std::ostream&
  operator << (std::ostream& os, const FloatComplexRowVector& a);
  friend std::istream&
  operator >> (std::istream& is, FloatComplexRowVector& a);
    */
  void resize (octave_idx_type n, const FloatComplex& rfv = FloatComplex (0));
  void clear (octave_idx_type n);

};