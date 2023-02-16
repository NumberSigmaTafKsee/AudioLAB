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

class ComplexRowVector : public MArray<Complex>
{
  friend class ComplexColumnVector;

public:

  ComplexRowVector (void);
  explicit ComplexRowVector (octave_idx_type n);
  explicit ComplexRowVector (const dim_vector& dv);
  ComplexRowVector (octave_idx_type n, const Complex& val);
  ComplexRowVector (const ComplexRowVector& a);
  ComplexRowVector (const MArray<Complex>& a);
  ComplexRowVector (const Array<Complex>& a);
  explicit ComplexRowVector (const RowVector& a);
  ComplexRowVector& operator = (const ComplexRowVector& a);

  bool operator == (const ComplexRowVector& a) const;
  bool operator != (const ComplexRowVector& a) const;

  // destructive insert/delete/reorder operations

  ComplexRowVector&
  insert (const RowVector& a, octave_idx_type c);
  ComplexRowVector&
  insert (const ComplexRowVector& a, octave_idx_type c);

  ComplexRowVector& fill (double val);
  ComplexRowVector& fill (const Complex& val);
  ComplexRowVector&
  fill (double val, octave_idx_type c1, octave_idx_type c2);
  ComplexRowVector&
  fill (const Complex& val, octave_idx_type c1, octave_idx_type c2);

  ComplexRowVector append (const RowVector& a) const;
  ComplexRowVector append (const ComplexRowVector& a) const;

  ComplexColumnVector hermitian (void) const;
  ComplexColumnVector transpose (void) const;

  friend ComplexRowVector conj (const ComplexRowVector& a);

  // resize is the destructive equivalent for this one

  ComplexRowVector
  extract (octave_idx_type c1, octave_idx_type c2) const;

  ComplexRowVector
  extract_n (octave_idx_type c1, octave_idx_type n) const;

  // row vector by row vector -> row vector operations

  ComplexRowVector& operator += (const RowVector& a);
  ComplexRowVector& operator -= (const RowVector& a);

  // row vector by matrix -> row vector

  friend ComplexRowVector
  operator * (const ComplexRowVector& a, const ComplexMatrix& b);

  friend ComplexRowVector
  operator * (const RowVector& a, const ComplexMatrix& b);

  // other operations

  Complex min (void) const;
  Complex max (void) const;

  // i/o
  /*
  friend std::ostream&
  operator << (std::ostream& os, const ComplexRowVector& a);
  friend std::istream&
  operator >> (std::istream& is, ComplexRowVector& a);
  */
  void resize (octave_idx_type n, const Complex& rfv = Complex (0));

  void clear (octave_idx_type n);
};
