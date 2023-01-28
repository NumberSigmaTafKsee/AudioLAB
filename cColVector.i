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

class ComplexColumnVector : public MArray<Complex>
{
  friend class ComplexMatrix;
  friend class ComplexRowVector;

public:

  ComplexColumnVector (void);
  explicit ComplexColumnVector (octave_idx_type n);
  explicit ComplexColumnVector (const dim_vector& dv);
  ComplexColumnVector (octave_idx_type n, const Complex& val);
  ComplexColumnVector (const ComplexColumnVector& a);
  ComplexColumnVector (const MArray<Complex>& a);
  ComplexColumnVector (const Array<Complex>& a);
  explicit ComplexColumnVector (const ColumnVector& a);
  ComplexColumnVector& operator = (const ComplexColumnVector& a);

  bool operator == (const ComplexColumnVector& a) const;
  bool operator != (const ComplexColumnVector& a) const;

  // destructive insert/delete/reorder operations

  ComplexColumnVector&
  insert (const ColumnVector& a, octave_idx_type r);
  ComplexColumnVector&
  insert (const ComplexColumnVector& a, octave_idx_type r);

  ComplexColumnVector& fill (double val);
  ComplexColumnVector& fill (const Complex& val);
  ComplexColumnVector&
  fill (double val, octave_idx_type r1, octave_idx_type r2);
  ComplexColumnVector&
  fill (const Complex& val, octave_idx_type r1, octave_idx_type r2);

  ComplexColumnVector stack (const ColumnVector& a) const;
  ComplexColumnVector stack (const ComplexColumnVector& a) const;

  ComplexRowVector hermitian (void) const;
  ComplexRowVector transpose (void) const;

  friend ComplexColumnVector conj (const ComplexColumnVector& a);

  // resize is the destructive equivalent for this one

  ComplexColumnVector
  extract (octave_idx_type r1, octave_idx_type r2) const;

  ComplexColumnVector
  extract_n (octave_idx_type r1, octave_idx_type n) const;

  // column vector by column vector -> column vector operations

  ComplexColumnVector& operator += (const ColumnVector& a);
  ComplexColumnVector& operator -= (const ColumnVector& a);

  // matrix by column vector -> column vector operations

  friend ComplexColumnVector operator * (const ComplexMatrix& a,
                                                    const ColumnVector& b);

  friend ComplexColumnVector operator * (const ComplexMatrix& a,
                                                    const ComplexColumnVector& b);

  // matrix by column vector -> column vector operations

  friend ComplexColumnVector operator * (const Matrix& a,
                                                    const ComplexColumnVector& b);

  // diagonal matrix by column vector -> column vector operations

  friend ComplexColumnVector operator * (const DiagMatrix& a,
                                                    const ComplexColumnVector& b);

  friend ComplexColumnVector operator * (const ComplexDiagMatrix& a,
                                                    const ColumnVector& b);

  friend ComplexColumnVector operator * (const ComplexDiagMatrix& a,
                                                    const ComplexColumnVector& b);

  // other operations

  Complex min (void) const;
  Complex max (void) const;

  ColumnVector abs (void) const;

  // i/o
  /* wont matter in script language
  friend std::ostream& operator << (std::ostream& os,
                                               const ComplexColumnVector& a);
  friend std::istream& operator >> (std::istream& is,
                                               ComplexColumnVector& a);
  */
  void resize (octave_idx_type n, const Complex& rfv = Complex (0));
  void clear (octave_idx_type n);
};
