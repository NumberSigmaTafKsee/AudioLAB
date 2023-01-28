////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1993-2022 The Octave Project Developers
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


class ColumnVector : public MArray<double>
{
public:

  ColumnVector (void);

  explicit ColumnVector (octave_idx_type n);
  explicit ColumnVector (const dim_vector& dv);

  ColumnVector (octave_idx_type n, double val);
  ColumnVector (const ColumnVector& a);
  ColumnVector (const MArray<double>& a);
  ColumnVector (const Array<double>& a);

  ColumnVector& operator = (const ColumnVector& a);

  bool operator == (const ColumnVector& a) const;
  bool operator != (const ColumnVector& a) const;

  // destructive insert/delete/reorder operations

  ColumnVector& insert (const ColumnVector& a, octave_idx_type r);

  ColumnVector& fill (double val);
  ColumnVector&
  fill (double val, octave_idx_type r1, octave_idx_type r2);

  ColumnVector stack (const ColumnVector& a) const;

  RowVector transpose (void) const;

  friend ColumnVector real (const ComplexColumnVector& a);
  friend ColumnVector imag (const ComplexColumnVector& a);

  // resize is the destructive equivalent for this one

  ColumnVector
  extract (octave_idx_type r1, octave_idx_type r2) const;

  ColumnVector
  extract_n (octave_idx_type r1, octave_idx_type n) const;

  // matrix by column vector -> column vector operations

  friend ColumnVector operator * (const Matrix& a,
                                             const ColumnVector& b);

  // diagonal matrix by column vector -> column vector operations

  friend ColumnVector operator * (const DiagMatrix& a,
                                             const ColumnVector& b);

  // other operations

  double min (void) const;
  double max (void) const;

  ColumnVector abs (void) const;

  // i/o
    /*
  friend std::ostream& operator << (std::ostream& os,
                                               const ColumnVector& a);
  friend std::istream& operator >> (std::istream& is,
                                               ColumnVector& a);
    */
  void resize (octave_idx_type n, const double& rfv = 0);

  void clear (octave_idx_type n);
};