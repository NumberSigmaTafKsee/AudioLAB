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

class RowVector : public MArray<double>
{
public:

  RowVector (void);
  explicit RowVector (octave_idx_type n);
  explicit RowVector (const dim_vector& dv);
  RowVector (octave_idx_type n, double val);
  RowVector (const RowVector& a);
  RowVector (const MArray<double>& a);
  RowVector (const Array<double>& a);

  RowVector& operator = (const RowVector& a);

  bool operator == (const RowVector& a) const;
  bool operator != (const RowVector& a) const;

  // destructive insert/delete/reorder operations
  RowVector& insert (const RowVector& a, octave_idx_type c);

  RowVector& fill (double val);
  RowVector& fill (double val, octave_idx_type c1, octave_idx_type c2);

  RowVector append (const RowVector& a) const;

  ColumnVector transpose (void) const;

  //friend RowVector real (const ComplexRowVector& a);
  //friend RowVector imag (const ComplexRowVector& a);

  // resize is the destructive equivalent for this one

  RowVector extract (octave_idx_type c1, octave_idx_type c2) const;

  RowVector extract_n (octave_idx_type c1, octave_idx_type n) const;

  // row vector by matrix -> row vector

  friend RowVector operator * (const RowVector& a, const Matrix& b);

  // other operations

  double min (void) const;
  double max (void) const;

  // i/o
    /*
  friend std::ostream& operator << (std::ostream& os,
                                               const RowVector& a);
  friend std::istream& operator >> (std::istream& is, RowVector& a);
    */
  void resize (octave_idx_type n, const double& rfv = 0);
  void clear (octave_idx_type n);
};