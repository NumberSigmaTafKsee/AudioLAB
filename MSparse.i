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

// Two dimensional sparse array with math ops.
template <typename T>
class MSparse : public Sparse<T>
{
public:

  MSparse (void);
  MSparse (octave_idx_type n, octave_idx_type m);
  MSparse (const dim_vector& dv, octave_idx_type nz = 0);
  MSparse (const MSparse<T>& a);
  MSparse (const MSparse<T>& a, const dim_vector& dv);
  MSparse (const Sparse<T>& a);

  template <typename U>
  MSparse (const Sparse<U>& a);
  MSparse (const Array<T>& a, const octave::idx_vector& r, const octave::idx_vector& c,
           octave_idx_type nr = -1, octave_idx_type nc = -1,
           bool sum_terms = true, octave_idx_type nzm = -1);
  explicit MSparse (octave_idx_type r, octave_idx_type c, T val);
  explicit MSparse (const PermMatrix& a) ;
  MSparse (octave_idx_type r, octave_idx_type c, octave_idx_type num_nz);
  ~MSparse (void) = default;

  MSparse<T>& operator = (const MSparse<T>& a);
  MSparse<T>& insert (const Sparse<T>& a, octave_idx_type r, octave_idx_type c);

  MSparse<T>& insert (const Sparse<T>& a, const Array<octave_idx_type>& indx);
  MSparse<T> transpose (void) const;
  MSparse<T> squeeze (void) const;
  MSparse<T> reshape (const dim_vector& new_dims) const;
  MSparse<T> permute (const Array<octave_idx_type>& vec, bool inv = false) const;
  MSparse<T> ipermute (const Array<octave_idx_type>& vec) const;
  MSparse<T> diag (octave_idx_type k = 0) const;
  
};