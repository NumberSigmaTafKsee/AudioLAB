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


template <typename T>
class MArray : public Array<T>
{
public:

  MArray (void);
  explicit MArray (const dim_vector& dv);
  explicit MArray (const dim_vector& dv, const T& val);
  MArray (const MArray<T>& a) : Array<T> (a);
  template <typename U>
  MArray (const Array<U>& a) : Array<T> (a);
  ~MArray (void) = default;
  MArray<T>& operator = (const MArray<T>& a);
  MArray<T> reshape (const dim_vector& new_dims) const;
  MArray<T> permute (const Array<octave_idx_type>& vec,
                     bool inv = false) const;
  MArray<T> ipermute (const Array<octave_idx_type>& vec) const;
  MArray squeeze (void) const;
  MArray<T> transpose (void) const;
  MArray<T> hermitian (T (*fcn) (const T&) = nullptr) const;
  void idx_add (const octave::idx_vector& idx, T val);
  void idx_add (const octave::idx_vector& idx, const MArray<T>& vals);
  
  void idx_min (const octave::idx_vector& idx, const MArray<T>& vals);

  void idx_max (const octave::idx_vector& idx, const MArray<T>& vals);

  void idx_add_nd (const octave::idx_vector& idx, const MArray<T>& vals,
              int dim = -1);

  void changesign (void);
};
