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

class dim_vector
{
public:


  template <typename... Ints>
  dim_vector (const octave_idx_type r, const octave_idx_type c,
              Ints... lengths);
  
  octave_idx_type& xelem (int i);
  //octave_idx_type xelem (int i) const;
  octave_idx_type& elem (int i);
  //octave_idx_type elem (int i) const;
  void chop_trailing_singletons (void);
  void chop_all_singletons (void);


  static octave_idx_type dim_max (void);

  explicit dim_vector (void);
  //dim_vector (const dim_vector& dv);
  dim_vector (dim_vector&& dv);
  static dim_vector alloc (int n);
  dim_vector& operator = (const dim_vector& dv);

  dim_vector& operator = (dim_vector&& dv);

  ~dim_vector (void);

  octave_idx_type ndims (void) const;
  int length (void) const;
  octave_idx_type& operator () (int i);
  //octave_idx_type operator () (int i) const;
  void resize (int n, int fill_value = 0);
  std::string str (char sep = 'x') const;
  bool all_zero (void) const;
  bool empty_2d (void) const;
  bool zero_by_zero (void) const;
  bool any_zero (void) const;
  int num_ones (void) const;
  bool all_ones (void) const;
  octave_idx_type numel (int n = 0) const;
  octave_idx_type safe_numel (void) const;
  bool any_neg (void) const;
  dim_vector squeeze (void) const;
  bool concat (const dim_vector& dvb, int dim);
  bool hvcat (const dim_vector& dvb, int dim);
  dim_vector redim (int n) const;
  dim_vector as_column (void) const;
  dim_vector as_row (void) const;
  bool isvector (void) const;
  bool is_nd_vector (void) const;;
  dim_vector make_nd_vector (octave_idx_type n) const;
  int first_non_singleton (int def = 0) const;
  octave_idx_type compute_index (const octave_idx_type *idx) const;
  octave_idx_type compute_index (const octave_idx_type *idx, int nidx) const;
  int increment_index (octave_idx_type *idx, int start = 0) const;
  dim_vector cumulative (void) const;
  octave_idx_type cum_compute_index (const octave_idx_type *idx) const;
  Array<octave_idx_type> as_array (void) const;

};