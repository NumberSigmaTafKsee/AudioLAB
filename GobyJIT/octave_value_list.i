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
class octave_value_list
{
public:

octave_value_list (void) = default;

explicit octave_value_list (octave_idx_type n);
octave_value_list (octave_idx_type n, const octave_value& val);
octave_value_list (const octave_value& tc);

octave_value_list (const Array<octave_value>& a);
octave_value_list (const Cell& c);
octave_value_list (const octave_value_list& obj);
//octave_value_list (octave_value_list&& obj);
// Concatenation constructors.
octave_value_list (const std::list<octave_value>&);
octave_value_list (const std::list<octave_value_list>&);
~octave_value_list (void);

octave_value_list& operator = (const octave_value_list& obj);
//octave_value_list& operator = (octave_value_list&& obj);

//Array<octave_value> array_value (void) const;
//Cell cell_value (void) const;

// Assignment will resize on range errors.
//octave_value& operator ();  
//const octave_value& operator () (octave_idx_type n) const;

%extend {
    octave_value __getitem__(size_t i) { return (*$self)(i); }
    void __setitem__(size_t i, octave_value& v) { (*$self)(i) = v; }
    void __setitem__(size_t i, double v) { (*$self)(i) = v; }
    void __setitem__(size_t i, const std::string& v) { (*$self)(i) = v; }
    void __setitem__(size_t i, bool v) { (*$self)(i) = v; }
}

octave_idx_type length (void) const;

bool empty (void) const;

void resize (octave_idx_type n, const octave_value& rfv = octave_value ());
octave_value_list& prepend (const octave_value& val);
octave_value_list& append (const octave_value& val);
octave_value_list& append (const octave_value_list& lst);
octave_value_list& reverse (void);
octave_value_list
slice (octave_idx_type offset, octave_idx_type len, bool tags = false) const;
octave_value_list
splice (octave_idx_type offset, octave_idx_type len,
        const octave_value_list& lst = octave_value_list ()) const;

bool all_strings_p (void) const;
bool all_scalars (void) const;
bool any_cell (void) const;
bool has_magic_colon (void) const;
string_vector make_argv (const std::string& = "") const;
void stash_name_tags (const string_vector& nm);
string_vector name_tags (void) const;
void make_storable_values (void);
octave_value& xelem (octave_idx_type i);
void clear (void);
};
