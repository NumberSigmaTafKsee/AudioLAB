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
#include <octave/ov-cell.h>
%}

class octave_cell : public octave_base_matrix<Cell>
{
public:

octave_cell (void);
octave_cell (const Cell& c);
octave_cell (const Array<std::string>& str);
octave_cell (const octave_cell& c);
~octave_cell (void) = default;

octave_base_value * clone (void) const;
octave_base_value * empty_clone (void) const;

void break_closure_cycles (const std::shared_ptr<octave::stack_frame>& frame);

octave_value subsref (const std::string& type,
                        const std::list<octave_value_list>& idx);
octave_value_list subsref (const std::string& type,
                            const std::list<octave_value_list>& idx,
                            int nargout);

octave_value subsref (const std::string& type,
                        const std::list<octave_value_list>& idx,
                        bool auto_add);

octave_value subsasgn (const std::string& type,
                        const std::list<octave_value_list>& idx,
                        const octave_value& rhs);

// FIXME: should we import the functions from the base class and
// overload them here, or should we use a different name so we don't
// have to do this?  Without the using declaration or a name change,
// the base class functions will be hidden.  That may be OK, but it
// can also cause some confusion.
using octave_base_value::assign;

void assign (const octave_value_list& idx, const Cell& rhs);

void assign (const octave_value_list& idx, const octave_value& rhs);

void delete_elements (const octave_value_list& idx);

std::size_t byte_size (void) const;

octave_value sort (octave_idx_type dim = 0, sortmode mode = ASCENDING) const;

octave_value sort (Array<octave_idx_type>& sidx, octave_idx_type dim = 0,
                    sortmode mode = ASCENDING) const;

sortmode issorted (sortmode mode = UNSORTED) const;

//Array<octave_idx_type> sort_rows_idx (sortmode mode = ASCENDING) const;

sortmode is_sorted_rows (sortmode mode = UNSORTED) const;

bool is_matrix_type (void) const;

bool isnumeric (void) const;

bool is_defined (void) const;

bool is_constant (void) const;

bool iscell (void) const;

builtin_type_t builtin_type (void) const;

bool iscellstr (void) const;

bool is_true (void) const;

//Cell cell_value (void) const;

octave_value_list list_value (void) const;

octave_value convert_to_str_internal (bool pad, bool, char type) const;

string_vector string_vector_value (bool pad = false) const;

//Array<std::string> cellstr_value (void) const;

//Array<std::string> cellstr_value (const char *fmt, ...) const;

bool print_as_scalar (void) const;

void print (std::ostream& os, bool pr_as_read_syntax = false);

void print_raw (std::ostream& os, bool pr_as_read_syntax = false) const;

bool print_name_tag (std::ostream& os, const std::string& name) const;

void short_disp (std::ostream& os) const;

bool save_ascii (std::ostream& os);

bool load_ascii (std::istream& is);

bool save_binary (std::ostream& os, bool save_as_floats);

bool load_binary (std::istream& is, bool swap,
                    octave::mach_info::float_format fmt);

bool save_hdf5 (octave_hdf5_id loc_id, const char *name, bool save_as_floats);

bool load_hdf5 (octave_hdf5_id loc_id, const char *name);

//octave_value map (unary_mapper_t umap) const;

mxArray * as_mxArray (bool interleaved) const;

// This function exists to support the MEX interface.
// You should not use it anywhere else.
const void * mex_get_data (void) const;
};
