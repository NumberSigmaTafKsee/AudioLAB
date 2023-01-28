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

class octave_value
{
public:

enum unary_op
{
    op_not,            // not
    op_uplus,          // uplus
    op_uminus,         // uminus
    op_transpose,      // transpose
    op_hermitian,      // ctranspose
    op_incr,
    op_decr,
    num_unary_ops,
    unknown_unary_op
};

enum binary_op
{
    op_add,            // plus
    op_sub,            // minus
    op_mul,            // mtimes
    op_div,            // mrdivide
    op_pow,            // mpower
    op_ldiv,           // mldivide
    op_lt,             // lt
    op_le,             // le
    op_eq,             // eq
    op_ge,             // ge
    op_gt,             // gt
    op_ne,             // ne
    op_el_mul,         // times
    op_el_div,         // rdivide
    op_el_pow,         // power
    op_el_ldiv,        // ldivide
    op_el_and,         // and
    op_el_or,          // or
    op_struct_ref,
    num_binary_ops,
    unknown_binary_op
};

enum compound_binary_op
{
    // ** compound operations **
    op_trans_mul,
    op_mul_trans,
    op_herm_mul,
    op_mul_herm,
    op_trans_ldiv,
    op_herm_ldiv,
    op_el_not_and,
    op_el_not_or,
    op_el_and_not,
    op_el_or_not,
    num_compound_binary_ops,
    unknown_compound_binary_op
};

enum assign_op
{
    op_asn_eq,
    op_add_eq,
    op_sub_eq,
    op_mul_eq,
    op_div_eq,
    op_ldiv_eq,
    op_pow_eq,
    op_el_mul_eq,
    op_el_div_eq,
    op_el_ldiv_eq,
    op_el_pow_eq,
    op_el_and_eq,
    op_el_or_eq,
    num_assign_ops,
    unknown_assign_op
};

static  binary_op assign_op_to_binary_op (assign_op);

static  assign_op binary_op_to_assign_op (binary_op);

static  std::string unary_op_as_string (unary_op);
static  std::string unary_op_fcn_name (unary_op);

static  std::string binary_op_as_string (binary_op);
static  std::string binary_op_fcn_name (binary_op);

//static  std::string binary_op_fcn_name (compound_binary_op);

static  std::string assign_op_as_string (assign_op);

static  octave_value
empty_conv (const std::string& type,
            const octave_value& rhs = octave_value ());

enum magic_colon { magic_colon_t };

octave_value (void)
    : m_rep (nil_rep ())
{
    m_rep->count++;
}
    /*
octave_value (short int i);
octave_value (unsigned short int i);
octave_value (int i);
octave_value (unsigned int i);
octave_value (long int i);
octave_value (unsigned long int i);
    */
// FIXME: These are kluges.  They turn into doubles internally, which will
// break for very large values.  We just use them to store things like
// 64-bit ino_t, etc, and hope that those values are never actually larger
// than can be represented exactly in a double.
/*
#if defined (OCTAVE_HAVE_LONG_LONG_INT)
octave_value (long long int i);
#endif
#if defined (OCTAVE_HAVE_UNSIGNED_LONG_LONG_INT)
octave_value (unsigned long long int i);
#endif
*/
octave_value (octave::sys::time t);
octave_value (double d);
//octave_value (float d);
octave_value (const Array<octave_value>& a,
                            bool is_cs_list = false);
octave_value (const Cell& c, bool is_cs_list = false);
octave_value (const Matrix& m,
                            const MatrixType& t = MatrixType ());
octave_value (const FloatMatrix& m,
                            const MatrixType& t = MatrixType ());
octave_value (const NDArray& nda);
octave_value (const FloatNDArray& nda);
octave_value (const Array<double>& m);
octave_value (const Array<float>& m);
octave_value (const DiagMatrix& d);
octave_value (const DiagArray2<double>& d);
octave_value (const DiagArray2<float>& d);
octave_value (const DiagArray2<Complex>& d);
octave_value (const DiagArray2<FloatComplex>& d);
octave_value (const FloatDiagMatrix& d);
octave_value (const RowVector& v);
octave_value (const FloatRowVector& v);
octave_value (const ColumnVector& v);
octave_value (const FloatColumnVector& v);
octave_value (const Complex& C);
octave_value (const FloatComplex& C);
octave_value (const ComplexMatrix& m,
                            const MatrixType& t = MatrixType ());
octave_value (const FloatComplexMatrix& m,
                            const MatrixType& t = MatrixType ());
octave_value (const ComplexNDArray& cnda);
octave_value (const FloatComplexNDArray& cnda);
octave_value (const Array<Complex>& m);
octave_value (const Array<FloatComplex>& m);
octave_value (const ComplexDiagMatrix& d);
octave_value (const FloatComplexDiagMatrix& d);
octave_value (const ComplexRowVector& v);
octave_value (const FloatComplexRowVector& v);
octave_value (const ComplexColumnVector& v);
octave_value (const FloatComplexColumnVector& v);
octave_value (const PermMatrix& p);
octave_value (bool b);
octave_value (const boolMatrix& bm,
                            const MatrixType& t = MatrixType ());
octave_value (const boolNDArray& bnda);
octave_value (const Array<bool>& bnda);
octave_value (char c, char type = '\'');
//octave_value (const char *s, char type = '\'');
octave_value (const std::string& s, char type = '\'');
octave_value (const string_vector& s, char type = '\'');
octave_value (const charMatrix& chm,  char type = '\'');
octave_value (const charNDArray& chnda, char type = '\'');
octave_value (const Array<char>& chnda, char type = '\'');


octave_value (const SparseMatrix& m,
                            const MatrixType& t = MatrixType ());
octave_value (const Sparse<double>& m,
                            const MatrixType& t = MatrixType ());
octave_value (const SparseComplexMatrix& m,
                            const MatrixType& t = MatrixType ());
octave_value (const Sparse<Complex>& m,
                            const MatrixType& t = MatrixType ());
octave_value (const SparseBoolMatrix& bm,
                            const MatrixType& t = MatrixType ());
octave_value (const Sparse<bool>& m,
                            const MatrixType& t = MatrixType ());

octave_value (const octave_int8& i);
octave_value (const octave_int16& i);
octave_value (const octave_int32& i);
octave_value (const octave_int64& i);
octave_value (const octave_uint8& i);
octave_value (const octave_uint16& i);
octave_value (const octave_uint32& i);
octave_value (const octave_uint64& i);


octave_value (const int8NDArray& inda);
octave_value (const Array<octave_int8>& inda);
octave_value (const int16NDArray& inda);
octave_value (const Array<octave_int16>& inda);
octave_value (const int32NDArray& inda);
octave_value (const Array<octave_int32>& inda);
octave_value (const int64NDArray& inda);
octave_value (const Array<octave_int64>& inda);
octave_value (const uint8NDArray& inda);
octave_value (const Array<octave_uint8>& inda);
octave_value (const uint16NDArray& inda);
octave_value (const Array<octave_uint16>& inda);
octave_value (const uint32NDArray& inda);
octave_value (const Array<octave_uint32>& inda);
octave_value (const uint64NDArray& inda);
octave_value (const Array<octave_uint64>& inda);
octave_value (const Array<octave_idx_type>& inda,
                            bool zero_based = false,
                            bool cache_index = false);
octave_value (const Array<std::string>& cellstr);
octave_value (const octave::idx_vector& idx, bool lazy = true);


public:

#if defined (OCTAVE_PROVIDE_DEPRECATED_SYMBOLS)
OCTAVE_DEPRECATED (7, "use 'octave_value (range<double>&)' instead")
octave_value (double base, double limit, double inc);
OCTAVE_DEPRECATED (7, "use 'octave_value (range<double>&)' instead")
octave_value (const Range& r, bool force_range = false);
#endif

octave_value (const octave::range<double>& r,
                            bool force_range = false);

// For now, disable all but range<double>.

#if 0

octave_value (const octave::range<float>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_int8>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_int16>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_int32>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_int64>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_uint8>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_uint16>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_uint32>& r,
                            bool force_range = false);

octave_value (const octave::range<octave_uint64>& r,
                            bool force_range = false);

octave_value (const octave::range<char>& r, char type,
                            bool force_range = false);
#endif

octave_value (const octave_map& m);
octave_value (const octave_scalar_map& m);
octave_value (const std::map<std::string, octave_value>&);
octave_value (const octave_map& m, const std::string& id,
                            const std::list<std::string>& plist);
octave_value (const octave_scalar_map& m, const std::string& id,
                            const std::list<std::string>& plist);

// This one is explicit because it can cause some trouble to
// accidentally create a cs-list when one was not intended.
explicit  octave_value (const octave_value_list& m);

//octave_value (octave_value::magic_colon);

octave_value (octave_base_value *new_rep, bool borrow = false);

// Copy constructor.

octave_value (const octave_value& a);
//octave_value (octave_value&& a);

// This should only be called for derived types.

octave_base_value * clone (void) const;

octave_base_value * empty_clone (void) const;

// Delete the representation of this constant if the count drops to zero.

~octave_value (void);

void make_unique (void);
void make_unique (int obsolete_copies);
void break_closure_cycles (const std::shared_ptr<octave::stack_frame>&);

// Simple assignment.

octave_value& operator = (const octave_value& a);
octave_value& operator = (octave_value&& a);
octave_idx_type get_count (void) const;

octave_base_value::type_conv_info numeric_conversion_function (void) const;
octave_base_value::type_conv_info numeric_demotion_function (void) const;
void maybe_mutate (void);

octave_value squeeze (void) const;
octave_value full_value (void) const;
octave_value as_double (void) const;
octave_value as_single (void) const;

octave_value as_int8 (void) const;
octave_value as_int16 (void) const;
octave_value as_int32 (void) const;
octave_value as_int64 (void) const;

octave_value as_uint8 (void) const;
octave_value as_uint16 (void) const;
octave_value as_uint32 (void) const;
octave_value as_uint64 (void) const;

octave_base_value * try_narrowing_conversion (void);
Matrix size (void);
octave_idx_type xnumel (const octave_value_list& idx);
octave_value
single_subsref (const std::string& type, const octave_value_list& idx);

octave_value subsref (const std::string& type,
                        const std::list<octave_value_list>& idx);
octave_value subsref (const std::string& type,
                        const std::list<octave_value_list>& idx,
                        bool auto_add);
octave_value_list
subsref (const std::string& type, const std::list<octave_value_list>& idx,
        int nargout);

octave_value
next_subsref (const std::string& type,
                const std::list<octave_value_list>& idx, std::size_t skip = 1);

octave_value_list
next_subsref (int nargout, const std::string& type,
                const std::list<octave_value_list>& idx, std::size_t skip = 1);

octave_value
next_subsref (bool auto_add, const std::string& type,
                const std::list<octave_value_list>& idx, std::size_t skip = 1);

octave_value index_op (const octave_value_list& idx, bool resize_ok = false);

#if defined (OCTAVE_PROVIDE_DEPRECATED_SYMBOLS)
OCTAVE_DEPRECATED (7, "use 'octave_value::index_op' instead")
octave_value do_index_op (const octave_value_list& idx,
                            bool resize_ok = false);
#endif

octave_value
subsasgn (const std::string& type, const std::list<octave_value_list>& idx,
            const octave_value& rhs);

octave_value
undef_subsasgn (const std::string& type,
                const std::list<octave_value_list>& idx,
                const octave_value& rhs);

octave_value&
assign (assign_op op, const std::string& type,
        const std::list<octave_value_list>& idx, const octave_value& rhs);

octave_value& assign (assign_op, const octave_value& rhs);

octave::idx_vector index_vector (bool require_integers = false) const;
// Size.

dim_vector dims (void) const;

std::string get_dims_str (void) const;

octave_idx_type rows (void) const;

octave_idx_type columns (void) const;

octave_idx_type length (void) const;

int ndims (void) const;

bool all_zero_dims (void) const;

// Are the dimensions of this constant zero by zero?
bool is_zero_by_zero (void) const;
octave_idx_type numel (void) const;
std::size_t byte_size (void) const;
octave_idx_type nnz (void) const;
octave_idx_type nzmax (void) const;
octave_idx_type nfields (void) const;

octave_value reshape (const dim_vector& dv) const;
octave_value permute (const Array<int>& vec, bool inv = false) const;
octave_value ipermute (const Array<int>& vec) const;
octave_value resize (const dim_vector& dv, bool fill = false) const;
MatrixType matrix_type (void) const;
MatrixType matrix_type (const MatrixType& typ) const;

bool is_defined (void) const;
bool is_undefined (void) const;
bool is_legacy_object (void) const;
bool isempty (void) const;
bool iscell (void) const;
bool iscellstr (void) const;
bool is_real_scalar (void) const;
bool is_real_matrix (void) const;
bool is_complex_scalar (void) const;
bool is_complex_matrix (void) const;
bool is_bool_scalar (void) const;
bool is_bool_matrix (void) const;
bool is_char_matrix (void) const;
bool is_diag_matrix (void) const;
bool is_perm_matrix (void) const;
bool is_string (void) const;
bool is_sq_string (void) const;
bool is_dq_string (void) const;
bool is_range (void) const;
bool isstruct (void) const;
bool is_classdef_meta (void) const;
bool is_classdef_object (void) const;
bool is_classdef_superclass_ref (void) const;
bool is_package (void) const;
bool isobject (void) const;
bool isjava (void) const;
bool is_cs_list (void) const;
bool is_magic_colon (void) const;
bool is_magic_int (void) const;
bool isnull (void) const;
octave_value all (int dim = 0) const;
octave_value any (int dim = 0) const;
builtin_type_t builtin_type (void) const;
bool is_double_type (void) const;
bool is_single_type (void) const;
bool isfloat (void) const;
bool is_int8_type (void) const;
bool is_int16_type (void) const;
bool is_int32_type (void) const;
bool is_int64_type (void) const;
bool is_uint8_type (void) const;
bool is_uint16_type (void) const;
bool is_uint32_type (void) const;
bool is_uint64_type (void) const;
bool isinteger (void) const;
bool islogical (void) const;
bool isreal (void) const;
bool iscomplex (void) const;
bool is_scalar_type (void) const;
bool is_matrix_type (void) const;
bool isnumeric (void) const;
bool issparse (void) const;
bool is_true (void) const;
bool is_equal (const octave_value&) const;
bool is_constant (void) const;
bool is_function_handle (void) const;
bool is_anonymous_function (void) const;
bool is_inline_function (void) const;
bool is_function (void) const;
bool is_user_script (void) const;
bool is_user_function (void) const;
bool is_user_code (void) const;
bool is_builtin_function (void) const;
bool is_dld_function (void) const;
bool is_mex_function (void) const;
void erase_subfunctions (void);
octave_value eval (void);
short int
short_value (bool req_int = false, bool frc_str_conv = false) const;
unsigned short int
ushort_value (bool req_int = false, bool frc_str_conv = false) const;
int int_value (bool req_int = false, bool frc_str_conv = false) const;
unsigned int
uint_value (bool req_int = false, bool frc_str_conv = false) const;
int nint_value (bool frc_str_conv = false) const;
long int
long_value (bool req_int = false, bool frc_str_conv = false) const;
unsigned long int
ulong_value (bool req_int = false, bool frc_str_conv = false) const;
int64_t
int64_value (bool req_int = false, bool frc_str_conv = false) const;
uint64_t
uint64_value (bool req_int = false, bool frc_str_conv = false) const;
octave_idx_type
idx_type_value (bool req_int = false, bool frc_str_conv = false) const;

double double_value (bool frc_str_conv = false) const;
float float_value (bool frc_str_conv = false) const;
double scalar_value (bool frc_str_conv = false) const;
float float_scalar_value (bool frc_str_conv = false) const;
Matrix matrix_value (bool frc_str_conv = false) const;
FloatMatrix float_matrix_value (bool frc_str_conv = false) const;
NDArray array_value (bool frc_str_conv = false) const;
FloatNDArray float_array_value (bool frc_str_conv = false) const;
Complex complex_value (bool frc_str_conv = false) const;
FloatComplex float_complex_value (bool frc_str_conv = false) const;
ComplexMatrix complex_matrix_value (bool frc_str_conv = false) const;
FloatComplexMatrix
float_complex_matrix_value (bool frc_str_conv = false) const;
ComplexNDArray complex_array_value (bool frc_str_conv = false) const;
FloatComplexNDArray
float_complex_array_value (bool frc_str_conv = false) const;
bool bool_value (bool warn = false) const;
boolMatrix bool_matrix_value (bool warn = false) const;
boolNDArray bool_array_value (bool warn = false) const;
charMatrix char_matrix_value (bool frc_str_conv = false) const;
charNDArray char_array_value (bool frc_str_conv = false) const;
SparseMatrix sparse_matrix_value (bool frc_str_conv = false) const;
SparseComplexMatrix
sparse_complex_matrix_value (bool frc_str_conv = false) const;
SparseBoolMatrix sparse_bool_matrix_value (bool warn = false) const;
DiagMatrix diag_matrix_value (bool force = false) const;
FloatDiagMatrix float_diag_matrix_value (bool force = false) const;
ComplexDiagMatrix complex_diag_matrix_value (bool force = false) const;
FloatComplexDiagMatrix
float_complex_diag_matrix_value (bool force = false) const;
PermMatrix perm_matrix_value (void) const;
octave_int8 int8_scalar_value (void) const;
octave_int16 int16_scalar_value (void) const;
octave_int32 int32_scalar_value (void) const;
octave_int64 int64_scalar_value (void) const;
octave_uint8 uint8_scalar_value (void) const;
octave_uint16 uint16_scalar_value (void) const;
octave_uint32 uint32_scalar_value (void) const;
octave_uint64 uint64_scalar_value (void) const;
int8NDArray int8_array_value (void) const;
int16NDArray int16_array_value (void) const;
int32NDArray int32_array_value (void) const;
int64NDArray int64_array_value (void) const;
uint8NDArray uint8_array_value (void) const;
uint16NDArray uint16_array_value (void) const;
uint32NDArray uint32_array_value (void) const;
uint64NDArray uint64_array_value (void) const;
std::string string_value (bool force = false) const;
string_vector string_vector_value (bool pad = false) const;
Cell cell_value (void) const;
Array<std::string> cellstr_value (void) const;
octave::range<double> range_value (void) const;
octave_map map_value (void) const;
octave_scalar_map scalar_map_value (void) const;
string_vector map_keys (void) const;
std::size_t nparents (void) const;
std::list<std::string> parent_class_name_list (void) const;
string_vector parent_class_names (void) const;
octave_base_value *
find_parent_class (const std::string& parent_class_name);
bool is_instance_of (const std::string& cls_name) const;
octave_classdef *
classdef_object_value (bool silent = false) const;

octave_function *
function_value (bool silent = false) const;

octave_user_function *
user_function_value (bool silent = false) const;

octave_user_script *
user_script_value (bool silent = false) const;

octave_user_code * user_code_value (bool silent = false) const;

octave_fcn_handle *
fcn_handle_value (bool silent = false) const;

octave_value_list list_value (void) const;

ColumnVector
column_vector_value (bool frc_str_conv = false,
                    bool frc_vec_conv = false) const;

ComplexColumnVector
complex_column_vector_value (bool frc_str_conv = false,
                            bool frc_vec_conv = false) const;

RowVector
row_vector_value (bool frc_str_conv = false,
                    bool frc_vec_conv = false) const;

ComplexRowVector
complex_row_vector_value (bool frc_str_conv = false,
                            bool frc_vec_conv = false) const;

FloatColumnVector
float_column_vector_value (bool frc_str_conv = false,
                            bool frc_vec_conv = false) const;

FloatComplexColumnVector
float_complex_column_vector_value (bool frc_str_conv = false,
                                    bool frc_vec_conv = false) const;

FloatRowVector
float_row_vector_value (bool frc_str_conv = false,
                        bool frc_vec_conv = false) const;

FloatComplexRowVector
float_complex_row_vector_value (bool frc_str_conv = false,
                                bool frc_vec_conv = false) const;

Array<int>
int_vector_value (bool req_int = false,
                    bool frc_str_conv = false,
                    bool frc_vec_conv = false) const;

Array<octave_idx_type>
octave_idx_type_vector_value (bool req_int = false,
                                bool frc_str_conv = false,
                                bool frc_vec_conv = false) const;

Array<double>
vector_value (bool frc_str_conv = false,
                bool frc_vec_conv = false) const;

Array<Complex>
complex_vector_value (bool frc_str_conv = false,
                        bool frc_vec_conv = false) const;

Array<float>
float_vector_value (bool frc_str_conv = false,
                    bool frc_vec_conv = false) const;

Array<FloatComplex>
float_complex_vector_value (bool frc_str_conv = false,
                            bool frc_vec_conv = false) const;

// Extract values of specific types without any implicit type conversions.
// Throw an error if an object is the wrong type for the requested value
// extraction.
//
// These functions are intended to provide a simple way to extract values of
// specific types and display error messages that are more meaningful than
// the generic "error: wrong type argument 'cell'" message.
/*
short int xshort_value (const char *fmt, ...) const;

unsigned short int xushort_value (const char *fmt, ...) const;

int xint_value (const char *fmt, ...) const;

unsigned int xuint_value (const char *fmt, ...) const;

int xnint_value (const char *fmt, ...) const;

long int xlong_value (const char *fmt, ...) const;

unsigned long int xulong_value (const char *fmt, ...) const;

int64_t xint64_value (const char *fmt, ...) const;

uint64_t xuint64_value (const char *fmt, ...) const;

octave_idx_type xidx_type_value (const char *fmt, ...) const;

double xdouble_value (const char *fmt, ...) const;

float xfloat_value (const char *fmt, ...) const;

double xscalar_value (const char *fmt, ...) const;

float xfloat_scalar_value (const char *fmt, ...) const;

Matrix xmatrix_value (const char *fmt, ...) const;

FloatMatrix xfloat_matrix_value (const char *fmt, ...) const;

NDArray xarray_value (const char *fmt, ...) const;

FloatNDArray xfloat_array_value (const char *fmt, ...) const;

Complex xcomplex_value (const char *fmt, ...) const;

FloatComplex xfloat_complex_value (const char *fmt, ...) const;

ComplexMatrix
xcomplex_matrix_value (const char *fmt, ...) const;

FloatComplexMatrix
xfloat_complex_matrix_value (const char *fmt, ...) const;

ComplexNDArray
xcomplex_array_value (const char *fmt, ...) const;

FloatComplexNDArray
xfloat_complex_array_value (const char *fmt, ...) const;

bool xbool_value (const char *fmt, ...) const;

boolMatrix xbool_matrix_value (const char *fmt, ...) const;

boolNDArray xbool_array_value (const char *fmt, ...) const;

charMatrix xchar_matrix_value (const char *fmt, ...) const;

charNDArray xchar_array_value (const char *fmt, ...) const;

SparseMatrix xsparse_matrix_value (const char *fmt, ...) const;

SparseComplexMatrix
xsparse_complex_matrix_value (const char *fmt, ...) const;

SparseBoolMatrix
xsparse_bool_matrix_value (const char *fmt, ...) const;

DiagMatrix xdiag_matrix_value (const char *fmt, ...) const;

FloatDiagMatrix
xfloat_diag_matrix_value (const char *fmt, ...) const;

ComplexDiagMatrix
xcomplex_diag_matrix_value (const char *fmt, ...) const;

FloatComplexDiagMatrix
xfloat_complex_diag_matrix_value (const char *fmt, ...) const;

PermMatrix xperm_matrix_value (const char *fmt, ...) const;

octave_int8 xint8_scalar_value (const char *fmt, ...) const;

octave_int16 xint16_scalar_value (const char *fmt, ...) const;

octave_int32 xint32_scalar_value (const char *fmt, ...) const;

octave_int64 xint64_scalar_value (const char *fmt, ...) const;

octave_uint8 xuint8_scalar_value (const char *fmt, ...) const;

octave_uint16 xuint16_scalar_value (const char *fmt, ...) const;

octave_uint32 xuint32_scalar_value (const char *fmt, ...) const;

octave_uint64 xuint64_scalar_value (const char *fmt, ...) const;

int8NDArray xint8_array_value (const char *fmt, ...) const;

int16NDArray xint16_array_value (const char *fmt, ...) const;

int32NDArray xint32_array_value (const char *fmt, ...) const;

int64NDArray xint64_array_value (const char *fmt, ...) const;

uint8NDArray xuint8_array_value (const char *fmt, ...) const;

uint16NDArray xuint16_array_value (const char *fmt, ...) const;

uint32NDArray xuint32_array_value (const char *fmt, ...) const;

uint64NDArray xuint64_array_value (const char *fmt, ...) const;

std::string xstring_value (const char *fmt, ...) const;

string_vector xstring_vector_value (const char *fmt, ...) const;

Cell xcell_value (const char *fmt, ...) const;

Array<std::string> xcellstr_value (const char *fmt, ...) const;

octave::range<double>
xrange_value (const char *fmt, ...) const;

// For now, disable all but range<double>.

#if 0

octave::range<float>
xfloat_range_value (const char *fmt, ...) const;

octave::range<octave_int8>
xint8_range_value (const char *fmt, ...) const;

octave::range<octave_int16>
xint16_range_value (const char *fmt, ...) const;

octave::range<octave_int32>
xint32_range_value (const char *fmt, ...) const;

octave::range<octave_int64>
xint64_range_value (const char *fmt, ...) const;

octave::range<octave_uint8>
xuint8_range_value (const char *fmt, ...) const;

octave::range<octave_uint16>
xuint16_range_value (const char *fmt, ...) const;

octave::range<octave_uint32>
xuint32_range_value (const char *fmt, ...) const;

octave::range<octave_uint64>
xuint64_range_value (const char *fmt, ...) const;

#endif

octave_map xmap_value (const char *fmt, ...) const;

octave_scalar_map
xscalar_map_value (const char *fmt, ...) const;

ColumnVector xcolumn_vector_value (const char *fmt, ...) const;

ComplexColumnVector
xcomplex_column_vector_value (const char *fmt, ...) const;

RowVector xrow_vector_value (const char *fmt, ...) const;

ComplexRowVector
xcomplex_row_vector_value (const char *fmt, ...) const;

FloatColumnVector
xfloat_column_vector_value (const char *fmt, ...) const;

FloatComplexColumnVector
xfloat_complex_column_vector_value (const char *fmt, ...) const;

FloatRowVector
xfloat_row_vector_value (const char *fmt, ...) const;

FloatComplexRowVector
xfloat_complex_row_vector_value (const char *fmt, ...) const;

Array<int> xint_vector_value (const char *fmt, ...) const;

Array<octave_idx_type>
xoctave_idx_type_vector_value (const char *fmt, ...) const;

Array<double> xvector_value (const char *fmt, ...) const;

Array<Complex>
xcomplex_vector_value (const char *fmt, ...) const;

Array<float> xfloat_vector_value (const char *fmt, ...) const;

Array<FloatComplex>
xfloat_complex_vector_value (const char *fmt, ...) const;

octave_function * xfunction_value (const char *fmt, ...) const;

octave_user_function *
xuser_function_value (const char *fmt, ...) const;

octave_user_script *
xuser_script_value (const char *fmt, ...) const;

octave_user_code *
xuser_code_value (const char *fmt, ...) const;

octave_fcn_handle *
xfcn_handle_value (const char *fmt, ...) const;

octave_value_list xlist_value (const char *fmt, ...) const;
*/

// Possibly economize a lazy-indexed value.

void maybe_economize (void);

octave_value storable_value (void) const;

void make_storable_value (void);

octave_value convert_to_str (bool pad = false, bool force = false,
                            char type = '\'') const;
octave_value
convert_to_str_internal (bool pad, bool force, char type) const;
void convert_to_row_or_column_vector (void);
bool print_as_scalar (void) const;
void print (std::ostream& os, bool pr_as_read_syntax = false);
void print_raw (std::ostream& os, bool pr_as_read_syntax = false) const;
bool print_name_tag (std::ostream& os, const std::string& name) const;
void print_with_name (std::ostream& os, const std::string& name) const;
void short_disp (std::ostream& os) const { m_rep->short_disp (os); }

//float_display_format get_edit_display_format (void) const;

std::string edit_display (const float_display_format& fmt,
                            octave_idx_type i, octave_idx_type j) const;

int type_id (void) const;

std::string type_name (void) const;

std::string class_name (void) const;

// Unary operations that are member functions.  There are also some
// non-member functions for unary and binary operations declared
// below, outside of the octave_value class declaration.

octave_value& non_const_unary_op (unary_op op);

octave_value&
non_const_unary_op (unary_op op, const std::string& type,
                    const std::list<octave_value_list>& idx);

#if defined (OCTAVE_PROVIDE_DEPRECATED_SYMBOLS)
OCTAVE_DEPRECATED (7, "use 'octave_value::non_const_unary_op' instead")
octave_value& do_non_const_unary_op (unary_op op, const std::string& type,
                                    const std::list<octave_value_list>& idx);
#endif

const octave_base_value& get_rep (void) const;

bool is_copy_of (const octave_value& val) const;

void
print_info (std::ostream& os, const std::string& prefix = "") const;

bool save_ascii (std::ostream& os);

bool load_ascii (std::istream& is);

bool save_binary (std::ostream& os, bool save_as_floats);
bool load_binary (std::istream& is, bool swap,
                                octave::mach_info::float_format fmt);

bool save_hdf5 (octave_hdf5_id loc_id, const char *name,
                bool save_as_floats);
bool load_hdf5 (octave_hdf5_id loc_id, const char *name);

int
write (octave::stream& os, int block_size,
        oct_data_conv::data_type output_type, int skip,
        octave::mach_info::float_format flt_fmt) const;

octave_base_value * internal_rep (void) const;

// These functions exist to support the MEX interface.
// You should not use them anywhere else.

const void *
mex_get_data (mxClassID class_id = mxUNKNOWN_CLASS,
                mxComplexity complexity = mxREAL) const;

const octave_idx_type * mex_get_ir (void) const;
const octave_idx_type *
mex_get_jc (void) const;

mxArray * as_mxArray (bool interleaved = false) const;
octave_value diag (octave_idx_type k = 0) const;
octave_value diag (octave_idx_type m, octave_idx_type n) const;
octave_value sort (octave_idx_type dim = 0, sortmode mode = ASCENDING) const;
octave_value sort (Array<octave_idx_type>& sidx, octave_idx_type dim = 0,
                    sortmode mode = ASCENDING) const;
sortmode issorted (sortmode mode = UNSORTED) const;

Array<octave_idx_type> sort_rows_idx (sortmode mode = ASCENDING) const;
sortmode is_sorted_rows (sortmode mode = UNSORTED) const;
void lock (void);

void unlock (void);

bool islocked (void) const;

void call_object_destructor (void);

octave_value dump (void) const;

octave_value abs();
octave_value acos();
octave_value acosh();
octave_value angle();
octave_value arg();
octave_value asin();
octave_value asinh();
octave_value atan();
octave_value atanh();
octave_value cbrt();
octave_value ceil();
octave_value conj();
octave_value cos();
octave_value cosh();
octave_value erf();
octave_value erfinv();
octave_value erfcinv();
octave_value erfc();
octave_value erfcx();
octave_value erfi();
octave_value dawson();
octave_value exp();
octave_value expm1();
octave_value isfinite();
octave_value fix();
octave_value floor();
octave_value gamma();
octave_value imag();
octave_value isinf();
octave_value isna();
octave_value isnan();
octave_value lgamma();
octave_value log();
octave_value log2();
octave_value log10();
octave_value log1p();
octave_value real();
octave_value round();
octave_value roundb();
octave_value signum();
octave_value sin();
octave_value sinh();
octave_value sqrt();
octave_value tan();
octave_value tanh();

// These functions are prefixed with X to avoid potential macro conflicts.

octave_value xisalnum();
octave_value xisalpha();
octave_value xisascii();
octave_value xiscntrl();
octave_value xisdigit();
octave_value xisgraph();
octave_value xislower();
octave_value xisprint();
octave_value xispunct();
octave_value xisspace();
octave_value xisupper();
octave_value xisxdigit();
octave_value xsignbit();
octave_value xtolower();
octave_value xtoupper();



octave_value map (octave_base_value::unary_mapper_t umap) const;
octave_value
fast_elem_extract (octave_idx_type n) const;
bool
fast_elem_insert (octave_idx_type n, const octave_value& x);

};




