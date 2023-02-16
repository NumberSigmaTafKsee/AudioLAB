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

class SparseComplexMatrix : public MSparse<Complex>
{
public:

  // Corresponding dense matrix type for this sparse matrix type.
  typedef ComplexMatrix dense_matrix_type;

  typedef void (*solve_singularity_handler) (double rcond);

  SparseComplexMatrix (void) : MSparse<Complex> ();
  SparseComplexMatrix (octave_idx_type r,
                       octave_idx_type c);
  SparseComplexMatrix (const dim_vector& dv, octave_idx_type nz = 0);

  explicit SparseComplexMatrix (octave_idx_type r, octave_idx_type c,
                                Complex val);
  SparseComplexMatrix (octave_idx_type r, octave_idx_type c, double val);
  SparseComplexMatrix (const SparseComplexMatrix& a);
  SparseComplexMatrix (const SparseComplexMatrix& a, const dim_vector& dv);
  SparseComplexMatrix (const MSparse<Complex>& a) : MSparse<Complex> (a);
  SparseComplexMatrix (const Sparse<Complex>& a) : MSparse<Complex> (a);
  explicit SparseComplexMatrix (const ComplexMatrix& a);
  explicit SparseComplexMatrix (const ComplexNDArray& a);
  SparseComplexMatrix (const Array<Complex>& a, const octave::idx_vector& r,
                       const octave::idx_vector& c, octave_idx_type nr = -1,
                       octave_idx_type nc = -1, bool sum_terms = true,
                       octave_idx_type nzm = -1);

  explicit SparseComplexMatrix (const SparseMatrix& a);

  explicit SparseComplexMatrix (const SparseBoolMatrix& a);

  explicit SparseComplexMatrix (const ComplexDiagMatrix& a);

  SparseComplexMatrix (octave_idx_type r, octave_idx_type c,
                       octave_idx_type num_nz);
  SparseComplexMatrix& operator = (const SparseComplexMatrix& a);

  bool operator == (const SparseComplexMatrix& a) const;
  bool operator != (const SparseComplexMatrix& a) const;

  bool ishermitian (void) const;

  SparseComplexMatrix max (int dim = -1) const;
  SparseComplexMatrix
  max (Array<octave_idx_type>& index, int dim = -1) const;
  SparseComplexMatrix min (int dim = -1) const;
  SparseComplexMatrix
  min (Array<octave_idx_type>& index, int dim = -1) const;

  SparseComplexMatrix&
  insert (const SparseComplexMatrix& a, octave_idx_type r, octave_idx_type c);
  SparseComplexMatrix&
  insert (const SparseMatrix& a, octave_idx_type r, octave_idx_type c);
  SparseComplexMatrix&
  insert (const SparseComplexMatrix& a, const Array<octave_idx_type>& indx);
  SparseComplexMatrix&
  insert (const SparseMatrix& a, const Array<octave_idx_type>& indx);

  SparseComplexMatrix
  concat (const SparseComplexMatrix& rb, const Array<octave_idx_type>& ra_idx);
  SparseComplexMatrix
  concat (const SparseMatrix& rb, const Array<octave_idx_type>& ra_idx);

  ComplexMatrix matrix_value (void) const;

  SparseComplexMatrix hermitian (void) const;  // complex conjugate transpose
  SparseComplexMatrix transpose (void) const
  { return MSparse<Complex>::transpose (); }

  friend SparseComplexMatrix conj (const SparseComplexMatrix& a);

  // extract row or column i.

  ComplexRowVector row (octave_idx_type i) const;

  ComplexColumnVector column (octave_idx_type i) const;


public:
  SparseComplexMatrix inverse (void) const;
  SparseComplexMatrix inverse (MatrixType& mattype) const;
  SparseComplexMatrix
  inverse (MatrixType& mattype, octave_idx_type& info) const;
  SparseComplexMatrix
  inverse (MatrixType& mattype, octave_idx_type& info, double& rcond,
           bool force = false, bool calc_cond = true) const;

  ComplexDET determinant (void) const;
  ComplexDET determinant (octave_idx_type& info) const;
  ComplexDET
  determinant (octave_idx_type& info, double& rcond,
               bool calc_cond = true) const;

public:
  // Generic interface to solver with no probing of type
  ComplexMatrix solve (MatrixType& mattype, const Matrix& b) const;
  ComplexMatrix
  solve (MatrixType& mattype, const Matrix& b, octave_idx_type& info) const;
  ComplexMatrix
  solve (MatrixType& mattype, const Matrix& b, octave_idx_type& info,
         double& rcond) const;
  ComplexMatrix
  solve (MatrixType& mattype, const Matrix& b, octave_idx_type& info,
         double& rcond, solve_singularity_handler sing_handler,
         bool singular_fallback = true) const;

  ComplexMatrix
  solve (MatrixType& mattype, const ComplexMatrix& b) const;
  ComplexMatrix
  solve (MatrixType& mattype, const ComplexMatrix& b,
         octave_idx_type& info) const;
  ComplexMatrix
  solve (MatrixType& mattype, const ComplexMatrix& b,
         octave_idx_type& info, double& rcond) const;
  ComplexMatrix
  solve (MatrixType& mattype, const ComplexMatrix& b,
         octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler,
         bool singular_fallback = true) const;

  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseMatrix& b) const;
  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseMatrix& b,
         octave_idx_type& info) const;
  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseMatrix& b, octave_idx_type& info,
         double& rcond) const;
  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseMatrix& b, octave_idx_type& info,
         double& rcond, solve_singularity_handler sing_handler,
         bool singular_fallback = true) const;

  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseComplexMatrix& b) const;
  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseComplexMatrix& b,
         octave_idx_type& info) const;
  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseComplexMatrix& b,
         octave_idx_type& info, double& rcond) const;
  SparseComplexMatrix
  solve (MatrixType& mattype, const SparseComplexMatrix& b,
         octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler,
         bool singular_fallback = true) const;

  ComplexColumnVector
  solve (MatrixType& mattype, const ColumnVector& b) const;
  ComplexColumnVector
  solve (MatrixType& mattype, const ColumnVector& b,
         octave_idx_type& info) const;
  ComplexColumnVector
  solve (MatrixType& mattype, const ColumnVector& b,
         octave_idx_type& info, double& rcond) const;
  ComplexColumnVector
  solve (MatrixType& mattype, const ColumnVector& b,
         octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  ComplexColumnVector
  solve (MatrixType& mattype, const ComplexColumnVector& b) const;
  ComplexColumnVector
  solve (MatrixType& mattype, const ComplexColumnVector& b,
         octave_idx_type& info) const;
  ComplexColumnVector
  solve (MatrixType& mattype, const ComplexColumnVector& b,
         octave_idx_type& info, double& rcond) const;
  ComplexColumnVector
  solve (MatrixType& mattype, const ComplexColumnVector& b,
         octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  // Generic interface to solver with probing of type
  ComplexMatrix solve (const Matrix& b) const;
  ComplexMatrix
  solve (const Matrix& b, octave_idx_type& info) const;
  ComplexMatrix
  solve (const Matrix& b, octave_idx_type& info, double& rcond) const;
  ComplexMatrix
  solve (const Matrix& b, octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  ComplexMatrix solve (const ComplexMatrix& b) const;
  ComplexMatrix
  solve (const ComplexMatrix& b, octave_idx_type& info) const;
  ComplexMatrix
  solve (const ComplexMatrix& b, octave_idx_type& info, double& rcond) const;
  ComplexMatrix
  solve (const ComplexMatrix& b, octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  SparseComplexMatrix solve (const SparseMatrix& b) const;
  SparseComplexMatrix
  solve (const SparseMatrix& b, octave_idx_type& info) const;
  SparseComplexMatrix
  solve (const SparseMatrix& b, octave_idx_type& info, double& rcond) const;
  SparseComplexMatrix
  solve (const SparseMatrix& b, octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  SparseComplexMatrix solve (const SparseComplexMatrix& b) const;
  SparseComplexMatrix
  solve (const SparseComplexMatrix& b, octave_idx_type& info) const;
  SparseComplexMatrix
  solve (const SparseComplexMatrix& b, octave_idx_type& info,
         double& rcond) const;
  SparseComplexMatrix
  solve (const SparseComplexMatrix& b, octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  ComplexColumnVector solve (const ColumnVector& b) const;
  ComplexColumnVector
  solve (const ColumnVector& b, octave_idx_type& info) const;
  ComplexColumnVector
  solve (const ColumnVector& b, octave_idx_type& info, double& rcond) const;
  ComplexColumnVector
  solve (const ColumnVector& b, octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  ComplexColumnVector solve (const ComplexColumnVector& b) const;
  ComplexColumnVector
  solve (const ComplexColumnVector& b, octave_idx_type& info) const;
  ComplexColumnVector
  solve (const ComplexColumnVector& b, octave_idx_type& info,
         double& rcond) const;
  ComplexColumnVector
  solve (const ComplexColumnVector& b, octave_idx_type& info, double& rcond,
         solve_singularity_handler sing_handler) const;

  SparseComplexMatrix squeeze (void) const;

  SparseComplexMatrix reshape (const dim_vector& new_dims) const;

  SparseComplexMatrix
  permute (const Array<octave_idx_type>& vec, bool inv = false) const;

  SparseComplexMatrix
  ipermute (const Array<octave_idx_type>& vec) const;

  bool any_element_is_nan (void) const;
  bool any_element_is_inf_or_nan (void) const;
  bool all_elements_are_real (void) const;
  bool all_integers (double& max_val, double& min_val) const;
  bool too_large_for_float (void) const;

  SparseBoolMatrix operator ! (void) const;

  SparseBoolMatrix all (int dim = -1) const;
  SparseBoolMatrix any (int dim = -1) const;

  SparseComplexMatrix cumprod (int dim = -1) const;
  SparseComplexMatrix cumsum (int dim = -1) const;
  SparseComplexMatrix prod (int dim = -1) const;
  SparseComplexMatrix sum (int dim = -1) const;
  SparseComplexMatrix sumsq (int dim = -1) const;
  SparseMatrix abs (void) const;

  SparseComplexMatrix diag (octave_idx_type k = 0) const;

  // i/o
  //friend std::ostream&
  //operator << (std::ostream& os, const SparseComplexMatrix& a);
  //friend std::istream&
  //operator >> (std::istream& is, SparseComplexMatrix& a);
};