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


class boolMatrix
{
public:

  boolMatrix (void);
  boolMatrix (const boolMatrix& a);
  boolMatrix& operator = (const boolMatrix& a);
  ~boolMatrix (void);

  boolMatrix (octave_idx_type r, octave_idx_type c);
  boolMatrix (octave_idx_type r, octave_idx_type c, bool val);
  boolMatrix (const dim_vector& dv) : boolNDArray (dv.redim (2));
  boolMatrix (const dim_vector& dv, bool val);
  boolMatrix (const Array<bool>& a) : boolNDArray (a.as_matrix ());

  bool operator == (const boolMatrix& a) const;
  bool operator != (const boolMatrix& a) const;

  boolMatrix transpose (void) const;

  // destructive insert/delete/reorder operations

  boolMatrix&
  insert (const boolMatrix& a, octave_idx_type r, octave_idx_type c);

  // unary operations

  boolMatrix operator ! (void) const;

  // other operations

  boolMatrix diag (octave_idx_type k = 0) const;


  void resize (octave_idx_type nr, octave_idx_type nc, bool rfv = false);
};