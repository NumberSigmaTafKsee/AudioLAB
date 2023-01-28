
#include <vector>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h> // fabs

using namespace std;

class PolyFit
{

public: 
	enum SATUS{
		SUCCESS = 0,
		ERROR = 1
	};

	/*
	* Given a set of 2D points, stored in vectors x and y,
	* finds polynomial coeficients of a polynomial of order
	* pol_order that minimizes the sum of square distances,
	* using LUP-decomposition.
	*
	* INPUT:
	*	- 	Vectors x and y, of equal length corresponding
	*			to a 2D dataset.
	*	- 	Integer pol_order corresponding to the desired order
	*	 		of the polynomial.
	* OUTPUT:
	*	-	Integer, 0 if sucess, 1 if failure.
	*	-	A vector poly_coef, of length pol_order+1,
	* 			with estimated polynomial coeficients. The
	*			coeficients are ordered as:
	*			poly_coef[0] + poly_coef[1] x + ...
	*			 + poly_coef[pol_order] x^pol_order.
	* NOTES:
	*	-	Prints message if matrix is close to degenerate
	*/
	static PolyFit::SATUS polyFit(
		vector<double>&outPolCoef, 
		const vector<double> &inX,
		const vector<double>& inY, 
		const unsigned int& inPolOrder
	);

private:

	/* 
	* INPUT: A - vector of vector, corresponding to a square matrix having dimension N
	*        Tol - small tolerance number to detect failure when the matrix is near degenerate
	* OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
	*        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
	*        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
	*        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
	* SOURCE: https://en.wikipedia.org/wiki/LU_decomposition
	*/
	static PolyFit::SATUS LUPDecompose(
		vector<vector<double> >& A,
		const unsigned int& N,
		const double& Tol, 
		vector<unsigned int> &P
	);

	/* Solves the linear system Ax=b, using LUP decomposition.
	*  INPUT: The square matrix A, of size N, and vector b.
	*  OUTPUT: Vector x, the solution of Ax = b.
	*  SOURCE: https://en.wikipedia.org/wiki/LU_decomposition
	*/
	static PolyFit::SATUS LUPSolve(
		vector<double>& x,
		vector<vector<double> >& A,
		const vector<double>& b,
		const unsigned int& N
	);
};


PolyFit::SATUS PolyFit::LUPDecompose(
	vector<vector<double> >& A,
	const unsigned int& N,
	const double& Tol,
	vector<unsigned int> &P
)
{
	unsigned int i, j, k, imax;
	double maxA, absA;

	for (i = 0; i <= N; i++)
		P[i] = i; 

	for (i = 0; i < N; i++) 
	{
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = fabs(A[k][i])) > maxA) 
			{
				maxA = absA;
				imax = k;
			}

		if (maxA < Tol) return SATUS::ERROR; //failure, matrix is degenerate

		if (imax != i) 
		{
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			A[i].swap(A[imax]);
			P[N]++;
		}

		for (j = i + 1; j < N; j++) 
		{
			A[j][i] /= A[i][i];

			for (k = i + 1; k < N; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}
	return SATUS::SUCCESS;  //decomposition done 
}


PolyFit::SATUS PolyFit::LUPSolve(
	vector<double>& x,
	vector<vector<double> >& A,
	const vector<double>& b,
	const unsigned int& N
)
{
	const double tol = 1e-12;
	vector<unsigned int> P(N + 1,0);

	if(LUPDecompose(A, N, tol, P) == SATUS::SUCCESS)
	{
		for (unsigned int i = 0; i < N; i++) 
		{
			x[i] = b[P[i]];

			for (unsigned int k = 0; k < i; k++)
				x[i] -= A[i][k] * x[k];
		}

		for (int i = N - 1; i >= 0; i--) 
		{
			for (unsigned int k = i + 1; k < N; k++)
				x[i] -= A[i][k] * x[k];

			x[i] = x[i] / A[i][i];
		}
	}
	else
	{

		return PolyFit::SATUS::ERROR;
	}
	return PolyFit::SATUS::SUCCESS;
}

PolyFit::SATUS PolyFit::polyFit(
	vector<double>& outPolCoef, 
    const vector<double>& inX,
	const vector<double>& inY, 
	const unsigned int& inPolOrder
)
{
	const unsigned int nX = inX.size();
	const unsigned int nCoef = inPolOrder + 1;
	const unsigned int nPower = nCoef + inPolOrder;

	outPolCoef.resize(nCoef,0.0);

	// The normal matrices.
	// (X' X) * poly_coef = (X' y).
	vector<vector<double> > XX(nCoef, vector<double>(nCoef, 0));
	vector<double > XY(nCoef, 0);
	
	vector<double > xxSum(nPower, 0);

	// Calculating all powers necessary for constructing
	// the normal matrix X' X, and X' y.
	for (unsigned int i = 0; i < nX; i++)
	{
		double accumulatedPower = 1;
		for (unsigned int j = 0; j < nPower; j++)
		{
			if (j < nCoef)
			{
				XY[j] += inY[i] * accumulatedPower;
			}
			xxSum[j] += accumulatedPower;
			accumulatedPower *= inX[i];
		}
	}
	for (unsigned int i = 0; i < nCoef; i++)
	{
		for (unsigned int j = 0; j < nCoef; j++)
		{
			XX[i][j] = xxSum[i + j];
			XX[j][i] = XX[i][j];
		}
	}

	const SATUS errCode = LUPSolve(outPolCoef,XX, XY, nCoef);
	return errCode;
}

using namespace std;
int main() {
    const int order = 3;
    const vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    const vector<double> y = {-7.5740, 5.2055, 11.3065, 9.4567, 7.0056, 1.2991, 20.3567, 52.3649, 105.1265, 196.0492};
    vector<double> polCoef;
    
    const PolyFit::SATUS ec = PolyFit::polyFit(polCoef,x,y,order);
    if(ec == PolyFit::SATUS::SUCCESS)
    {
        cout << "Pol coeficient (c_0 + c_1 * x + c_2 x^2 + ... + c_n * x^n" << endl;
        cout << "c = [  ";
        for(const auto& c : polCoef)
            cout << c << " ";
        cout << " ]" << endl;;
    }



    return 0;
}