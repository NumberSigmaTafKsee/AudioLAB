#include <iostream>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

struct Point
{
	double x, y;
	Point(double new_x, double new_y)
	{
		x = new_x;
		y = new_y;
	}
};

vector<float> GaussianElimination(vector< vector<float> > A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k = i; k < n + 1; k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; k++) {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n + 1; j++) {
                if (i == j) {
                    A[k][j] = 0;
                }
                else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<float> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = A[i][n] / A[i][i];
        for (int k = i - 1; k >= 0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    cout << "Coefficients: " << endl;
    for (int i = 0; i < A.size(); i++)
        cout << "A " << i << " = " << x[i] << endl;
    return x;
}

double Approx(double (*basis)(double x, double degree), float x, vector<float>& a) {
	double result = 0;
	for (int i = 0; i < a.size(); i++)
		result += basis(x, i) * a[i];
	return result;
}

vector<Point> LeastSquaresPolynomialApprox(const vector<Point>& XY, unsigned short degree, float step) {
	vector<vector<float>> AB(degree + 1);
	for (int i = 0; i < degree + 1; i++) {//row iterator
		AB[i].resize(AB.size() + 1);
		for (int k = 0; k < degree + 1; k++)//column iterator, except B col
			for (int j = 0; j < XY.size(); j++)//point iterator
				AB[i][k] += pow(XY[j].x, i) * pow(XY[j].x, k);
		for (int k = 0; k < XY.size(); k++)
			AB[i][AB.size()] += XY[k].y * pow(XY[k].x, i);

	}
	auto A = GaussianElimination(AB);// A coeff
	//make an approximation
	vector<Point> ps;
	ps.reserve((XY[XY.size() - 1].x - XY[0].x) / step + 1.5);
	int j = 0;
    
	for (double i = XY[0].x; i <= XY[XY.size() - 1].x; i += step, j++)
		ps.push_back({ i,Approx(pow,i,A) }); 

   // for (auto it = XY.begin(); it < XY.end() - 1; it++)
    //    ps.push_back({ (*it).x, Approx(pow, (*it).x, A) });

	ps.push_back({ XY[XY.size() - 1].x,Approx(pow,XY[XY.size() - 1].x,A) });

    cout << endl << "Linear system: ";
    for (int i = 0; i < degree + 1; i++)
    {
        cout << endl;
        for (int k = 0; k < degree + 2; k++)
            cout << AB[i][k] << " ";
    }
    cout << endl << endl;
	return ps;
}

void renormalizationX(vector<Point>& XY, float left, float right) {
    for (auto it = XY.begin(); it < XY.end(); it++) {
        (*it).x = 2 * (((*it).x - left) / (right - left)) - 1;
    }
}

void normalizationX(vector<Point>& XY, float left, float right) {
    for (auto it = XY.begin(); it < XY.end(); it++)
    {
        (*it).x = (right + left) / 2.0 + (right - left) / 2.0 * (*it).x;
    }
}

double ChebyshevPolynomial(double x, unsigned int degree) {
    if (degree == 0) return 1.0;
    if (degree == 1) return x;
    return 2.0 * x * ChebyshevPolynomial(x, degree - 1) - ChebyshevPolynomial(x, degree - 2);
}

double ChebyshevApprox(float x, vector<float>& a) {
    double result = 0;
    for (unsigned int i = 0; i < a.size(); i++)
        result += ChebyshevPolynomial(x, i) * a[i];
    return result;
}

vector<Point> chebyshevPolyAppr(vector<Point>& XY, unsigned short degree, float step) {
    vector<float> A(degree);
    float left = XY[0].x, right = XY[XY.size() - 1].x;
    renormalizationX(XY, left, right);
    //step = 2 * ((step - left) / (right - left)) - 1;
    cout << endl;
    for (unsigned int i = 0; i < degree; i++) {
        float numerator = 0, denominator = 0;
        for (unsigned int j = 0; j < A.size(); j++) {
            float T = ChebyshevPolynomial(XY[j].x, i);     
            numerator += XY[j].y * T;
            denominator += T * T;
        }
        A[i] = numerator / denominator;
    }

    cout << endl << "Coefficients: " << endl;
    for (int i = 0; i < A.size(); i++)
        cout << "A " << i << " = " << A[i] << endl;
    cout << " " << endl;

    vector<Point> ps;
    ps.reserve((XY[XY.size() - 1].x - XY[0].x) / step + 1.5);
    int j = 0;
    double temp = 0;
    for (double i = XY[0].x; i <= XY[XY.size() - 1].x; i += step, j++)
        ps.push_back({ i, ChebyshevApprox(i, A) });
    //for (auto it = XY.begin(); it < XY.end() - 1; it++)
     //   ps.push_back({ (*it).x, ChebyshevApprox((*it).x, A) });

    ps.push_back({ 1,ChebyshevApprox(1, A) });
    normalizationX(ps, left, right);
    normalizationX(XY, left, right);
    return ps;
}
