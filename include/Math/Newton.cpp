/*
 * Copyright (c) 2020 Samuel Prevost.
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <iostream>
#include <assert.h>
#include <vector>
// Scaling of numeric_limit::epsilon(), if below 12 the resulting epsilon
// tends to be too small for accurate gradient computation, but a scaling too big is
// risky too !
#define EPS_SCALING 18

#include <Eigen/Dense>

// 0 : no verbose except what test is running
// 1 : a bit of verbose
// 2 : print grad and hessian at starting point
#define VERBOSE 1

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;


namespace Newton {
    /**
     * Returns a function which approximates the gradient of func at the inputted point.
     *
     * @param func the function to approximate the gradient of
     * @return a function which approximate the gradient of the input
     */
    function<VectorXd(VectorXd)> grad_apprx(function<double(VectorXd)> func);

    /**
     * Returns a function which approximates the hessian of func at the inputted point.
     *
     * @param func the function to approximate the hessian of
     * @return a function which approximate the hessian of the input
     */
    function<MatrixXd(VectorXd)> hessian_apprx(function<double(VectorXd)> func);

    /**
     * Finds a stationary point of func that is closest to the inputted x.
     * x is modified in place at each iteration.
     *
     * @param x the stating point of the newton method which will update at each iteration.
     * @param func the function which to find a stationnary point of.
     * @param max_it the maximum number of iterations to perform, set to -1 to stop when the gradient is low enough.
     * @return the number of iterations that have been needed to find an acceptable stationary point.
     */
    int newton(VectorXd& x, const function<double(VectorXd)>& func, int max_it);

    /**
     * See int newton(VectorXd& x, function<double(VectorXd)> func) for complete doc
     * Instead of approximating the gradient and hessian functions, you can provide them directly if you
     * know their analytical expression, making the search way faster.
     *
     * @param x the stating point of the newton method which will update at each iteration.
     * @param grad the function used to compute the gradient at each step
     * @param hess the function used to compute the hessian matrix at each step
     * @param max_it the maximum number of iterations to perform, set to -1 to stop when the gradient is low enough.
     * @return the number of iterations that have been needed to find an acceptable stationary point.
     */
    int newton(VectorXd& x, const function<VectorXd(VectorXd)>& grad, const function<MatrixXd(VectorXd)>& hess, int max_it);
}

// Absolute and relative tolerance for equality check
#define ATOL 1e-8
#define RTOL 1e-5

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace Tests {

    struct TestBundle {
        string name;
        function<double(VectorXd)> func;
        int max_it;
        VectorXd starting_x;
        VectorXd target_x;
    };

    /**
     * Checks that |v1 - v2| <= atol + rtol * |v2| element-wise
     * @param v1 vector to compare
     * @param v2 vector to compare
     * @return true if formula is verified for each element, false otherwise.
     */
    bool equality(VectorXd v1, VectorXd v2);

    TestBundle cos_1d_eq();
    TestBundle cos_2d_eq();
    TestBundle square();

    vector<TestBundle> all_bundles();

}

namespace Newton {

    /* Computes only the ith component of the gradient */
    static double grad_apprx(const VectorXd& v, function<double(VectorXd)> func, int i){
        double eps = v(i)*sqrt(numeric_limits<double>::epsilon())*EPS_SCALING;
        if (abs(v(i)) <= numeric_limits<double>::epsilon()){
            eps = numeric_limits<double>::epsilon()*EPS_SCALING;
        }
        VectorXd x_more = v, x_less = v;
        // Add/Sub a diff to the component we want the grad of
        x_more(i) += eps;
        x_less(i) -= eps;
        // Central difference scheme to cancel one or more terms of the Taylor series
        return (func(x_more) - func(x_less))/(2*eps);
    }

    static VectorXd grad_apprx(const VectorXd& v, function<double(VectorXd)> func){
        int size = v.size();
        VectorXd g(size);
        for (int i = 0; i < size; ++i){
            g(i) = grad_apprx(v, func, i);
        }
        return g;
    }

    static MatrixXd hessian_apprx(const VectorXd& v, function<double(VectorXd)> func){
        MatrixXd hess = MatrixXd::Zero(v.size(), v.size());
        // Will compute df/didj for all pairs i j
        for (int i = 0; i < v.size(); ++i) {
            // [&] denote lexical closure hence the use of std::function instead
            // of function pointers (which can't store the extra context)

            // This lambda returns a scalar which is the eval of
            // the grad wrt to the i-th component of func at v2
            auto grad_xi = [&](const VectorXd& v2){
                return grad_apprx(v2, func, i);
            };
            // We then compute the gradient of the gradient of the i-th component along all other dims
            // hence having a vector of df/did1, df/did2, df/did3, df/did4, df/did5...
            // which we ofc put in the ith line of our hessian
            hess.row(i) = grad_apprx(v, grad_xi);
        }
        return hess;
    }

    function<VectorXd(VectorXd)> grad_apprx(function<double(VectorXd)> func){
        return [&](const VectorXd& v){
            return grad_apprx(v, func);
        };
    }

    function<MatrixXd(VectorXd)> hessian_apprx(function<double(VectorXd)> func){
        return [&](const VectorXd& v){
            return hessian_apprx(v, func);
        };
    }

    int newton(VectorXd& x, const function<double(VectorXd)>& func, int max_it){
        return newton(x, grad_apprx(func), hessian_apprx(func), max_it);
    }

    int newton(VectorXd& x, const function<VectorXd(VectorXd)>& grad, const function<MatrixXd(VectorXd)>& hess, int max_it){
        // x = current
        // A = hess
        // X = next
        // b = hess.x - grad(x)
        MatrixXd A;
        VectorXd b, gradVect = grad(x);
        int it = 0;
        if (max_it <= 0){
            max_it = INT_MAX;
        }
        while (gradVect.norm() > numeric_limits<double>::epsilon() && it < max_it){
            /* In theory, the following should work
             * x -= hess(x).inverse() * gradVect;
             * However, inversing the Hessian is often impossible because the determinant is almost zero
             * hence quite unstable (just print out hess(x).determinant()).
             * Using the .inverse() will work when optimising stuff like x^2+y^2, but not on
             * sin(x) - (x^4)/4 where it will just be zero.
             */
            A = hess(x);
            b = A*x - gradVect;
            x = A.colPivHouseholderQr().solve(b);
            gradVect = grad(x);
            it++;
        }
        return it;
    }
}

namespace Tests {

    vector<TestBundle> all_bundles() {
        vector<TestBundle> bundles;

        bundles.push_back(square());
        bundles.push_back(cos_1d_eq());
        bundles.push_back(cos_2d_eq());

        return bundles;
    }

    bool equality(VectorXd v1, VectorXd v2){
        if (v1.size() != v2.size()){
            return false;
        }
        for (int i = 0; i < v1.size(); ++i) {
            if (abs(v1(i)-v2(i)) > ATOL + RTOL * abs(v2(i))){
                return false;
            }
        }
        return true;
    }

    double func_rosenbrock(VectorXd v){
        // v dim is 3
        double x = v(0), y = v(1), z = v(2);
        // Newton cannot optimize this since it has no stationary point
        // but it's still nice to try out
        return 100*(y - x*x) + pow(1 - x, 2)
               + 100*(z - y*y) + pow(1 - y, 2);
    }

    static double func_cos_1d_eq(VectorXd v){
        // v dim is 1
        double x = v(0);
        // From https://www.wikiwand.com/en/Newton%27s_method#/Solution_of_cos(x)_=_x3
        // Models cos(x) - x^3 = 0
        // Integral of cos(x) - x^3 is sin(x) - 1/4 * x^4
        // When the derivative of sin(x) - 1/4 * x^4 = 0, we have
        // cos(x) - x^3 and hence the equation is solved
        return sin(x) - 0.25*x*x*x*x;
    }

    TestBundle cos_1d_eq(){
        VectorXd starting_x(1), target_x(1);
        starting_x << 3;
        target_x << 0.865474;
        return {"cos_1d_eq", func_cos_1d_eq, 20, starting_x, target_x};
    }

    static double func_cos_2d_eq(VectorXd v) {
        // v dim is 2
        double x = v(0), y = v(1);
        // Tried to make the problem a bit more interesting
        // and the solution is non trivial, nice !
        return y * (sin(x) - 0.25 * x * x * x * x);
    }

    TestBundle cos_2d_eq(){
        VectorXd starting_x(2), target_x(2);
        starting_x << 0.87, 10;
        target_x << 0, 0;
        return {"cos_2d_eq", func_cos_2d_eq, 50, starting_x, target_x};
    }

    static double func_square(VectorXd v) {
        // v dim is 2
        double x = v(0), y = v(1);
        // The humble 3D parabola
        return x*x + y*y;
    }

    TestBundle square(){
        VectorXd starting_x(2), target_x(2);
        starting_x << 1000, 1000;
        target_x << 0, 0;
        return {"square", func_square, 10, starting_x, target_x};
    }
}


using namespace Newton;
using namespace Tests;

int main() {

    auto bundles = all_bundles();
    for (TestBundle test : bundles){
        cout << "Running test \"" << test.name << "\"" << endl;
        function<double(VectorXd)> func = test.func;
        VectorXd x = test.starting_x;
        if (VERBOSE) {
            cout << "\tStarting x = [" << x.transpose() << "]" << endl;
            cout << "\tf(x) = " << func(x) << endl;
            if(VERBOSE > 1) {
                cout << "\tgrad(f(x)) = [" << grad_apprx(func)(x).transpose() << "]" << endl;
                cout << "\thess(f(x)) = " << endl << hessian_apprx(func)(x) << endl;
            }
        }

        int it = newton(x, func, max(test.max_it + 10, 1000));

        if (VERBOSE){
            cout << "\tTook " << it << " iterations" << endl;
            cout << "\tNewton result = [" << x.transpose() << "]" << endl;
            cout << "\tFunc = " << func(x) << endl;
            if(VERBOSE > 1) {
                cout << "\tGrad at result = [" << grad_apprx(func)(x).transpose() << "]" << endl;
                cout << "\tHess at result = " << endl << hessian_apprx(func)(x) << endl;
            }
        }

        // We shouldn't exceed the max it count !
        assert(it <= test.max_it);
        // X should be what we expect
        assert(equality(x, test.target_x));
        // Gradient should be null at stationnary point
        assert(equality(grad_apprx(func)(x), VectorXd::Zero(x.size())));
    }
    return 0;
}