#pragma once

#include <vector>
#include <cmath>
#include <iostream>

class Utils {

public:
    static double array_sum(double arr[], int len);

    static double *array_pow(double arr[], int len, int power);

    static double *array_multiplication(double arr1[], double arr2[], int len);

    static double *array_diff(double arr1[], double arr2[], int len);

};



class LinearRegression {

public:

    // First feature
    double *x;

    // Target feature
    double *y;

    // Number of training examples
    int m;

    // The theta coefficients
    double *theta;

    /**
     * Create a new instance from the given data set.
     */
    LinearRegression(double x[], double y[], int m);

    /**
     * Train the model with the supplied parameters.
     *
     * @param alpha         The learning rate, e.g. 0.01.
     * @param iterations    The number of gradient descent steps to do.
     */
    void train(double alpha, int iterations);

    /**
     * Try to predict y, given an x.
     */
    double predict(double x);

private:

    /**
     * Compute the cost J.
     */
    static double compute_cost(double x[], double y[], double theta[], int m);

    /**
     * Compute the hypothesis.
     */
    static double h(double x, double theta[]);

    /**
     * Calculate the target feature from the other ones.
     */
    static double *calculate_predictions(double x[], double theta[], int m);

    /**
     * Performs gradient descent to learn theta by taking num_items gradient steps with learning rate alpha.
     */
    static double *gradient_descent(double x[], double y[], double alpha, int iters, double *J, int m);

};


double Utils::array_sum(double arr[], int len) {
    double s = 0;

    for (int i = 0; i < len; ++i) {
        s += arr[i];
    }

    return s;
}

double *Utils::array_pow(double arr[], int len, int power) {
    double *arr2 = new double[len];

    for (int i = 0; i < len; ++i) {
        arr2[i] = pow(arr[i], power);
    }

    return arr2;
}

double *Utils::array_multiplication(double arr1[], double arr2[], int len) {
    double *arr = new double[len];

    for (int i = 0; i < len; ++i) {
        arr[i] = arr1[i] * arr2[i];
    }

    return arr;
}

double *Utils::array_diff(double arr1[], double arr2[], int len) {
    double *arr = new double[len];

    for (int i = 0; i < len; ++i) {
        arr[i] = arr1[i] - arr2[i];
    }

    return arr;
}


using namespace std;

LinearRegression::LinearRegression(double x[], double y[], int m) {
    this->x = x;
    this->y = y;
    this->m = m;
}

void LinearRegression::train(double alpha, int iterations) {
    double *J = new double[iterations];
    this->theta = gradient_descent(x, y, alpha, iterations, J, m);

    cout << "J = ";
    for (int i = 0; i < iterations; ++i) {
        cout << J[i] << ' ';
    }
    cout << endl << "Theta: " << theta[0] << " " << theta[1] << endl;
}

double LinearRegression::predict(double x) {
    return h(x, theta);
}

double LinearRegression::compute_cost(double x[], double y[], double theta[], int m) {
    double *predictions = calculate_predictions(x, theta, m);
    double *diff = Utils::array_diff(predictions, y, m);
    double *sq_errors = Utils::array_pow(diff, m, 2);
    return (1.0 / (2 * m)) * Utils::array_sum(sq_errors, m);
}

double LinearRegression::h(double x, double theta[]) {
    return theta[0] + theta[1] * x;
}

double *LinearRegression::calculate_predictions(double x[], double theta[], int m) {
    double *predictions = new double[m];

    // calculate h for each training example
    for (int i = 0; i < m; ++i) {
        predictions[i] = h(x[i], theta);
    }

    return predictions;
}

double *LinearRegression::gradient_descent(double x[], double y[], double alpha, int iters, double *J, int m) {
    double *theta = new double[2];
    theta[0] = 1;
    theta[1] = 1;

    for (int i = 0; i < iters; ++i) {
        double *predictions = calculate_predictions(x, theta, m);
        double *diff = Utils::array_diff(predictions, y, m);

        double *errors_x1 = diff;
        double *errors_x2 = Utils::array_multiplication(diff, x, m);

        theta[0] = theta[0] - alpha * (1.0 / m) * Utils::array_sum(errors_x1, m);
        theta[1] = theta[1] - alpha * (1.0 / m) * Utils::array_sum(errors_x2, m);

        J[i] = compute_cost(x, y, theta, m);
    }

    return theta;
}