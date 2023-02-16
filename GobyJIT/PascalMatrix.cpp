#include "Eigen/Eigen"
#include <iostream>

double factorial(double x)
{
    if (x == 0)
        return 1;
    return x * factorial(x - 1);
}
double binomial(double n, double k)
{
    if(n == 0 || k == 0) return 1;
    return factorial(n) / (factorial(k) * factorial(fabs(n - k)));
}

Eigen::MatrixXd LPPascalMatrix(int n)
{
    Eigen::MatrixXd matrix(n,n);
    for(size_t j = 0; j < n; j++)
    {
        matrix(0,j) = 1;
        matrix(n-1,j) = pow(-1,j);
    }  
    
    for(size_t i = 1; i < n; i++)
    {
        matrix(i-1,0) = binomial(n-1,i-1);
        matrix(i-1,n-1) = pow(-1,i+1)*binomial(n-1,i-1);
    } 
    
    for(size_t row = 1; row < n-1; row++)
    {        
        for(size_t i = 1; i < n-1; i++ )
        {
            
            matrix(i,row) = matrix(i,row-1) - matrix(i-1,row-1) - matrix(i-1,row);
        }
    }
    
    return matrix;
}

Eigen::MatrixXd HPPascalMatrix(int n)
{
    Eigen::MatrixXd matrix(n,n);
    for(size_t j = 0; j < n; j++)
    {
        matrix(0,j) = 1;
        matrix(n-1,j) = pow(-1,j);
    }  
    
    for(size_t i = 1; i < n; i++)
    {
        matrix(i-1,n-1) = binomial(n-1,i-1);
        matrix(i-1,0) = pow(-1,i+1)*binomial(n-1,i-1);
    } 
    
    for(size_t row = 1; row < n-1; row++)
    {        
        for(size_t i = 1; i < n-1; i++ )
        {
            
            matrix(i,row) = matrix(i,row-1) + matrix(i-1,row-1) + matrix(i-1,row);
        }
    }
    
    return matrix;
}
// Convert H(s) => H(z)
// it just builds the frequency vector
Eigen::VectorXd convertS(double cf, Eigen::VectorXd & c)
{
    Eigen::VectorXd out(c.rows());        
    out(0) = 1;            
    for(size_t i = 1; i < c.rows(); i++)
    {
        out(i) = c(i)*pow(cf,(double)i);        
    }    
    return out;
}

// convert H(z) => H(s)
Eigen::VectorXd convertD(double cf, Eigen::VectorXd & c)
{
    Eigen::VectorXd out(5);        
    out(0) = 1;        
    for(size_t i = 1; i < c.rows(); i++)
    {
        out(i) = c(i)/pow(cf,(double)i);        
    }
    
    return out;
}
void Lowpass()
{
    Eigen::MatrixXd m = LPPascalMatrix(5);
    
    //std::cout << m << std::endl;
    
    Eigen::VectorXd c(5);
    
    double cf = 1/tan(M_PI*200.0/1000.0);
    //double cf = tan(M_PI*200.0/1000.0);
    c << 1,0,0,0,0;    
    c = convertS(cf,c);
    std::cout << m*c << std::endl;
    
    c << 1,2.6131,3.4142,2.6131,1;
    
    //std::cout << c << std::endl;
    c = convertS(cf,c);

    //c << 1,3.5967,6.4680,6.8136,3.5889;    
    std::cout << m*c << std::endl;
}
void Highpass()
{
    Eigen::MatrixXd m = HPPascalMatrix(5);
    
    //std::cout << m << std::endl;
    
    Eigen::VectorXd c(5);
    
    
    double cf = tan(M_PI*200.0/1000.0);
    c << 1,0,0,0,0;    
    c = convertS(cf,c);
    std::cout << m*c << std::endl;
    
    c << 1,2.6131,3.4142,2.6131,1;       
    c = convertS(cf,c);

    //c << 1,3.5967,6.4680,6.8136,3.5889;    
    std::cout << m*c << std::endl;
    
}


void Bandpass()
{
    double cv = 1/tan(M_PI*3000/10000);
    double tl = tan(M_PI*1000/10000);
    double U = cv / (1-cv*tl);
    double L = tl / (1-cv*tl);
    
    
    int M = 5;

    // we only need to solve biquads not the whole fucking polynomial

    Eigen::MatrixXd HP(M,M);
    HP = HPPascalMatrix(M);
    Eigen::MatrixXd LP(M,M);
    LP = LPPascalMatrix(M);

    // (U+L)^2 = (U+L)(U+L) = U^2 + 2UL + L^2
    Eigen::VectorXd A(M),B(M);
        
    A << 1,1.412,1,0,0;
    
    std::cout << HP*A << std::endl;
    
}
int main()
{    
    Bandpass();
    
}