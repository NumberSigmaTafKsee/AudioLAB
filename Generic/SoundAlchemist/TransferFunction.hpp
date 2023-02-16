

struct TransferFunction
{
    std::vector<double> num,den;
    std::vector<double> X,Y;
    double x,y;

    TransferFunction(size_t n, size_t d) {
        num.resize(n);
        den.resize(n);
        X.resize(n);
        Y.resize(d);
    }

    void setNumCoefficient(size_t i, double c) {
        num[i] = c;
    }
    void setDenCoefficient(size_t i, double c) {
        den[i] = c;
    }

    double Tick(double I)
    {
        x = I;
        y = num[0]*x;
        for(size_t i = 1; i < num.size(); i++)
            y += num[i] * X[i];
        for(size_t i = 0; i < den.size(); i++)
            y -= den[i] * Y[i];            
        for(size_t i = 1; i < num.size(); i++)
            X[i] = X[i-1];
        for(size_t i = 1; i < den.size(); i++)
            Y[i] = Y[i-1];
        X[0] = x;
        Y[0] = y;
        return y;
    }
};

