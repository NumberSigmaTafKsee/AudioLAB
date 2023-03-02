#include "Octopus.hpp"
#include <Eigen/Core>
#include <initializer_list>

using namespace Octave;
Octopus interp;              
void sequence(VectorXf & v, int inc)
{
    for(size_t i = 0; i < v.size(1); i++)
        v(i) = i*inc;
}
void sequence(Eigen::VectorXf & v, int inc)
{
    for(size_t i = 0; i < v.rows(); i++)
        v(i) = i*inc;
}

struct value_type
{
    int i;
    float f;
    double d;
    std::string s;
    std::complex<float> cf;
    std::complex<double> cd;
    VectorXf vf;
    VectorXd vd;
    VectorXcf vcf;
    VectorXcd vcd;
    MatrixXf mf;
    MatrixXd md;
    MatrixXcf mcf;
    MatrixXcd mcd;

    enum {
        INT,
        FLOAT,
        DOUBLE,
        COMPLEX_FLOAT,
        COMPLEX_DOUBLE,
        STRING,
        VECTORXF,
        VECTORXD,
        VECTORXCF,
        VECTORXCD,
        MATRIXXF,
        MATRIXXD,
        MATRIXXCF,
        MATRIXXCD,
    };
    int type;

    value_type(int i) {
        this->i = i;
        type    = INT;
    }
    value_type(float f) {
        this->f = f;
        type    = FLOAT;
    }
    value_type(double v) {
        this->d = v;
        type    = DOUBLE;
    }
    value_type(const std::string &s) {
        this->s = s;
        type    = STRING;
    }
};

ValueList arguments(std::initializer_list<value_type> & inputs)
{
    ValueList r;
    int n = 0;
    for(auto i = inputs.begin(); i != inputs.end(); i++)
    {
        switch(i->type)
        {
        case value_type::INT: r(n) = i->i; break;
        case value_type::FLOAT: r(n) = i->f; break;
        case value_type::DOUBLE: r(n) = i->d; break;
        case value_type::COMPLEX_FLOAT: r(n) = i->cf; break;
        case value_type::COMPLEX_DOUBLE: r(n) = i->cd; break;
        case value_type::STRING: r(n) = i->s; break;
        case value_type::VECTORXF: r(n) = i->vf; break;
        case value_type::VECTORXD: r(n) = i->vd; break;
        case value_type::VECTORXCF: r(n) = i->vcf; break;
        case value_type::VECTORXCD: r(n) = i->vcd; break;
        case value_type::MATRIXXF: r(n) = i->mf; break;
        case value_type::MATRIXXD: r(n) = i->md; break;
        case value_type::MATRIXXCF: r(n) = i->mcf; break;
        case value_type::MATRIXXCD: r(n) = i->mcd; break;
        }
        n++;        
    }
    return r;
}
ValueList arguments(size_t n, value_type v1, ... )
{

}

void proof1()
{ 
    Eigen::VectorXf v(16);
    VectorXf  a;
    VectorXcf r;
    sequence(v,1);    
    ValueList l,l1;
    l(0) = Octavate(v);
    l = fft(l);
    r = l(0).float_complex_row_vector_value();
    plot(l);        
    pause();
    
    l(0) = r;
    l = ifft(l);
    plot(l);    
    pause();

    a = l(0).float_row_vector_value();  
}
void proof2()
{ 
    ValueList l;
    Eigen::MatrixXf a(3,3),b(3,3);
    a.fill(1);
    b.fill(2);
    MatrixXf x1 = Octavate(a);
    MatrixXf x2 = Octavate(b);    
    std::cout << x1*x2 << std::endl;
}
void proof3()
{
    ValueList values;
    values(0) = 4;
    values(1) = 1.0;
    values(2) = 's';
    values = Octave::butter(values,2);
    values = Octave::tf2sos(values,2);
    MatrixXd m = values(0).matrix_value();
    std::cout << m << std::endl;    
}

int main(int argc, char * argv[])
{       
    Matrix a(3,3),b(3,3),r;
    a.fill(1.0);
    b.fill(2.0);
    ValueList v;
    v(0) = a;
    v(1) = b;
    v = interp.eval("hadamard",v,1);
    std::cout << v(0).matrix_value() << std::endl;
}