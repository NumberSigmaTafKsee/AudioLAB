#include "octopus.hpp"
#include <Eigen/Core>
#include <initializer_list>


using namespace Octopus;




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
    

/*
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

void proof4()
{
    Matrix a(3,3),b(3,3),r;
    a.fill(1.0);
    b.fill(2.0);
    ValueList v;
    v(0) = a;
    v(1) = b;
    v = interp.eval("hadamard",v,1);
    std::cout << v(0).matrix_value() << std::endl;

    std::string g = "function x = foo(a,b) x=a+b; endfunction";
    interp.interpreter->eval(g.c_str(),0);       
    interp.eval_string("foo(1,2)");
}
void proof()
{
    OctopusValueList ls;
    OctopusMatrixXf  a(3,3),b(3,3),c;
    a.fill(1);
    b.fill(2);
    OctopusFunction func(interp,"foo","function x = foo(a,b) x=a+b; endfunction;");
    ls(0) = a;
    ls(1) = b;
    ls = func(ls,1);
    c = ls(0).float_matrix_value();
    std::cout << c << std::endl;
}
*/

int main(int argc, char * argv[])
{           
    OctopusValueList values,r,v;    
    values.set(0,0.0);
    values.set(1,0.001);
    values.set(2,5.0);
    r = Octopus::octave_linspace(values,1);
    Octopus::octave_display(r);
    v = Octopus::octave_chirp(r,1);
    Octopus::octave_plot(v);
    Octopus::octave_pause(v);    
}