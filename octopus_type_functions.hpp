#pragma once
// I want it to work with Octopus::Functions directly
// this is old
namespace Octopus
{   
///////////////////////////
// MatrixXd
///////////////////////////
    void display(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;        
        l = Functions::octave_display.eval(l,0);          
    }
    OctopusMatrixXf cos(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cos.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf sin(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sin.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf tan(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tan.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf acos(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acos.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf asin(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asin.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf atan(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atan.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf cosh(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cosh.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf sinh(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sinh.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf tanh(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tanh.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf acosh(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acosh.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf asinh(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asinh.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf atanh(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atanh.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf exp(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_exp.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf log(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf log10(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log10.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }
    OctopusMatrixXf pow(const OctopusMatrixXf & a, double b)
    {
        OctopusValueList l;
        l(0) = a;
        l(1) = b;
        l = Functions::octave_power.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }    
    OctopusMatrixXf sqrt(const OctopusMatrixXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sqrt.eval(l,1);  
        return OctopusMatrixXf(l(0).float_matrix_value());
    }    

///////////////////////////
// MatrixXd
///////////////////////////
    void display(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;        
        l = Functions::octave_display.eval(l,0);          
    }
    OctopusMatrixXd cos(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cos.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd sin(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sin.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd tan(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tan.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd acos(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acos.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd asin(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asin.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd atan(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atan.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd cosh(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cosh.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd sinh(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sinh.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd tanh(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tanh.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd acosh(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acosh.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd asinh(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asinh.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd atanh(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atanh.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd exp(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_exp.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd log(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd log10(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log10.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }
    OctopusMatrixXd pow(const OctopusMatrixXd & a, double b)
    {
        OctopusValueList l;
        l(0) = a;
        l(1) = b;
        l = Functions::octave_power.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }    
    OctopusMatrixXd sqrt(const OctopusMatrixXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sqrt.eval(l,1);  
        return OctopusMatrixXd(l(0).float_matrix_value());
    }


///////////////////////////
// RowVectorXf
///////////////////////////
    void display(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;        
        l = Functions::octave_display.eval(l,0);          
    }

    OctopusRowVectorXf cos(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cos.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf sin(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sin.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf tan(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tan.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf acos(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acos.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf asin(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asin.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf atan(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atan.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf cosh(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cosh.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf sinh(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sinh.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf tanh(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tanh.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf acosh(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acosh.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf asinh(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asinh.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf atanh(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atanh.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf exp(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_exp.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf log(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf log10(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log10.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }
    OctopusRowVectorXf pow(const OctopusRowVectorXf & a, double b)
    {
        OctopusValueList l;
        l(0) = a;
        l(1) = b;
        l = Functions::octave_power.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }    
    OctopusRowVectorXf sqrt(const OctopusRowVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sqrt.eval(l,1);  
        return OctopusRowVectorXf(l(0).float_matrix_value());
    }    


///////////////////////////
// RowVectorXd
///////////////////////////
    void display(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;        
        l = Functions::octave_display.eval(l,0);          
    }

    OctopusRowVectorXd cos(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cos.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd sin(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sin.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd tan(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tan.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd acos(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acos.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd asin(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asin.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd atan(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atan.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd cosh(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cosh.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd sinh(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sinh.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd tanh(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tanh.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd acosh(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acosh.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd asinh(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asinh.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd atanh(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atanh.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd exp(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_exp.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd log(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd log10(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log10.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }
    OctopusRowVectorXd pow(const OctopusRowVectorXd & a, double b)
    {
        OctopusValueList l;
        l(0) = a;
        l(1) = b;
        l = Functions::octave_power.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }    
    OctopusRowVectorXd sqrt(const OctopusRowVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sqrt.eval(l,1);  
        return OctopusRowVectorXd(l(0).matrix_value());
    }    

///////////////////////////
// ColVectorXf
///////////////////////////
    void display(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;        
        l = Functions::octave_display.eval(l,0);          
    }

    OctopusColVectorXf cos(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cos.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf sin(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sin.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf tan(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tan.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf acos(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acos.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf asin(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asin.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf atan(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atan.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf cosh(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cosh.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf sinh(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sinh.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf tanh(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tanh.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf acosh(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acosh.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf asinh(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asinh.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf atanh(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atanh.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf exp(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_exp.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf log(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf log10(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log10.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }
    OctopusColVectorXf pow(const OctopusColVectorXf & a, double b)
    {
        OctopusValueList l;
        l(0) = a;
        l(1) = b;
        l = Functions::octave_power.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }    
    OctopusColVectorXf sqrt(const OctopusColVectorXf & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sqrt.eval(l,1);  
        return OctopusColVectorXf(l(0).float_column_vector_value());
    }        

///////////////////////////
// ColVectorXd
///////////////////////////
    void display(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;        
        l = Functions::octave_display.eval(l,0);          
    }

    OctopusColVectorXd cos(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cos.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd sin(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sin.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd tan(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tan.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd acos(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acos.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd asin(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asin.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd atan(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atan.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd cosh(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_cosh.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd sinh(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sinh.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd tanh(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_tanh.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd acosh(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_acosh.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd asinh(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_asinh.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd atanh(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_atanh.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd exp(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_exp.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd log(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd log10(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_log10.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }
    OctopusColVectorXd pow(const OctopusColVectorXd & a, double b)
    {
        OctopusValueList l;
        l(0) = a;
        l(1) = b;
        l = Functions::octave_power.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }    
    OctopusColVectorXd sqrt(const OctopusColVectorXd & a)
    {
        OctopusValueList l;
        l(0) = a;
        l = Functions::octave_sqrt.eval(l,1);  
        return OctopusColVectorXd(l(0).column_vector_value());
    }        


    
    OctopusRowVectorXf chirp(float start, float inc, float end)
    {
        OctopusRowVectorXf v;
        OctopusValueList l,l1;
        l(0) = start;
        l(1) = inc;
        l(2) = end;
        l = Functions::octave_linspace.eval(l,1);
        v = l(0).float_row_vector_value();
        l1(0) = v;
        l = Functions::octave_chirp(l1,1);
        v = l(0).float_row_vector_value();        
        return v;
    }

}
