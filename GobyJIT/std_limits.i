%{
#include <limits>
%}

namespace std {
    template<class T>
    struct numeric_limits
    {
        //bool is_specialized();
        T min();
        T max();
        //int digits();
        //int digits10();
        //bool is_signed();
        //bool is_integer();
        //bool is_exact();
        //int radix();
        T epsilon();
        T round_error();
        //int min_exponent();
        //int min_exponent10();
        //int max_exponent();
        //int max_exponent10();
        //bool has_infinity();
        //bool has_quiet_NaN();
        //bool has_signaling_NaN();
        //bool has_denorm_loss();
        T infinity();
        T quiet_NaN();
        T signaling_NaN();
        T denorm_min();
        //bool is_iec5599();
        //bool is_bounded();
        //bool is_modulo();
        //bool traps();
        //bool tinyness_before();
    };
}