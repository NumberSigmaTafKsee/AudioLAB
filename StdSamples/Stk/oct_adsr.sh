swig -Iinclude/stk -octave -c++ stk_adsr.i
gcc -fmax-errors=1 -Iinclude/stk -I/usr/local/include/octave-7.3.0 -O2 -fPIC -march=native -mavx2 -shared -o oct_adsr.so stk_adsr_wrap.cxx lib/libstk.a -lstdc++ -lm -L/usr/local/lib/octave/7.3.0 -loctave -loctinterp
mkoctfile -Iinclude/stk -o adsr stk_adsr_wrap.cxx -Llib -lstk
