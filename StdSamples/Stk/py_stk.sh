swig -Iinclude/stk -python -c++ stk.i 
gcc -fmax-errors=1 -Iinclude/stk -I/usr/local/include/octave-7.3.0 -I/opt/intel/oneapi/intelpython/python3.9/include/python3.9 \
-O2 -fPIC -march=native -mavx2 -shared -o _stk.so stk_wrap.cxx lib/libstk.a \
-lstdc++ -lm -L/opt/intel/oneapi/intelpython/python3.9/lib -lpython3.9 -lpulse -lpulse-simple -lasound -ljack

