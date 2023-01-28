swig -lua -c++ -Iinclude audiodsp_delay.i
gcc -Iinclude -fmax-errors=1 -std=c++17 -I. \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Delay.so audiodsp_delay_wrap.cxx  \
-lstdc++ -lm -lluajit
