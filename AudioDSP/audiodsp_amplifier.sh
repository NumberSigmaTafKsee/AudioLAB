swig -lua -c++ -Iinclude audiodsp_amplifier.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o audiodsp_amplifier.so audiodsp_amplifier_wrap.cxx \
-lstdc++ -lm -lluajit
