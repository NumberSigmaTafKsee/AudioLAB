swig -lua -c++ -Iinclude audiodsp_gamma.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o gamma.so audiodsp_gamma_wrap.cxx lib/libGamma.a -lstdc++ -lm -lportaudio -lluajit
