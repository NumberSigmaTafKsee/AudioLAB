swig -lua -c++ src/BiQuad.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o BiQuad.so src/BiQuad_wrap.cxx src/BiQuad.cpp -lstdc++ -lm -lluajit
