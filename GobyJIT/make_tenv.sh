swig -lua -c++ -Iinclude src/TEnvelope.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o TEnvelope.so src/TEnvelope_wrap.cxx -lstdc++ -lm -lluajit
