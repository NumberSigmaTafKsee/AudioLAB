swig -lua -c++ spuceremezfir.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spuceremezfir.so spuceremezfir_wrap.cxx -lstdc++ -lm -lluajit -lspuce
