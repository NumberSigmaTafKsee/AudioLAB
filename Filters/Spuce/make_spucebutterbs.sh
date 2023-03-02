swig -lua -c++ iirspucebutterbs.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spucebutterbs.so iirspucebutterbs_wrap.cxx -lstdc++ -lm -lluajit -lspuce
