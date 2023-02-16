swig -lua -c++ spuceaudioeq.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spuceaudioeq.so spuceaudioeq_wrap.cxx -lstdc++ -lm -lluajit -lspuce
