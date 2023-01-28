swig -python -c++ -I/usr/include/aubio aubio.i
gcc -O2 -fPIC -march=native -mavx2 -shared -I/usr/local/include/python3.9 -I/usr/include/aubio -o  _aubio.so aubio_wrap.cxx -lstdc++ -lm -lpython3.9 -laubio

