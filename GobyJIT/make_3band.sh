swig -lua src/3band.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o 3band.so src/3band_wrap.c src/3band.c -lm -lluajit
