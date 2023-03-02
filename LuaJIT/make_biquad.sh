swig -lua biquad.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o biquad.so biquad_wrap.c biquad.c -lm -lluajit
