swig -lua -Iinclude src/xtract.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o xtract.so src/xtract_wrap.c -lm -lluajit
