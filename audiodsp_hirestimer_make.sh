swig -lua -c++ -Iinclude audiodsp_hirestimer.i
gcc -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC -shared \
-o audiodsp_hirestimer.so audiodsp_hirestimer_wrap.cxx -lstdc++ -lm -lluajit
