swig -lua -c++ expr.i
gcc -O2 -fPIC -shared -o expr.so expr.cpp expr_wrap.cxx -lstdc++ -lm -lluajit
