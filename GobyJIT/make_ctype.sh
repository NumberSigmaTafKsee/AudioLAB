swig -lua -c++ CType.i
gcc -Iexpr -O2 -fPIC -shared -o CType.so CType_wrap.cxx expr.cpp -lstdc++ -lm -lluajit
