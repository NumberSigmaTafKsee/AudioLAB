swig -Iinclude -lua -c++ re2.i
gcc -Iinclude -O2 -fPIC -shared -o re2.so re2_wrap.cxx lib/libre2.a -lstdc++ -lm -lluajit -lpthread
