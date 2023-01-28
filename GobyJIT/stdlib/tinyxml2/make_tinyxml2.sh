swig -lua -c++ tinyxml2.i
gcc -O2 -fPIC -shared -o tinyxml2.so tinyxml2_wrap.cxx tinyxml2.cpp -lstdc++ -lm -lluajit
