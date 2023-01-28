swig -lua -c++ src/LuaStd.i
gcc -fpermissive -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o LuaStd.so src/LuaStd_wrap.cxx -lstdc++ -lm -lluajit
