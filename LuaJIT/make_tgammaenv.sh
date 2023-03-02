swig -lua -c++ -Iinclude src/tgammaenv.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o tGammaenv.so src/tgammaenv_wrap.cxx -lstdc++ -lm -lluajit
