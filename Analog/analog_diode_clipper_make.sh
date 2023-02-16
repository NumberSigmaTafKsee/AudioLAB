swig -lua -c++ -Iinclude analog_diode_clipper.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o analog_diode_clipper.so \
analog_diode_clipper_wrap.cxx -lstdc++ -lm -lluajit
