kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_analog.i
gcc -fmax-errors=1 -std=c++17 -Iinclude  \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Analog.so $kit/audiodsp_analog_wrap.cxx \
-lstdc++ -lm -lluajit
