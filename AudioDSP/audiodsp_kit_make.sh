kit="Kits/Analog"
swig -lua -c++ -Iinclude $kit/analog_kit.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o analog_kit.so $kit/analog_kit_wrap.cxx \
-lstdc++ -lm -lluajit
