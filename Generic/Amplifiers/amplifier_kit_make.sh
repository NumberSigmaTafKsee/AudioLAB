swig -lua -c++ -Iinclude Kits/Analog/amplifier_kit.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o amplifier_kit.so Kits/Analog/amplifier_kit_wrap.cxx \
-lstdc++ -lm -lluajit
