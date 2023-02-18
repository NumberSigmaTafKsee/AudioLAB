swig -lua -c++ -Iinclude Kits/AudioDSP/audiodsp_voltagecontrol.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -IKits/AudioDSP \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o VoltageControl.so Kits/AudioDSP/audiodsp_voltagecontrol_wrap.cxx \
-lstdc++ -lm -lluajit
