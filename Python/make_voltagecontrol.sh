swig -lua -c++ -Iinclude audiodsp_voltagecontrol.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I/usr/local/include/lilv-0 \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o VoltageControl.so audiodsp_voltagecontrol_wrap.cxx \
-lstdc++ -lm -lluajit
