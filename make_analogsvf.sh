gcc -std=c++17 -fmax-errors=1 -O2 -march=native -mavx2 -mfma -Iinclude -I/usr/local/include -I/usr/local/include/luajit-2.1 \
-fPIC -shared -o plugins/AnalogSVF.lv2/AnalogSVF.so AnalogSVF.cpp -lstdc++ -lm -llv2-gui -llv2-plugin
