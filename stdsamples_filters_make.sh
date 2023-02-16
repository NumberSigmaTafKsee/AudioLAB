kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_filters.i
gcc -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Filters.so $kit/audiodsp_filters_wrap.cxx  \
-lstdc++ -lm -lluajit
