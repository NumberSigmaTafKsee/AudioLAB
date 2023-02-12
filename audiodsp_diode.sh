kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude -Iinclude/DSPFilters audiodsp_diode.i
gcc -fmax-errors=1 -std=c++17 -Iinclude/DSPFilters -Iinclude -O2 -fPIC \
-march=native -mavx2 -shared -o audiodsp_diode.so audiodsp_diode_wrap.cxx -lstdc++ -lm -lluajit 
