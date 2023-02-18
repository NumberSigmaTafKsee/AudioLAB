swig -lua -c++ -Iinclude VCA.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o VCA.so \
VCA_wrap.cxx -lstdc++ -lm -lluajit

swig -lua -c++ -Iinclude VCF.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o VCF.so \
VCF_wrap.cxx -lstdc++ -lm -lluajit

swig -lua -c++ -Iinclude VCO.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o VCO.so \
VCO_wrap.cxx -lstdc++ -lm -lluajit