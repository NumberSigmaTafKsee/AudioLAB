    swig -lua -c++ -Iinclude analog_vcf.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_vcf.so analog_vcf_wrap.cxx -lstdc++ -lm -lluajit 
    