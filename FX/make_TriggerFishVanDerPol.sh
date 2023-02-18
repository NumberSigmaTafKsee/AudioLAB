    swig -lua -c++ -I../include -I../include TriggerFishVanDerPol.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o TriggerFishVanDerPol.so TriggerFishVanDerPol_wrap.cxx -lstdc++ -lm -lluajit 
    