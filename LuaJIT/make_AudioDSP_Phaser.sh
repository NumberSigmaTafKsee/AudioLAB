    swig -lua -c++ -Iinclude AudioDSP_Phaser.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_Phaser.so AudioDSP_Phaser_wrap.cxx -lstdc++ -lm -lluajit    
    