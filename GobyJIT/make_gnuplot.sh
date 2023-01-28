swig -lua -c++ -IGNUPlot src/GNUPlot.i
gcc -O2 -IGNUPlot -fPIC -march=native -mavx2 -shared -o GNUPlot.so src/GNUPlot_wrap.cxx src/GNUPlot.cpp -lstdc++ -lm -lluajit
