swig -DPYTHON -python -c++ -Iinclude PyPlot.i
gcc -std=c++11 -Wfatal-errors -Iinclude -I/usr/local/include/python3.9 -fpermissive -O2 -fPIC -shared \
-o _PyPlot.so PyPlot_wrap.cxx -lstdc++ -lm -L/usr/local/lib -lpython3.9 -lutil -lrt -ldl
