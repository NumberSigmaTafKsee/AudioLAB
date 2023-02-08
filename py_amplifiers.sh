swig -python -c++ -I/usr/local/include/python3.9 -Iinclude Amplifiers.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -I/usr/local/include/python3.9\
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o _Amplifiers.so Amplifiers_wrap.cxx \
-lstdc++ -lm -lpython3.9 -lutil -lrt -ldl
