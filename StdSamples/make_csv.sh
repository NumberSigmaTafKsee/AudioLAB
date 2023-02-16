swig -Iinclude/Csv -lua -c++ csv.i
gcc -Iinclude/Csv -O2 -fPIC -shared -o csv.so csv_wrap.cxx -lstdc++ -lm -lluajit
