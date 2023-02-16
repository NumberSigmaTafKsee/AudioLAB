swig -Iinclude/Csv -lua -c++ src/csv.i
gcc -Iinclude/Csv -O2 -fPIC -shared -o csv.so src/csv_wrap.cxx -lstdc++ -lm -lluajit
