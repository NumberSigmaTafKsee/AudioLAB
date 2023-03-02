gcc -Wfatal-errors -I/usr/local/include/octave-7.3.0  -o octopus octopus.cpp octopus_wrap.cxx \
-lstdc++ -lm -lluajit -L/usr/local/lib/octave/7.3.0 -loctave -loctinterp
