% abg.m - dB to absolute amplitude units 
%
% usage: G = abg(Gdb)
%
% computes G = 10^(Gdb/20), for a vector Gdb

function G = abg(Gdb)

G = 10.^(Gdb/20);



