[z,p,k] = butter(9,300/500,'high');
sos = zp2sos(z,p,k);
disp(sos);
