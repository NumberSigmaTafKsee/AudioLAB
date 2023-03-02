syms z a1 a2 a3
f =  (1+a1*z^-1+a2*z^-2+a3*z^-3)/(1-(0.18*z^-1)+(0.81*z^-2))
[nf, df] = numden(f)
tfdata(f)
