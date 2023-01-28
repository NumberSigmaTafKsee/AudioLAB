t = 2*pi*linspace(0,1,1024)';
y = sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3);

data1 = abs(y); # Positive values
[pks idx] = findpeaks(data1);

data2 = y; # Double-sided
[pks2 idx2] = findpeaks(data2,"DoubleSided");
[pks3 idx3] = findpeaks(data2,"DoubleSided","MinPeakHeight",0.5);

subplot(1,2,1)
plot(t,data1,t(idx),data1(idx),'xm')
axis tight
subplot(1,2,2)
plot(t,data2,t(idx2),data2(idx2),"xm;>2*std;",t(idx3),data2(idx3),"or;>0.1;")
axis tight
legend("Location","NorthOutside","Orientation","horizontal")
