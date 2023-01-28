require('se')
a = se.FloatMatrix(10,10)
for i=0,9 do
for j=0,9 do
a[i][j]=i*10 + j 
end
end
r = a*a
r:print()
