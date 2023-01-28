require('vf')
a = vf.Matrix(128,128)
t = os.clock()
for i=1,100000 do 
r = a*a
--collectgarbage()
end 
print(os.clock()-t)
