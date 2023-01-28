require('luapa')
err = luapa.Pa_Initialize()
print(err)
n = luapa.GetNumDevices()
print(n)
for i=1,n do 
	print(luapa.GetDeviceName(i-1))
end
