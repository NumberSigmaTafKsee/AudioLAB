require('luapm')

f = function(c,d1,d2) print(c,d1,d2) end 

luapm.set_note_on_func(f)
luapm.set_note_off_func(f)

n = luapm.GetNumDevices()
for i=1,n do 
	print(i-1,luapm.GetDeviceName(i-1))
end 

luapm.InitMidi(1,3,3)
while 1 do

end
