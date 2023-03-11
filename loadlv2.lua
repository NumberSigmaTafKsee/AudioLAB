require('audiosystem')
require('lv2plugin')
lv2plugins = lv2plugin.LV2Plugins()
lv2luajit = lv2plugins:LoadPlugin(arg[1])

numpress =0
function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

function noise(input,output,frames)     	
    lv2luajit:ProcessBlock(frames,output,output)
    return 0
end

function note_on(c,n,v)    
    local f = math.pow(2.0, (n-69)/12)*440.0                     
end
function note_off(c,n,v)        
    numpress = numpress-1
    if(numpress <= 0) then 
        numpress = 0;        
    end
    print(numpress)
end

function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 102) then        
        fc = d2/127                
        lv2luajit:setControl(2,fc*20000)
    elseif(d1 == 103) then
        q = d2/127                
        lv2luajit:setControl(3,q*10.0)
    end
end



audiosystem.Init()
audiosystem.Pm_Initialize()


audiosystem.set_note_on_func(note_on)
audiosystem.set_note_off_func(note_off)
audiosystem.set_control_change_func(control)


for i=0,audiosystem.GetNumMidiDevices()-1 do 
    print(i,audiosystem.GetMidiDeviceName(i))
end

audiosystem.set_audio_func(noise)

device=14
audiosystem.Pa_Initialize()
for i=0,audiosystem.GetNumAudioDevices()-1 do 
    if( audiosystem.GetAudioDeviceName(i) == 'jack') then        
        device = i         
    end    
    print(i,audiosystem.GetAudioDeviceName(i))
end

audiosystem.InitMidiDevice(1,3,3)
audiosystem.InitAudioDevice(device,-1,1,44100,256)
audiosystem.RunAudio()
audiosystem.StopAudio()
