require('luapa')
require('luapm')
require('LuaPolyBLEP')
require('tauadsr')
require('bitcrush')
require('bigverb')

osc1 = LuaPolyBLEP.PolyBLEP(44100)
osc2 = LuaPolyBLEP.PolyBLEP(44100)
adsr = tauadsr.TauADSR_New(44100)
gate = 0
bit  = bitcrush.BitCrush_New(44100)
bverb = bigverb.BigVerb_New(44100)


function noise(input,output,frames)
    luapm.LockMidi()
    for i = 0,frames-1 do
        local x = osc1:getAndInc() + osc2:getAndInc()
        x = x * 0.5
        local a = tauadsr.TauADSR_Tick(adsr,gate)           
        x = bitcrush.BitCrush_Tick(bit,x * a)
        bigverb.BigVerb_Tick(bverb,x,x)        
        luapa.float_set(output,i,bverb.bigverb.out1)
        luapa.float_set(output,i+1,bverb.bigverb.out2)                            
        luapa.float_set(output,2*i,x)
        luapa.float_set(output,2*i+1,x)
    end   
    luapm.UnlockMidi()
end 
function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end
function note_on(c,n,v)    
    local f = math.pow(2.0, (n-69)/12)*440.0                
    osc1:setFrequency(f)
    gate = 1
end
function note_off(c,n,v)
    gate = 0    
end
-- app callback, midi handling and logic
-- isAudioRunning shuld be atomic
-- either block audio or wait until finished
-- this is run every 10ms, or can be changed in portaudio.i
function callback()
    
end 

luapm.set_note_on_func(note_on)
luapm.set_note_off_func(note_off)
luapm.InitMidi(1,3,3)
for i=0,luapm.GetNumDevices()-1 do 
    print(i,luapm.GetDeviceName(i))
end
luapa.set_audio_func(noise)
--luapa.set_callback_func(callback)
luapa.Pa_Initialize()
for i=0,luapa.GetNumDevices()-1 do 
    print(i,luapa.GetDeviceName(i))
end
-- use jack on my system it is device 6
-- use pulse on my system it is device 14
-- no input is used (-1)
luapa.InitAudioDevice(10,-1,2,44100,64)
luapa.RunAudio()
luapa.StopAudio()