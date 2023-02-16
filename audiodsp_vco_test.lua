require('analog_blit_oscillators')
require('stdsamples')
require('plot')
require('audiodsp_hirestimer')
timer = audiodsp_hirestimer.HiResTimer()
x = analog_blit_oscillators.Blit2SawOscillator_f64(44100)
m = stdsamples.sample_vector(256)
timer:Start()
for j=1,10 do
    x:ProcessSIMD(256,nil,m:data())
end
time1 = timer:Stop()
print("time1=",time1,"/10=",time1/10)

