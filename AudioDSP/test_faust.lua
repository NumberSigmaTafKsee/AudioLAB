require('audiodsp_faustfx')
require('plot')

faust = audiodsp_faustfx.FaustFX("ar.dsp")
faust:setControl("attack",0.002)
faust:setControl("release",0.002)
faust:setControl("trigger",1.0)
v1 = audiodsp_faustfx.float_vector(128)
v2 = audiodsp_faustfx.float_vector(128)
faust:Run(128,v1:data(),v1:data())
faust:setControl("trigger",0.0)
faust:Run(128,v2:data(),v2:data())
v = audiodsp_faustfx.float_vector(256)
for i=1,128 do v[i] = v1[i] end
for i=128,256 do v[i] = v2[i-127] end
p = plot.Plot_Float()
p:setstyle("lines")
p:plot_x(v:data(),v:size(),"ar")
os.execute("sleep 16;")


	
