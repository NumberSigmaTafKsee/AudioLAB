import("envelopes.lib");
import("stdfaust.lib");

at = hslider("attack",1.0,0.0,10.0,0.1);
rt = hslider("release",1.0,0.0,10.0,0.1);
t  = button("trigger");
process = en.ar(at,rt,t);
	
