require('cdll')

--[[
cdll.LoadDLL("./cenv.so");
pc = cdll.CreateEnvironment();
cdll.Learn(pc);
cdll.Run(pc);
cdll.CloseEnvironment(pc);
]]--

cdll.LoadDLL("./cenv.so")
f = io.open('test.json')
s = f:read("*a")
pc = cdll.CreateEnvironmentJSON(s);
cdll.Learn(pc);
cdll.Run(pc);
cdll.CloseEnvironment(pc);

