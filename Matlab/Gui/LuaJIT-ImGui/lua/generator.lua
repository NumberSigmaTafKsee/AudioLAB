package.path = package.path..";../cimgui/generator/?.lua"
local cpp2ffi = require"cpp2ffi"
local save_data = cpp2ffi.save_data
local read_data = cpp2ffi.read_data
local location = cpp2ffi.location
----utility functions
local function get_cdefs(gccline,locat,cdef)
	cdef = cdef or {}
	local pipe,err = io.popen(gccline,"r")
	if not pipe then error("could not execute gcc "..err) end

	for line in location(pipe,{locat}) do
		line = line:gsub("extern __attribute__%(%(dllexport%)%)%s*","")
		line = line:gsub("extern __declspec%(dllexport%)%s*","")
		if line~="" then table.insert(cdef,line) end
	end
	pipe:close()
	return cdef
end

local function get_all_cdefs(sources)
	local cdefs = {}
	for i,v in ipairs(sources) do
		print("get cdefs from",v)
		cdefs = get_cdefs([[gcc -E -DCIMGUI_DEFINE_ENUMS_AND_STRUCTS -I "../cimgui" ]].."../"..v.."/"..v..".h",v,cdefs)
	end
	return cdefs
end

--------------------------------------------------------
--first get cdefs
print"get cdefs"
local sources = {"cimgui", "cimplot", "cimguizmo", "cimguizmo_quat", "cimnodes","cimnodes_r"}
local cdefs = get_all_cdefs(sources)

print"get cimgui_impl cdefs"
cdefs = get_cdefs([[gcc -E -DCIMGUI_API="" -DCIMGUI_USE_OPENGL3 -DCIMGUI_USE_SDL -DCIMGUI_USE_GLFW -DCIMGUI_USE_OPENGL2 ../cimgui/generator/output/cimgui_impl.h]],"cimgui_impl",cdefs)

table.insert(cdefs,1,"typedef void FILE;")

----- create imgui/cdefs.lua
print"save cdefs.lua"
--time_t is needed for implot_internal
local time_t_cdefs = [=[
--time_t needed for implot_internal
local ffi = require"ffi"
local IS_64_BIT = ffi.abi('64bit')
local cdecl = ""
cdecl = cdecl .. [[
typedef struct tm
{
   int tm_sec   ;
   int tm_min   ;
   int tm_hour  ;
   int tm_mday  ;
   int tm_mon   ;
   int tm_year  ;
   int tm_wday  ;
   int tm_yday  ;
   int tm_isdst ;
} tm;
]]
if ffi.os == "Windows" then
	if IS_64_BIT then
		cdecl = cdecl..[[typedef __int64 time_t;]]
	else
		cdecl = cdecl..[[typedef __int32 time_t;]]
	end
else
	cdecl = cdecl..[[typedef size_t time_t;]]
end

]=]
table.insert(cdefs,1,time_t_cdefs)
table.insert(cdefs,2,"--[[ BEGIN AUTOGENERATED SEGMENT ]]\n cdecl = cdecl .. [[")
table.insert(cdefs,"]]\n--[[ END AUTOGENERATED SEGMENT ]]\n")
local hstrfile = read_data"./imgui_base_cdefs.lua"
save_data("./imgui/cdefs.lua",table.concat(cdefs,"\n"), hstrfile)

----- generate imgui/glfw.lua
print"save glfw.lua"
local class_gen = require"class_gen"
local classes = class_gen(sources)
local iniclass = "local cimguimodule = 'cimgui_glfw' --set imgui directory location\n"
local base = read_data("./imgui_base.lua")
save_data("./imgui/glfw.lua",iniclass, base, classes)

----- generate imgui/sdl.lua
print"save sdl.lua"
local iniclass = "local cimguimodule = 'cimgui_sdl' --set imgui directory location\n"
save_data("./imgui/sdl.lua",iniclass, base, classes)

print"-----------------------------done generation"