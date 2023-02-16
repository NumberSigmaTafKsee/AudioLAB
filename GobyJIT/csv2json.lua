#!/usr/bin/env th
local json = require 'cjson'

local opt = lapp[[
Convert a CSV file to JSON.
Output paths ending in / are interpreted as a directories. In this case, each row
will become a separate file in the output directory. Otherwise, the output is a
single JSON file containing the array of objects.
This script buffers input/output and thus should be safe on large files.
   <input>  (string)           The input CSV file
   <output> (optional string)  Path to JSON output, stdout by default.
   --map    (optional string)  Path to lua file with mapping function (see README)
   --no-header                 No CSV column headers; treat rows as unlabeled array
]]

local CHUNKSIZE = 2^13 -- 8kb
local dirOut = opt.output and
   (stringx.endswith(opt.output, '/') or path.isdir(opt.output))
local fileOut = not dirOut
local mapFunc = opt.map and require(opt.map) or (function(d) return d end)
local headers -- used to hold the column headers
local lineIdx = 0


if fileOut then
   opt.output = opt.output and
      io.open(opt.output, 'w') or io.output(opt.output)
   opt.output:write('[\n')
end

if dirOut then
   opt.output = path.abspath(path.expanduser(opt.output))
   os.execute('mkdir -p "' .. opt.output .. '"')
end

opt.input = io.open(opt.input, 'r')

local readChunk = function()
   local lines, rest = opt.input:read(CHUNKSIZE, '*line')
   if not lines then return nil end
   if rest then lines = lines .. rest end
   lines = stringx.strip(lines, '\n')
   return stringx.split(lines, '\n')
end

local splitLine = function(line)
   return stringx.split(stringx.strip(line), ',')
end

local decodeLine = function(line, idx)
   local cols = splitLine(line)
   local res = cols
   if headers then
      res = {}
      for i,val in ipairs(cols) do
         local k = headers[i] or '__UNKNOWN_HEADER__'..i
         if val ~= '' then
            res[k] = val
         end
      end
   end
   return mapFunc(res, idx)
end

local writeData = function(data, idx)
   if not data then return end
   if dirOut then
      local opath = path.join(opt.output,
         (data.__filename or tostring(idx)) .. '.json')
      data.__filename = nil
      local ofile = io.open(opath, 'w')
      ofile:write(json.encode(data))
      ofile:close()
   else
      if idx > 0 then opt.output:write(',\n') end
      opt.output:write(json.encode(data))
   end
end

while true do
   local lines = readChunk()
   if not lines then break end

   for _,line in ipairs(lines) do
      if lineIdx == 0 and not opt['no-header'] and not headers then
         -- Grab the headers from the first line instead of writing
         headers = splitLine(line)
      else
         local data = decodeLine(line, lineIdx)
         writeData(data, lineIdx)
         lineIdx = lineIdx + 1
      end
   end
end

if fileOut then
   opt.output:write('\n]')
   opt.output:close()
end