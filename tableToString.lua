function table.tostringex (t)
	local tableT = {"{"}

	function toStringF (k, v)
		if "number" == type(k) then
		elseif "string" == type(k) then
			table.insert(tableT, string.format("%s=", k))
		elseif "table" == type(k) then
			table.insert(tableT, string.format("table:%s=", k))
		elseif "function" == type(k) then
			table.insert(tableT, tostring(k))
		end

		if "table" == type(v) then
			table.insert(tableT, table.tostringex(v))
		elseif "number" == type(v) or "boolean" == type(v) then
			table.insert(tableT, tostring(v))
		else 
			table.insert(tableT, string.format("\"%s\"", tostring(v)))
		end

		table.insert(tableT, ",")
	end

	table.foreach(t, toStringF)
	table.remove(tableT, table.getn(tableT)) -- 删除逗号
	table.insert(tableT, "}")
	local tableStr = table.concat(tableT)
	return tableStr
end

-- 第二个参数可选，用于递归之间传值
function table.tostring (t, tabCountedP)
	local tableT = {"{\n"}
	if (nil == tabCountedP) then
		tabCounted = "    "
	else 
		tableT = {"\n", tabCounted, "{\n"}
	end

	function toStringF (k, v)
		table.insert(tableT, tabCounted)
		if "number" == type(k) then
			table.insert(tableT, string.format("[%s] = ", k))
		elseif "string" == type(k) then
			table.insert(tableT, string.format("\"%s\" = ", k))
		elseif "table" == type(k) then
			table.insert(tableT, string.format("table: %s = ", tostring(k)))
		end

		if "table" == type(v) then
			tabCounted = string.format("%s%s", tabCounted, "    ")
			returnStr = table.tostring(v, tabCounted)
			table.insert(tableT, returnStr)
			table.insert(tableT, ",") --分步添加，方便后面删除逗号
			table.insert(tableT, "\n")
		else
			table.insert(tableT, tostring(v))
			table.insert(tableT, ",")
			table.insert(tableT, "\n")
		end
	end

	table.foreach(t, toStringF)
	table.remove(tableT, table.getn(tableT) -1) -- 删除逗号
	tabCounted = string.sub(tabCounted, 1, -5) -- 缩进-1
	table.insert(tableT, string.format("%s}", tabCounted))
	local tableStr = table.concat(tableT)
	return tableStr, tabCounted
end

function stringToTable (stringP)
	local loadFunction = loadstring("return " .. stringP)
	local loadTable = loadFunction()
	return loadTable
end