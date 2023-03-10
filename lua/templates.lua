module('templates', package.seeall)

-- TODO: maybe automate this?
local translate = dofile(template_file)

--- Creates a deep copy of an object.
local function deepcopy(object)
	local lookup_table = {}
	local function _copy(object)
		if type(object) ~= "table" then
			return object
		elseif lookup_table[object] then
			return lookup_table[object]
		end
		local new_table = {}
		lookup_table[object] = new_table
		for index, value in pairs(object) do
			-- HACK: generate new ids for copied nodes
			if index == "id" then new_table.id = next_id()
			-- deep copy value
			else new_table[_copy(index)] = _copy(value) end
		end
		return new_table
	end
	return _copy(object)
end

-- cannot modify idindex directly while traversing it ->
-- new methods from template classes are added here first and then
-- added to idindex after the traversal
local idindex_add = {}

function contains(name)
	for _,t in pairs(translate[module_name] or {}) do
		for _,decl in pairs(t) do
			if decl == name then
				return true
			end
		end
	end
	return false
end

--- Return true if an instance of templated class should be created.
function should_copy(class)
	return translate[class.xarg.fullname] or
		(translate[module_name] and translate[module_name][class.xarg.fullname])
end

--- Creates instantiated copies of template class.
-- New classes are created according to the 'translate' table as deep copies, and
-- are inserted into the 'ret' table.
function create(class, ret)
	local temps = should_copy(class)

	local replace_in = {
		name = true,
		context = true,
		fullname = true,
		member_of = true,
		member_of_class = true,
		scope = true,
		type_base = true,
		type_name = true,
		return_type = true,
		defaultvalue = true,
	}

	local function is_disabled_func(name, list)
		for _,n in ipairs(list) do
			if name == n then
				return true
			end
		end
		return false
	end

	local function template_repare(obj, orig, new, disable_funcs)
		-- generate replace-map splite by ','
		--	orig: Key,T
		--	new: QString,QVariant
		--  	-> { Key = 'QString', T = 'QVariant' }
		local replaces = {}
		for o in orig:gmatch('[^,]+') do table.insert(replaces, o) end
		for n in new:gmatch('[^,]+') do
			local o = table.remove(replaces, 1)
			assert(o ~= nil, string.format('invalid replace tempalte: %s -> %s', orig, new))
			replaces[o] = n
		end

		for k,v in pairs(obj) do
			if replace_in[k] then
				-- local old = obj[k]

				obj[k] = obj[k]:gsub('%w+', function(o)
					return replaces[o] or o
				end)

				-- if old ~= obj[k] then
				-- 	print(' >> replace', orig, new, old, obj[k])
				-- end
			elseif k == 'member_template_parameters' then
				-- ignore
				obj[k] = nil
			elseif type(v) == "table" then
				template_repare(v, orig, new, disable_funcs)
			end
		end
		if obj.label and obj.label:match'^Function' and not is_disabled_func(obj.xarg.name, disable_funcs) then
			idindex_add[obj] = true -- will be copied to index, so that later it can be picked up by copy_functions
		end
	end

	local name = class.xarg.fullname
	for _, t in ipairs(temps) do
		-- extract name from [name, [disable function list] ]
		local disable_funcs = {}
		if type(t) == 'table' then
			t,disable_funcs = t[1],t[2]
		end
		local oclass, oparams = name:match('^(.+)<([^>]+)>$')
		local tclass, tparams = t:match('^(.+)<([^>]+)>$')
		if tclass == oclass then
			-- TODO: handle multiple template parameters
			local copy = deepcopy(class)
			template_repare(copy, oparams, tparams, disable_funcs)
			copy.xarg.cname = copy.xarg.fullname
				:gsub('[<>*]', '_')
				:gsub('::', '_LQT_')
				:gsub(',%s*', '_')
				:gsub('const ', '')
			ret[copy] = true
		else
			ignore(name, 'template not bound')
		end
	end
end

--- Append any created classes to the index
function finish(index)
	for f in pairs(idindex_add) do
		index[f] = true
	end
	idindex_add = {}
end
