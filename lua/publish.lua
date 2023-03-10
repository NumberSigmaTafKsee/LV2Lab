#!/usr/bin/env lua
-- amalgamate and render the markdown wiki pages as a single HTML file

-- some variables you might want to tweak:
local TITLE = "fltk4lua" -- title of the resulting HTML file
local MAIN = "Home" -- starting page of the wiki
local OUT = "index.html" -- output HTML file name


local md = require( "markdown" )

local function readall( fname, files )
  local file = assert( io.open( fname, "r" ) )
  local contents = file:read( "*a" )
  file:close()
  return contents:gsub( "%[([^%]%[]-)%]%(([%w%._]-)%)", function( text, link )
    if not files[ link ] then
      files[ link ] = true
      files[ #files+1 ] = link
    end
    return "["..text.."](#"..link..")"
  end )
end

local function make_html_name( s )
  return (s:gsub( "[^%w._]", "_" ))
end

local filelist = { MAIN, [ MAIN ] = true }
local contents = {}
for i,fname in ipairs( filelist ) do
  if i == 1 then
    contents[ i ] = md( readall( fname..".md", filelist ) ) ..
                    "<hr />\n\n"
  else
    local name = make_html_name( fname )
    contents[ i ] = '<a name="'..name..'" id="'..name..'"></a>'
                    .. md( readall( fname..".md", filelist ) )
  end
end

local header = [[
<!DOCTYPE html>
<html>
<head>
  <meta http-equiv="content-type" content="text/html; charset=UTF-8" />
  <title>]]..TITLE..[[</title>
  <style type="text/css">
body {
  margin: 30px;
  padding: 0px;
  font-size: 17px;
  font-family: helvetica,geneva,sans-serif;
  max-width: 50em;
}
p, pre, blockquote {
  margin-left: 20px;
  margin-right: 20px;
}
h1 {
  text-align: center;
  font-size: 36px;
  font-weight: bold;
}
h2, h3, h4, h5, h6 {
  margin-top: 35px;
  padding: 0px 40px 0px 40px;
  background-color: #E7E7E7;
  line-height: 30px;
  font-weight: normal;
}
h2 {
  font-size: 26px;
  font-weight: bold;
}
h3 {
  font-size: 22px;
}
h4 {
  font-size: 20px;
  font-weight: bold;
}
:target + h2, :target + h3 {
  padding: 2px 40px 2px 40px;
  border: solid #A0A0A0 2px;
  border-radius: 8px;
}
a {
  color: #000088;
  text-decoration: none;
  font-weight: bold;
}
a:hover {
  color: #FF0000;
  background-color: #E7E7E7;
  text-decoration: underline;
}
code {
  font-weight: bold;
}
pre {
  background-color: #F0F0F0;
  border-left: 1px black solid;
  padding: 5px 15px 5px 15px;
}
pre > code {
  font-weight: normal;
}
  </style>
</head>
<body>
]]
local footer = [[
</body>
</html>
]]

local out = assert( io.open( OUT, "w" ) )
out:write( header )
for _,c in ipairs( contents ) do
  out:write( c )
end
out:write( footer )
out:close()
print( "done!" )

