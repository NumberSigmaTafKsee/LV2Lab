require('audiosystem')
require('lv2plugin')
lv2plugins = lv2plugin.LV2Plugins()
lv2luajit = lv2plugins:LoadPlugin("urn:james5:lv2luajit")
