#!/usr/bin/lua
dofile(arg[0]:gsub('test[/\\].+', 'examples/init.lua'))

local QtCore = require 'qtcore'

local qa = QtCore.QCoreApplication.new(1, {'virt_test'})

qa:__addsignal('valueChanged(bool,double,QString,QObject*,int)')
qa:__addslot('setValue(bool,double,QString,QObject*,int)', function(self, arg1, arg2, arg3, arg4, arg5)
	print('setValue:', self, arg1, arg2, arg3:toStdString(), arg4, arg5)
end)

QtCore.QObject.connect(qa, '2valueChanged(bool,double,QString,QObject*,int)'
	, qa, '1setValue(bool,double,QString,QObject*,int)'
)

qa:connect('2valueChanged(bool,double,QString,QObject*,int)'
	, qa, '1setValue(bool,double,QString,QObject*,int)'
)

qa:connect('2destroyed(QObject*)', function()
end)
qa:connect('2applicationNameChanged()', function()
    print('applicationNameChanged: ', qa.applicationName():toStdString())
end)
table.foreach(qa:__methods(), print)

-- lqt-specific fields
local LQT_OBJMETASTRING = "Lqt MetaStringData"
local LQT_OBJMETADATA = "Lqt MetaData"
local LQT_OBJSLOTS = "Lqt Slots"
local LQT_OBJSIGS = "Lqt Signatures"

print(qa[LQT_OBJMETADATA])
print(qa[LQT_OBJMETASTRING])
print(qa['*' .. LQT_OBJMETADATA])
print(qa['*' .. LQT_OBJMETASTRING])

qa.setApplicationName('New Application Name')

local meta = qa:metaObject()
local signalIndex = meta:indexOfMethod('setValue(bool,double,QString,QObject*,int)')
local signal = meta:method(signalIndex)
print('signal:', meta, signalIndex, signal)

signal:invoke(qa, 'AutoConnection', 3.14, '7758521', qa, 7)

QtCore.QMetaObject.invokeMethod(qa
	, 'valueChanged'
	-- , 'setValue'
	, 'AutoConnection'
	, true
	, 3.14, '9.28'
	, { 'QObject*', qa }
	, { 'int', 7 }
)

qa:__emit('setValue', true, 3.1415926, '9.28wtf', { 'QObject*', qa }, { 'int', 7 })
-- qa:__emit('setValue', false, 3.1415926, '9.28wtf', { 'QObject*', nil }, { 'int', 7 })

-- qa:__emit('valueChanged', 1234, 5678)
