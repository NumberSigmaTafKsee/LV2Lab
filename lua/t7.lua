#!/usr/bin/lua
dofile(arg[0]:gsub('test[/\\].+', 'examples/init.lua'))

local QtCore = require 'qtcore'
local QtGui = require 'qtgui'
local QtWidgets = require 'qtwidgets'

local LCD_Range = function(...)
	local this = QtWidgets.QWidget.new(...)
	--print(this:metaObject():className(), this:metaObject():methodCount())
	--print(this:metaObject():className(), this:metaObject():methodCount())

	local lcd = QtWidgets.QLCDNumber.new()
	lcd:setSegmentStyle 'Filled'

	local slider = QtWidgets.QSlider.new'Horizontal'
	slider:setRange(0, 99)
	slider:setValue(0)

	this:__addsignal('valueChanged(int)')
	this:__addslot('setValue(int)', function(_, val) slider:setValue(val) end)

	-- slider:connect('2valueChanged(int)', lcd, '1display(int)')
	QtCore.QObject.connect(slider, '2valueChanged(int)', lcd, '1display(int)')
	QtCore.QObject.connect(slider, '2valueChanged(int)', this, '2valueChanged(int)')

	local layout = QtWidgets.QVBoxLayout.new()
	layout:addWidget(lcd)
	layout:addWidget(slider)
	this:setLayout(layout)
	return this
end

local new_MyWidget = function(...)
	local this = QtWidgets.QWidget.new(...)

	local quit = QtWidgets.QPushButton.new('Quit')
	quit:setFont(QtGui.QFont('Times', 18, 75))

	QtCore.QObject.connect(quit, '2clicked()', this, '1close()')

	local grid = QtWidgets.QGridLayout.new()
	local previousRange = nil
	for row = 1, 3 do
		for column = 1, 3 do
			local lcdrange = LCD_Range()
			grid:addWidget(lcdrange, row, column)
			if previousRange then
				QtCore.QObject.connect(lcdrange, '2valueChanged(int)',
					previousRange, '1setValue(int)')
			end
			previousRange = lcdrange
		end
	end

	local layout = QtWidgets.QVBoxLayout.new()
	layout:addWidget(quit)
	layout:addLayout(grid)
	this:setLayout(layout)
	return this
end

local app = QtWidgets.QApplication.new(1 + select('#', ...), {arg[0], ...})
app.__gc = app.delete -- take ownership of object

local widget = new_MyWidget()
widget:show()
print('gc begin')
collectgarbage()
print('gc end')

app.exec()


