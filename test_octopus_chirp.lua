--require('octopus')
--i = octopus.OctaveInterpreter()
v = octopus.OctopusValueList()

v:set(0,0.0)
v:set(1,1/8000)
v:set(2,1024.0)
r = octopus.octave_linspace(v,1)

r:set(1,200)
r:set(2,2)
r:set(3,500)
r:set(4,"logarithmic")
v = octopus.octave_chirp(r,1)

Plot(v)
--octopus.octave_plot(v)
--octopus.octave_pause(v)
--v:set(1,256)
--v:set(2,8000)
--octopus.specgram(v)
--t = v:get(0):getRowFloatVector()
--x = octopus.Octavate()
--r = v:get(0):getFloatRowVector()
--p = plot.Plot_Float()
--p:setstyle("lines")
--p:plot_x(r:data(),r:size(),"foo")
os.execute("sleep 10;")