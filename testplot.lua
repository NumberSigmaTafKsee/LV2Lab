require('stdsamples')
require('plot')
x = stdsamples.generate_sin(440.0,44100.0,1024)
p = plot.Plot_Float()
--[[
p:clear()
p:reset()
p:key(false)
p:border(3)
p:xrange(0,20)
p:yrange(0,11)
p:unset_yticks()
p:title("Sin")
p:set_xlabel("x")
p:set_ylabel("y")
--]]
p:setstyle("lines")
p:plot_x(x:data(),x:size(),"Sin")
os.execute("sleep 60")