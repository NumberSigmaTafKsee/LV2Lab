require('stk')
require('AudioTK')
require('Amplifiers')
require('audio_lfo')

lfo = audio_lfo.LFO(44100)
lfo:setRate(0.00005)

stk.Stk.setRawwavePath("Data/rawwaves")
stk.Stk.setSampleRate(44100)

chorus = stk.Chorus()
chorus:setModFrequency(2)

diode = AudioTK.MonoHalfTanhShaper()
delay = AudioTK.MonoFDNDelay(44100)

-- it needs to be oversampled to reduce the aliasing
clip  = Amplifiers.MorphClipper(Amplifiers.ClipFunction.SERPENT_CURVE,Amplifiers.ClipFunction.ERFMOIDER)

diode:setPort(AudioTK.MonoHalfTanhShaper.PORT_COEFF,10)

delay:setPort(AudioTK.MonoFDNDelay.PORT_INGAIN,1)
delay:setPort(AudioTK.MonoFDNDelay.PORT_OUTGAIN,1)
delay:setPort(AudioTK.MonoFDNDelay.PORT_FEEDBACK,0.9)
delay:setPort(AudioTK.MonoFDNDelay.PORT_DELAY,0.5*44100)

-- bowed doesn't work either ...
bow = stk.Bowed()
bow:clear()



fc = 440.0
q  = 0.5

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end

function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

v = AudioTK.double_vector(256)
v2 = AudioTK.double_vector(256)

function noise(input,output,frames)            
    local outbuf = new_buffer(output)         
    for i=0,frames-1 do        
        local t = lfo:tick()        
        outbuf[i] = botl:tick()
    end
    
    for i=0,frames-1 do
        v[i+1] = outbuf[i]
    end
    diode:ProcessBlock(frames,v:data(),v:data())
    delay:ProcessBlock(frames,v:data(),v2:data())
    for i=0,frames-1 do
        outbuf[i] = v2[i+1] + v[i+1]
    end
    
end 

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end

function note_on(c,n,v)        
    local f = midi_to_freq(n)
    botl:noteOn(f,1)
end

function note_off(c,n,v)          
    botl:noteOff(1)
end

function pitchbend(c,d1,d2)
    local x = 127*d2
    
end

function control(c,d1,d2)
    
    if(d1 == 1) then
        print(c,d1,d2)
    elseif(d1 == 102) then        
        fc = d2/127                
    elseif(d1 == 103) then
        q = d2/127                
    end
end


function randomize()   
    --[[
    for i=1,8 do     
        b3 = voices[i]
        b3:setRatio(0,math.random()*4);
        b3:setRatio(1,math.random()*4);
        b3:setRatio(2,math.random()*4);
        b3:setRatio(3,math.random()*4);
        b3:setGain(0,math.random()*4);
        b3:setGain(1,math.random()*4);
        b3:setGain(2,math.random()*4);
        b3:setGain(3,math.random()*4);
        b3:setModulationSpeed(math.random()*10);
        b3:setModulationDepth(math.random());
        --b3:setControl1(math.rand()*10);
        --b3:setControl2(math.rand()*10);
    end
    ]]
end
