require('stk')
require('AudioTK')
require('Amplifiers')
require('floatbuffer')

stk.setRawwavePath("Data/rawwaves")
stk.setSampleRate(44100)
adsr = stk.ADSR()
chorus = stk.Chorus()
chorus:setModFrequency(12)
diode = AudioTK.MonoHalfTanhShaper()
delay = AudioTK.MonoFDNDelay(44100)
-- it needs to be oversampled to reduce the aliasing
clip  = Amplifiers.MorphClipper(Amplifiers.ClipFunction.SERPENT_CURVE,Amplifiers.ClipFunction.ERFMOIDER)
diode:setPort(AudioTK.MonoHalfTanhShaper.PORT_COEFF,10)
delay:setPort(AudioTK.MonoFDNDelay.PORT_INGAIN,1)
delay:setPort(AudioTK.MonoFDNDelay.PORT_OUTGAIN,1)
delay:setPort(AudioTK.MonoFDNDelay.PORT_FEEDBACK,0.9)
delay:setPort(AudioTK.MonoFDNDelay.PORT_DELAY,0.5*44100)
adsr:setAllTimes(0.1,0.3,0.8,0.2)
b3   = stk.BeeThree()
envs   = {}
voices = {}
inuse  = {}
freqs  = {}
vels   = {}
voicer = stk.Voicer()

for i=1,8 do
    envs[i] = stk.ADSR()
    envs[i]:setAllTimes(0.1,0.3,0.8,0.2)
    voices[i] = stk.BeeThree()    
    voices[i]:noteOff(0)
    inuse[i] = false
    freqs[i] = 0
    vels[i] = 0
    voicer:addInstrument(voices[i])
end


fc = 440.0
q  = 0.5
num_pressed = 0

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end


function noise(input,output,frames)            
	outbuf = floatbuffer.FloatBuffer(frames,output)	
	for i=0,frames-1 do                
		local sample = voicer:tick()
		outbuf[i] = sample		
	end
end 
function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end

function note_on(c,n,v)    	
    voicer:noteOn(n,v)
    num_pressed=num_pressed+1
end

function note_off(c,n,v)        
    voicer:noteOff(n,v)
    num_pressed = num_pressed-1
end

function pitchbend(c,d1,d2)
    local x = 127*d2
    voicer:pitchBend(x)
end

function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 1) then
        voicer:controlChange(11,d2)
    elseif(d1 == 102) then        
        fc = d2/127                
    elseif(d1 == 103) then
        q = d2/127                
    end
end


function randomize()   
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
end
