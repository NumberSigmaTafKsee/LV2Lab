-- use xoscope to see the waveform on the device
require('audiosystem')
require('sndfile')
require('simpleresampler')

s = sndfile.SndFileReaderFloat("AcGtr.wav")
w = audiosystem.FloatWave()
v = sndfile.float_vector(s:size())
s:read(v:size(),v:data())
audiosystem.InitWave(w,v:data(),s:frames(),s:channels(),44100,true)
audiosystem.ReverseWave(w)

function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

function noise(input,output,frames)        
    local buffer = new_buffer(output)    
    for i=0,frames-1 do buffer[i]=0 end
end 

num_keys=0
function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end
function note_on(c,n,v)      
    local f = math.pow(2.0, (n-69)/12)*440.0            
    freq = f    
    num_keys = num_keys+1
end
function note_off(c,n,v) 
    num_keys = num_keys-1
    if(num_keys <= 0) then
        num_keys = 0;
    end       
end
function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 102) then
        filt:setCutoff(22050.0*d2/127)
    elseif(d1 == 103) then
        q = 10*d2/127
        if(q >= 0.5) then
            filt:setQ(q)
        end
    end
end
-- app callback, midi handling and logic
-- isAudioRunning shuld be atomic
-- either block audio or wait until finished
-- this is run every 10ms, or can be changed in portaudio.i
function callback()
    print('hi')
end 

audiosystem.Init()
audiosystem.Pm_Initialize()

audiosystem.set_note_on_func(note_on)
audiosystem.set_note_off_func(note_off)
audiosystem.set_control_change_func(control)

for i=0,audiosystem.GetNumMidiDevices()-1 do 
    print(i,audiosystem.GetMidiDeviceName(i))
end

audiosystem.set_audio_func(noise)
--audiosystem.set_callback_func(callback)
device=14
audiosystem.Pa_Initialize()
for i=0,audiosystem.GetNumAudioDevices()-1 do 
    if( audiosystem.GetAudioDeviceName(i) == 'jack') then        
        device = i 
        goto done
    end    
end
::done::

audiosystem.InitMidiDevice(1,3,3)
audiosystem.InitAudioDevice(device,-1,1,44100,256)
audiosystem.PlayWave(w)
audiosystem.RunAudio()
audiosystem.StopAudio()
