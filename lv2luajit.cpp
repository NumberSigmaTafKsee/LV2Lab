/* include libs */
#include "stddef.h"
#include "stdint.h"
#include "stdlib.h"
//#include "math.h"
#include "lv2.h"

#include "LuaJIT.hpp"

Lua::LuaJIT interpreter("lv2.lua");


template <class T>
class LinearFader
{
private:
    T destination_;
    uint32_t distance_;
    T value_;

public:
    LinearFader (const T destination) :
        destination_ (destination),
        distance_ (0),
        value_ (destination)
    {

    }

    void set (const T destination, const uint32_t distance)
    {
        destination_ = destination;
        distance_ = distance;
        if (distance == 0) value_ = destination;
    }

    T get () const
    {
        return value_;
    }

    void proceed ()
    {
        if (distance_ == 0) value_ = destination_;
        else
        {
            value_ += (destination_ - value_) * (1.0 / static_cast<double> (distance_));
            distance_ -= 1;
        }
    }
};

template <class T>
T limit (const T x, const T min, const T max)
{
    return (x < min ? min : (x > max ? max : x));
}


/* class definition */
struct LV2Lua
{
    float* audio_in_ptr;
    float* audio_out_ptr;
    const LV2_Atom_Sequence* midi_in_ptr;
    const LV2_Atom_Sequence* midi_out_ptr;    
    LV2_URID_Map* map ;
    Urids urids;
    double rate;
    
    LV2Lua(const double sampleRate, const LV2_Feature *const features) :    
    midi_in_ptr (nullptr),
    audio_out_ptr (nullptr),
    controlLevel (0.0f),
    map (nullptr),
    rate (sample_rate),
    {
		const char* missing = lv2_features_query
		(
			features,
			LV2_URID__map, &map, true,
			NULL
		);

		if (missing) throw std::invalid_argument ("Feature map not provided by the host. Can't instantiate LV2Lua");

		urids.midi_MidiEvent = map->map (map->handle, LV2_MIDI__MidiEvent);
    }
} ;



/* internal core methods */
static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    
    LV2Lua *m = new LV2Lua(sample_rate, features);
    return m;
}

enum PortGroups
{
    PORT_MIDI_IN    = 0,
    PORT_MIDI_OUT   = 1,
    PORT_AUDIO_IN   = 2,
    PORT_AUDIO_OUT  = 3,
    PORT_CONTROL    = 4,
    PORT_NR         = 5
};

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
    MyAmp* m = (MyAmp*) instance;
    if (!m) return;

    switch (port)
    {
    case PORT_AUDIO_IN:
        m->audio_in_ptr = (float*) data_location;
        break;

    case PORT_AUDIO_OUT:
        m->audio_out_ptr = (float*) data_location;
        break;

    case PORT_MIDI_OUT:
		midi_out_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
		break;
		
    case PORT_MIDI_IN:
        midi_in_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
        break;

	case PORT_CONTROL:
		break;
		
	case PORT_NR:
		break;
		
    default:
        break;
    }
}

static void activate (LV2_Handle instance)
{
    /* not needed here */
    
}

static void run (LV2_Handle instance, uint32_t sample_count)
{
    LV2Lua* m = (LV2Lua*) instance;
    if (!m) return;
    
    /* analyze incomming MIDI data */
    uint32_t last_frame = 0;
    LV2_ATOM_SEQUENCE_FOREACH (m->midi_in_ptr, ev)
    {
        /* play frames until event */
        const uint32_t frame = ev->time.frames;
        play (last_frame, frame);
        last_frame = frame;

        if (ev->body.type == urids.midi_MidiEvent)
        {
            const uint8_t* const msg = reinterpret_cast<const uint8_t*> (ev + 1);
            const uint8_t typ = lv2_midi_message_type (msg);

            switch (typ)
            {
            case LV2_MIDI_MSG_NOTE_ON:            
                interpreter.PushNumber(msg[0]);
                interpreter.PushNumber(msg[1]);
                interpreter.PushNumber(msg[2]);
                interpreter.Call("note_on");                
                break;

            case LV2_MIDI_MSG_NOTE_OFF:                
                interpreter.PushNumber(msg[1]);
                interpreter.PushNumber(msg[2]);
                interpreter.Call("note_off");                
                break;

            case LV2_MIDI_MSG_CONTROLLER:
                interpreter.PushNumber(msg[0]);
                interpreter.PushNumber(msg[1]);
                interpreter.PushNumber(msg[2]);
                interpreter.Call("control");                
                break;
            
            default:
                break;
            }
        }
    }

	// run the audio code
    interpreter.PushLightUserData(m->audio_in_ptr);
    interpreter.PushLightUserData(m->audio_out_ptr);
    interpreter.PushLightUserData(sample_count);
    interpreter.Call("noise");
    
}

static void deactivate (LV2_Handle instance)
{
    /* not needed here */
}

static void cleanup (LV2_Handle instance)
{
    LV2Lua * l = (LV2Lua*)instance;
    if(l) delete l;
}

static const void* extension_data (const char *uri)
{
    return NULL;
}

/* descriptor */
static LV2_Descriptor const descriptor =
{
    "urn:james5:lv2luajit",
    instantiate,
    connect_port,
    activate /* or NULL */,
    run,
    deactivate /* or NULL */,
    cleanup,
    extension_data /* or NULL */
};

/* interface */
const LV2_SYMBOL_EXPORT LV2_Descriptor* lv2_descriptor (uint32_t index)
{
    if (index == 0) return &descriptor;
    else return NULL;
}
