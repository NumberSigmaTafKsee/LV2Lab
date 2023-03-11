/* include libs */
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <array>
#include <string>

#include "lv2.h"
#include <lv2/atom/atom.h>
#include <lv2/urid/urid.h>
#include <lv2/midi/midi.h>
#include <lv2/core/lv2_util.h>
#include <lv2/atom/util.h>

#include "LuaJIT.hpp"


std::string bundlepath;

struct Urids
{
    LV2_URID midi_MidiEvent;
};

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
	Lua::LuaJIT interpreter;
	
    LV2Lua(const std::string& path, const double sampleRate, const LV2_Feature *const *features) :    
    midi_in_ptr (nullptr),
    midi_out_ptr (nullptr),
    audio_in_ptr (nullptr),    
    audio_out_ptr (nullptr),    
    map (nullptr),
    rate (sampleRate),
    interpreter((path +"/lv2.lua").c_str())
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
	std::cout << "LV2Luajit Instantiated" << std::endl;
    LV2Lua *m = new LV2Lua(bundle_path,sample_rate, features);
    bundlepath = bundle_path;
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
    LV2Lua* m = (LV2Lua*) instance;
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
		m->midi_out_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
		break;
		
    case PORT_MIDI_IN:
        m->midi_in_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
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
    Lua::LuaJIT &interpreter = m->interpreter;
    uint32_t last_frame = 0;
    /* analyze incomming MIDI data */    
    LV2_ATOM_SEQUENCE_FOREACH (m->midi_in_ptr, ev)
    {
        /* play frames until event */
        const uint32_t frame = ev->time.frames;        
		last_frame = frame;
		
        if (ev->body.type == m->urids.midi_MidiEvent)
        {
            const uint8_t* const msg = reinterpret_cast<const uint8_t*> (ev + 1);
            const uint8_t typ = lv2_midi_message_type (msg);
			printf("MIDI Msg %d %d\n",msg[0],msg[1]);
            switch (typ)
            {
            case LV2_MIDI_MSG_NOTE_ON:            				
				interpreter.GetGlobal("note_on");
				interpreter.PushNumber(0);
                interpreter.PushNumber(msg[1]);
				interpreter.PushNumber(msg[2]);				                
                interpreter.pcall(3,0);                                
                break;

            case LV2_MIDI_MSG_NOTE_OFF:   
				interpreter.GetGlobal("note_off");
				interpreter.PushNumber(0);
                interpreter.PushNumber(msg[1]);
				interpreter.PushNumber(msg[2]);				                
                interpreter.pcall(3,0);                                
                break;

            case LV2_MIDI_MSG_CONTROLLER:
                interpreter.GetGlobal("control");
                interpreter.PushNumber(0);
                interpreter.PushNumber(msg[1]);
				interpreter.PushNumber(msg[2]);				                
                interpreter.pcall(3,0);                                
                break;
            
            default:
                break;
            }
        }
    }

	// run the audio code
	interpreter.GetGlobal("noise");
    interpreter.PushLightUserData(m->audio_in_ptr);
    interpreter.PushLightUserData(m->audio_out_ptr);
    interpreter.PushNumber(sample_count);
    interpreter.pcall(3,0);    
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
LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor (uint32_t index)
{
    if (index == 0) return &descriptor;
    else return NULL;
}
