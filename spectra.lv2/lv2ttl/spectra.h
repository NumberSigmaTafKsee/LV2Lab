// generated by lv2ttl2c from
// http://gareus.org/oss/lv2/spectra#Mono

extern const LV2_Descriptor* lv2_descriptor(uint32_t index);
extern const LV2UI_Descriptor* lv2ui_descriptor(uint32_t index);

static const RtkLv2Description _plugin = {
	&lv2_descriptor,
	&lv2ui_descriptor
	, 0 // uint32_t dsp_descriptor_id
	, 0 // uint32_t gui_descriptor_id
	, "Spectr" // const char *plugin_human_id
	, (const struct LV2Port[7])
	{
		{ "control", ATOM_IN, nan, nan, nan, "GUI to plugin communication"},
		{ "notify", ATOM_OUT, nan, nan, nan, "Plugin to GUI communication"},
		{ "fftsize", CONTROL_IN, 4096.000000, 1024.000000, 16384.000000, "FFT Size"},
		{ "color", CONTROL_IN, 0.000000, 0.000000, 1.000000, "1/f scale"},
		{ "window", CONTROL_IN, 0.000000, 0.000000, 4.000000, "Window Function"},
		{ "in", AUDIO_IN, nan, nan, nan, "Audio Input"},
		{ "out", AUDIO_OUT, nan, nan, nan, "Audio Signal pass-thru"},
	}
	, 7 // uint32_t nports_total
	, 1 // uint32_t nports_audio_in
	, 1 // uint32_t nports_audio_out
	, 0 // uint32_t nports_midi_in
	, 0 // uint32_t nports_midi_out
	, 1 // uint32_t nports_atom_in
	, 1 // uint32_t nports_atom_out
	, 3 // uint32_t nports_ctrl
	, 3 // uint32_t nports_ctrl_in
	, 0 // uint32_t nports_ctrl_out
	, 33024 // uint32_t min_atom_bufsiz
	, false // bool send_time_info
	, UINT32_MAX // uint32_t latency_ctrl_port
};