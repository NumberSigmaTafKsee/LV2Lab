# SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
# SPDX-License-Identifier: CC0-1.0

project('mephisto.lv2', 'c', default_options : [
	'buildtype=release',
	'warning_level=3',
	'werror=false',
	'b_lto=true',
	'c_std=gnu11'])

d2tk = subproject('d2tk')
varchunk = subproject('varchunk')
ser_atom_lv2 = subproject('ser_atom.lv2')
timely_lv2 = subproject('timely.lv2')
props_lv2 = subproject('props.lv2')

lv2libdir = get_option('lv2libdir')
build_tests = get_option('build-tests')
use_fontconfig = get_option('use-fontconfig')

inst_dir = join_paths(lv2libdir, meson.project_name())

if get_option('use-backend-nanovg').enabled()
	d2tk_dep = d2tk.get_variable('d2tk_nanovg')
elif get_option('use-backend-cairo').enabled()
	d2tk_dep = d2tk.get_variable('d2tk_cairo')
else
	error('no valid UI backend given')
endif

varchunk_dep = varchunk.get_variable('varchunk')
ser_atom_lv2_dep = ser_atom_lv2.get_variable('ser_atom_lv2')
timely_lv2_dep = timely_lv2.get_variable('timely_lv2')
props_lv2_dep = props_lv2.get_variable('props_lv2')

source_root = meson.source_root()
build_root = meson.build_root()

add_project_arguments('-D_GNU_SOURCE', language : 'c')

conf_data = configuration_data()
cc = meson.get_compiler('c')

lv2_validate = find_program('lv2_validate', native : true, required : false)
sord_validate = find_program('sord_validate', native : true, required : false)
lv2lint = find_program('lv2lint', required : false)
reuse = find_program('reuse', required : false)
faust = find_program('faust', required : false)

if cc.has_member('UIGlue', 'addSoundFile',
		prefix : '#include <faust/gui/CInterface.h>')
	add_project_arguments('-DFAUST_ADDSOUNDFILE=addSoundFile', language : 'c')
else
	add_project_arguments('-DFAUST_ADDSOUNDFILE=addSoundfile', language : 'c')
endif

if cc.has_header(join_paths('faust', 'dsp', 'llvm-dsp-c.h'))
	add_project_arguments('-DFAUST_LLVM_DSP_C_H=<faust/dsp/llvm-dsp-c.h>', language : 'c')
else
	add_project_arguments('-DFAUST_LLVM_DSP_C_H=<faust/dsp/llvm-c-dsp.h>', language : 'c')
endif

m_dep = cc.find_library('m')
lv2_dep = dependency('lv2', version : '>=1.16.0')
faust_dep = cc.find_library('faust')

if cc.has_member('LV2UI_Request_Value', 'request',
		prefix : '#include <lv2/lv2plug.in/ns/extensions/ui/ui.h>')
	add_project_arguments('-D_LV2_HAS_REQUEST_VALUE', language : 'c')
	message('building with ui:requestValue support')
endif

dsp_deps = [m_dep, lv2_dep, faust_dep, varchunk_dep, timely_lv2_dep, props_lv2_dep]
ui_deps = [lv2_dep, d2tk_dep, ser_atom_lv2_dep, props_lv2_dep]

inc_dir = []

dsp_srcs = [join_paths('src', 'mephisto.c')]

ui_srcs = [join_paths('src', 'mephisto_ui.c')]

c_args = ['-fvisibility=hidden']

version = get_option('version').split('.')
add_project_arguments('-DMEPHISTO_VERSION="'
	+ version[0] + '.' + version[1] + '.' + version[2] + '"', language : 'c')
conf_data.set('MAJOR_VERSION', version[0])
conf_data.set('MINOR_VERSION', version[1])
conf_data.set('MICRO_VERSION', version[2])

mod = shared_module('mephisto', dsp_srcs,
	c_args : c_args,
	include_directories : inc_dir,
	name_prefix : '',
	dependencies : dsp_deps,
	install : true,
	install_dir : inst_dir)

ui = shared_module('mephisto_ui', ui_srcs,
	c_args : c_args,
	include_directories : inc_dir,
	name_prefix : '',
	dependencies : ui_deps,
	install : true,
	install_dir : inst_dir)

bank_filter_through_dsp = join_paths(
  meson.current_source_dir(), 'dsp', 'bank-filter_through.dsp')
bank_filter_gain_dsp = join_paths(
  meson.current_source_dir(), 'dsp', 'bank-filter_gain.dsp')
bank_time_lfo_dsp = join_paths(
  meson.current_source_dir(), 'dsp', 'bank-time_lfo.dsp')
bank_instrument_osc_dsp = join_paths(
  meson.current_source_dir(), 'dsp', 'bank-instrument_osc.dsp')
bank_analyzer_vu_meter_dsp = join_paths(
  meson.current_source_dir(), 'dsp', 'bank-analyzer_vu-meter.dsp')

bank_filter_through = run_command('cat', bank_filter_through_dsp,
  check : true)
conf_data.set('BANK-FILTER_THROUGH', bank_filter_through.stdout())

bank_filter_gain = run_command('cat', bank_filter_gain_dsp,
  check : true)
conf_data.set('BANK-FILTER_GAIN', bank_filter_gain.stdout())

bank_time_lfo = run_command('cat', bank_time_lfo_dsp,
  check : true)
conf_data.set('BANK-TIME_LFO', bank_time_lfo.stdout())

bank_instrument_osc = run_command('cat', bank_instrument_osc_dsp,
  check : true)
conf_data.set('BANK-INSTRUMENT_OSC', bank_instrument_osc.stdout())

bank_analyzer_vu_meter = run_command('cat', bank_analyzer_vu_meter_dsp,
  check : true)
conf_data.set('BANK-ANALYZER_VU-METER', bank_analyzer_vu_meter.stdout())

suffix = mod.full_path().strip().split('.')[-1]
conf_data.set('MODULE_SUFFIX', '.' + suffix)

manifest_ttl = configure_file(
	input : join_paths('ttl', 'manifest.ttl.in'),
	output : 'manifest.ttl',
	configuration : conf_data,
	install : true,
	install_dir : inst_dir)

dsp_ttl = configure_file(
	input : join_paths('ttl', 'mephisto.ttl.in'),
	output : 'mephisto.ttl',
	configuration : conf_data,
	install : true,
	install_dir : inst_dir)

ui_ttl = configure_file(
	input : join_paths('ttl', 'mephisto_ui.ttl'),
	output : 'mephisto_ui.ttl',
	copy: true,
	install : true,
	install_dir : inst_dir)

pset_ttl = configure_file(
	input : join_paths('ttl', 'presets.ttl.in'),
	output : 'presets.ttl',
	configuration : conf_data,
	install : true,
	install_dir : inst_dir)

alert_triangle_png = configure_file(
	input : join_paths('png', 'alert-triangle.png'),
	output : 'alert-triangle.png',
	copy : true,
	install : true,
	install_dir : inst_dir)
clipboard_png = configure_file(
	input : join_paths('png', 'clipboard.png'),
	output : 'clipboard.png',
	copy: true,
	install : true,
	install_dir : inst_dir)
copy_png = configure_file(
	input : join_paths('png', 'copy.png'),
	output : 'copy.png',
	copy: true,
	install : true,
	install_dir : inst_dir)
delete_png = configure_file(
	input : join_paths('png', 'delete.png'),
	output : 'delete.png',
	copy: true,
	install : true,
	install_dir : inst_dir)
eye_off_png = configure_file(
	input : join_paths('png', 'eye-off.png'),
	output : 'eye-off.png',
	copy: true,
	install : true,
	install_dir : inst_dir)
eye_png = configure_file(
	input : join_paths('png', 'eye.png'),
	output : 'eye.png',
	copy: true,
	install : true,
	install_dir : inst_dir)
save_png = configure_file(
	input : join_paths('png', 'save.png'),
	output : 'save.png',
	copy: true,
	install : true,
	install_dir : inst_dir)

if not use_fontconfig.enabled()
	fira_sans_bold_ttf = d2tk.get_variable('fira_sans_bold_ttf')
	fira_code_bold_ttf = d2tk.get_variable('fira_code_bold_ttf')
	fira_code_light_ttf = d2tk.get_variable('fira_code_light_ttf')
	fira_code_medium_ttf = d2tk.get_variable('fira_code_medium_ttf')
	fira_code_regular_ttf = d2tk.get_variable('fira_code_regular_ttf')

	configure_file(
		input : fira_sans_bold_ttf,
		output : 'FiraSans:bold.ttf',
		copy : true,
		install : true,
		install_dir : inst_dir)
	configure_file(
		input : fira_code_bold_ttf,
		output : 'FiraCode:bold.ttf',
		copy : true,
		install : true,
		install_dir : inst_dir)
	configure_file(
		input : fira_code_light_ttf,
		output : 'FiraCode:light.ttf',
		copy : true,
		install : true,
		install_dir : inst_dir)
	configure_file(
		input : fira_code_medium_ttf,
		output : 'FiraCode:medium.ttf',
		copy : true,
		install : true,
		install_dir : inst_dir)
	configure_file(
		input : fira_code_regular_ttf,
		output : 'FiraCode:regular.ttf',
		copy : true,
		install : true,
		install_dir : inst_dir)
endif

if build_tests
	if lv2_validate.found() and sord_validate.found()
		test('LV2 validate', lv2_validate,
			args : [manifest_ttl, dsp_ttl, ui_ttl, pset_ttl])
	endif

	if lv2lint.found()
		test('LV2 lint', lv2lint,
			args : ['-M', 'pack',
				'-E', 'warn',
				'-I', join_paths(build_root, ''),
				'http://open-music-kontrollers.ch/lv2/mephisto#audio_1x1',
				'http://open-music-kontrollers.ch/lv2/mephisto#audio_2x2',
				'http://open-music-kontrollers.ch/lv2/mephisto#audio_4x4',
				'http://open-music-kontrollers.ch/lv2/mephisto#audio_8x8',
				'http://open-music-kontrollers.ch/lv2/mephisto#cv_1x1',
				'http://open-music-kontrollers.ch/lv2/mephisto#cv_2x2',
				'http://open-music-kontrollers.ch/lv2/mephisto#cv_4x4',
				'http://open-music-kontrollers.ch/lv2/mephisto#cv_8x8'
			])
	endif

  if reuse.found() and not meson.is_cross_build()
    test('REUSE', reuse, args : [
      '--root', source_root,
      'lint'
    ])
  endif

	if faust.found() and not meson.is_cross_build()
		test('FAUST bank-filter_through', faust,
      args : [ bank_filter_through_dsp ])
		test('FAUST bank-filter_gain', faust,
      args : [ bank_filter_gain_dsp ])
		test('FAUST bank-time_lfo', faust,
      args : [ bank_time_lfo_dsp ])
		test('FAUST bank-instrument_osc', faust,
      args : [ bank_instrument_osc_dsp ])
		test('FAUST bank-analyzer_vu-meter', faust,
      args : [ bank_analyzer_vu_meter_dsp ])
	endif
endif
