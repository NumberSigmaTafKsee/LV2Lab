# SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
# SPDX-License-Identifier: CC0-1.0

project('d2tk', 'c', default_options : [
	'buildtype=release',
	'warning_level=3',
	'werror=false',
	'b_lto=true',
	'c_std=gnu11'])

static_link = false #meson.is_cross_build()

build_debug_overlay = get_option('build-debug-overlay')
build_examples = get_option('build-examples')
build_tests = get_option('build-tests')
build_doc = get_option('build-doc')

use_backend_cairo = get_option('use-backend-cairo')
use_backend_nanovg = get_option('use-backend-nanovg')
use_frontend_fbdev = get_option('use-frontend-fbdev')
use_frontend_pugl = get_option('use-frontend-pugl')
use_frontend_glfw = get_option('use-frontend-glfw')

use_evdev = get_option('use-evdev')
use_fontconfig = get_option('use-fontconfig')

grep = find_program('grep',
	native : true,
	required : use_fontconfig)
fc_list = find_program('fc-list',
	native : true,
	required : use_fontconfig)
check_for_font = find_program('check_for_font',
	native : true,
	required : use_fontconfig)
robodoc = find_program('robodoc',
	native : true,
	required : build_doc)
reuse = find_program('reuse',
  required : false)

cc = meson.get_compiler('c')

# mandatory dependencies
m_dep = cc.find_library('m')
thread_dep = dependency('threads')

# dependencies for backend_cairo/frontend_pugl
freetype_dep = dependency('freetype2',
	version : '>=18.0.0',
	static : static_link,
	required : use_backend_cairo)
pixman_dep = dependency('pixman-1',
	version : '>=0.34.0',
	static : static_link,
	required : use_backend_cairo)
cairo_dep = dependency('cairo',
	version : '>=1.14.0',
	static : static_link,
	required : use_backend_cairo)
cairo_deps = [freetype_dep, pixman_dep, cairo_dep]
if use_frontend_pugl.enabled()
	cairo_xlib_dep = dependency('cairo-xlib',
		version : '>=1.14.0',
		static : static_link,
		required : use_backend_cairo)
endif

# dependencies for frontend_fbdev
input_dep = dependency('libinput',
	version : '>=1.6.0',
	static : static_link,
	required : use_frontend_fbdev)
udev_dep = dependency('libudev',
	static : static_link,
	required : use_frontend_fbdev)
evdev_dep = dependency('libevdev',
	version : '>=1.5.0',
	static : static_link,
	required : use_frontend_fbdev)

# dependencies for backend_nanovg
glu_dep = dependency('glu',
	version : '>=9.0.0',
	static : static_link,
	required : use_backend_nanovg)
glew_dep = dependency('glew',
	version : '>=2.0.0',
	static : static_link,
	required : use_backend_nanovg)
if use_frontend_glfw.enabled()
	glfw_dep = dependency('glfw3',
		version : '>=3.3.0',
		static : static_link,
		required : use_backend_nanovg)
endif

# optional dependencies
util_dep = cc.find_library('util')
vterm_dep = dependency('vterm',
	version : '>=0.1',
	static : static_link)
if not use_frontend_fbdev.enabled()
	evdev_dep = dependency('libevdev',
		version : '>=1.5.0',
		static : static_link,
		required : use_evdev)
endif
fontconfig_dep = dependency('fontconfig',
	version : '>=2.0.0',
	static : static_link,
	required : use_fontconfig)

deps = [m_dep, thread_dep, evdev_dep, util_dep, vterm_dep, fontconfig_dep]
links = []

d2tk_inc = include_directories('')
pugl_inc = include_directories(join_paths('pugl', 'include'))
nanovg_inc = include_directories(join_paths('nanovg', 'src'))
linenoise_inc = include_directories('linenoise')
inc_dir = [d2tk_inc, pugl_inc, nanovg_inc, linenoise_inc]

version = get_option('version').split('.')

conf_data = configuration_data()
conf_data.set('MAJOR_VERSION', version[0])
conf_data.set('MINOR_VERSION', version[1])
conf_data.set('MICRO_VERSION', version[2])

add_project_arguments('-D_GNU_SOURCE', language : 'c')

if build_debug_overlay
	conf_data.set('D2TK_DEBUG', 1)
else
	conf_data.set('D2TK_DEBUG', 0)
endif

lib_srcs = [
	join_paths('src', 'hash.c'),
	join_paths('src', 'core.c'),
	join_paths('src', 'base.c'),
	join_paths('src', 'base_table.c'),
	join_paths('src', 'base_frame.c'),
	join_paths('src', 'base_layout.c'),
	join_paths('src', 'base_scrollbar.c'),
	join_paths('src', 'base_pane.c'),
	join_paths('src', 'base_cursor.c'),
	join_paths('src', 'base_button.c'),
	join_paths('src', 'base_image.c'),
	join_paths('src', 'base_bitmap.c'),
	join_paths('src', 'base_custom.c'),
	join_paths('src', 'base_meter.c'),
	join_paths('src', 'base_combo.c'),
	join_paths('src', 'base_textfield.c'),
	join_paths('src', 'base_label.c'),
	join_paths('src', 'base_separator.c'),
	join_paths('src', 'base_tooltip.c'),
	join_paths('src', 'base_link.c'),
	join_paths('src', 'base_dial.c'),
	join_paths('src', 'base_spinner.c'),
	join_paths('src', 'base_bar.c'),
	join_paths('src', 'base_wave.c'),
	join_paths('src', 'base_flowmatrix.c'),
	join_paths('src', 'base_pty.c'),
	join_paths('src', 'base_lineedit.c'),
	join_paths('src', 'util_spawn.c'),
	join_paths('linenoise', 'linenoise.c'),
	join_paths('linenoise', 'encodings', 'utf8.c')
]

if use_evdev.enabled()
	conf_data.set('D2TK_EVDEV', 1)
	lib_srcs += join_paths('src', 'base_vkb.c')
else
	conf_data.set('D2TK_EVDEV', 0)
endif

if input_dep.found() and input_dep.version().version_compare('>=1.15.0')
	conf_data.set('D2TK_INPUT_1_15', 1)
else
	conf_data.set('D2TK_INPUT_1_15', 0)
endif

if use_fontconfig.enabled()
	conf_data.set('D2TK_FONTCONFIG', 1)
else
	conf_data.set('D2TK_FONTCONFIG', 0)
endif

example_srcs = [
	join_paths('example', 'example.c')
]

example_cairo_srcs = [
	join_paths('example', 'custom_cairo.c')
]

example_nanovg_srcs = [
	join_paths('example', 'custom_nanovg.c')
]

example_pugl_srcs = [
	join_paths('example', 'd2tk_pugl.c')
]

example_fbdev_srcs = [
	join_paths('example', 'd2tk_fbdev.c')
]

example_glfw_srcs = [
	join_paths('example', 'd2tk_glfw.c')
]

pugl_srcs = [
	join_paths('src', 'frontend_pugl.c'),
	join_paths('pugl', 'src', 'implementation.c')
]

pugl_gl_srcs = []

pugl_cairo_srcs = []

nanovg_srcs = [
	join_paths('nanovg', 'src', 'nanovg.c'),
	join_paths('src', 'backend_nanovg.c')
]

cairo_srcs = [
	join_paths('src', 'backend_cairo.c')
]

fbdev_srcs = [
	join_paths('src', 'frontend_fbdev.c')
]

glfw_srcs = [
	join_paths('src', 'frontend_glfw.c')
]

test_core_srcs = [
	join_paths('test', 'core.c'),
	join_paths('test', 'mock.c')
]

test_base_srcs = [
	join_paths('test', 'base.c'),
	join_paths('test', 'mock.c')
]

c_args = ['-fvisibility=hidden',
	'-ffast-math']

if host_machine.system() == 'windows'
	deps += cc.find_library('opengl32', required : use_frontend_pugl)
	deps += cc.find_library('gdi32', required : use_frontend_pugl)
	deps += cc.find_library('ws2_32', required : true)
	pugl_srcs += join_paths('pugl', 'src', 'win.c')
	pugl_gl_srcs += join_paths('pugl', 'src', 'win_gl.c')
	pugl_cairo_srcs += join_paths('pugl', 'src', 'win_cairo.c')
elif host_machine.system() == 'darwin'
	add_languages('objc')
	links += ['-framework', 'OpenGL']
	links += ['-framework', 'Cocoa']
	pugl_srcs += join_paths('pugl','src', 'mac.m')
	pugl_gl_srcs += join_paths('pugl', 'src', 'mac_gl.m')
	pugl_cairo_srcs += join_paths('pugl', 'src', 'mac_cairo.m')
else
	deps += dependency('gl', required : use_frontend_pugl)
	deps += dependency('x11', version : '>=1.6.0', required : use_frontend_pugl)
	deps += dependency('xext', version : '>=1.3.0', required : use_frontend_pugl)
	pugl_srcs += join_paths('pugl', 'src', 'x11.c')
	pugl_srcs += join_paths('pugl', 'src', 'x11_stub.c')
	pugl_gl_srcs += join_paths('pugl', 'src', 'x11_gl.c')
	pugl_cairo_srcs += join_paths('pugl', 'src', 'x11_cairo.c')
endif

if use_backend_cairo.enabled()
	if use_frontend_pugl.enabled()
		d2tk_cairo = declare_dependency(
			compile_args : ['-DPUGL_HAVE_CAIRO', '-DPUGL_STATIC'],
			include_directories : inc_dir,
			dependencies : [deps, cairo_deps, cairo_xlib_dep],
			link_args : links,
			sources : [lib_srcs, cairo_srcs, pugl_srcs, pugl_cairo_srcs])

		if build_examples
			executable('d2tk.cairo', [example_srcs, example_pugl_srcs, example_cairo_srcs],
				c_args : c_args,
				include_directories : inc_dir,
				dependencies: d2tk_cairo,
				install : false)
		endif
	endif

	if use_frontend_fbdev.enabled()
		d2tk_fbdev = declare_dependency(
			include_directories : inc_dir,
			dependencies : [deps, cairo_deps, input_dep, udev_dep, evdev_dep],
			link_args : links,
			sources : [lib_srcs, cairo_srcs, fbdev_srcs])

		if build_examples
			executable('d2tk.fbdev', [example_srcs, example_fbdev_srcs, example_cairo_srcs],
				c_args : c_args,
				include_directories : inc_dir,
				dependencies: d2tk_fbdev,
				install : false)
		endif
	endif
endif

if use_backend_nanovg.enabled()
	if use_frontend_pugl.enabled()
		d2tk_nanovg = declare_dependency(
			compile_args : ['-DPUGL_STATIC'],
			include_directories : inc_dir,
			dependencies : [deps, glu_dep, glew_dep],
			link_args : links,
			sources : [lib_srcs, nanovg_srcs, pugl_srcs, pugl_gl_srcs])

		if build_examples
			executable('d2tk.nanovg', [example_srcs, example_pugl_srcs, example_nanovg_srcs],
				c_args : c_args,
				include_directories : inc_dir,
				dependencies: d2tk_nanovg,
				install : false)
		endif
	endif

	if use_frontend_glfw.enabled()
		d2tk_glfw = declare_dependency(
			include_directories : inc_dir,
			dependencies : [deps, glfw_dep, glew_dep],
			link_args : links,
			sources : [lib_srcs, nanovg_srcs, glfw_srcs])

		if build_examples
			executable('d2tk.glfw', [example_srcs, example_glfw_srcs, example_nanovg_srcs],
				c_args : c_args,
				include_directories : inc_dir,
				dependencies: d2tk_glfw,
				install : false)
		endif
	endif
endif

config_h = configure_file(
	input : join_paths('d2tk', 'config.h.in'),
	output : 'config.h',
	configuration : conf_data,
	install : false)

if not use_fontconfig.enabled()
	fira_sans_bold_ttf = configure_file(
		input : join_paths('ttf', 'FiraSans-Bold.ttf'),
		output : 'FiraSans:bold.ttf',
		copy : true,
		install : false)

	fira_code_bold_ttf = configure_file(
		input : join_paths('ttf', 'FiraCode-Bold.ttf'),
		output : 'FiraCode:bold.ttf',
		copy : true,
		install : false)

	fira_code_light_ttf = configure_file(
		input : join_paths('ttf', 'FiraCode-Light.ttf'),
		output : 'FiraCode:light.ttf',
		copy : true,
		install : false)

	fira_code_medium_ttf = configure_file(
		input : join_paths('ttf', 'FiraCode-Medium.ttf'),
		output : 'FiraCode:medium.ttf',
		copy : true,
		install : false)

	fira_code_regular_ttf = configure_file(
		input : join_paths('ttf', 'FiraCode-Regular.ttf'),
		output : 'FiraCode:regular.ttf',
		copy : true,
		install : false)
endif

if build_examples
	configure_file(
		input : join_paths('example', 'libre-arrow-circle-right.png'),
		output : 'libre-arrow-circle-right.png',
		copy : true,
		install : false)

	configure_file(
		input : join_paths('example', 'libre-gui-folder.png'),
		output : 'libre-gui-folder.png',
		copy : true,
		install : false)

	configure_file(
		input : join_paths('example', 'libre-gui-file.png'),
		output : 'libre-gui-file.png',
		copy : true,
		install : false)
endif

if build_tests
	test_core = executable('test.core', [test_core_srcs, lib_srcs],
		c_args : c_args,
		dependencies : deps,
		include_directories : inc_dir,
		install : false)

	test_base = executable('test.base', [test_base_srcs, lib_srcs],
		c_args : c_args,
		dependencies : deps,
		include_directories : inc_dir,
		install : false)

	test('Test core', test_core)
	test('Test base', test_base)

	if fc_list.found() and grep.found() and check_for_font.found()
		test('FiraSans-Bold.ttf', check_for_font, args : ['FiraSans-Bold.ttf'])
		test('FiraCode-Light.ttf', check_for_font, args : ['FiraCode-Light.tt'])
		test('FiraCode-Regular.ttf', check_for_font, args : ['FiraCode-Regular.ttf'])
		test('FiraCode-Medium.ttf', check_for_font, args : ['FiraCode-Medium.ttf'])
		test('FiraCode-Bold.ttf', check_for_font, args : ['FiraCode-Bold.ttf'])
	endif
endif

if reuse.found()
  test('REUSE', reuse, args : [
    '--root', meson.current_source_dir(),
    'lint'
  ])
endif

if build_doc
	run_command('rm', '-rf', 'doc',
    check : true)

	run_command(robodoc,
		'--src', './d2tk',
		'--doc', 'doc',
		'--multidoc', '--troff', '--nosort', '--nodesc', '--cmode',
		'--compress', 'gzip',
    check : true)
endif
