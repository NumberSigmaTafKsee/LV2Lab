# SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
# SPDX-License-Identifier: CC0-1.0

project('timely.lv2', 'c', default_options : [
	'buildtype=release',
	'warning_level=3',
	'werror=false',
	'b_lto=true',
	'c_std=c11'])

lv2libdir = get_option('lv2libdir')

inst_dir = join_paths(lv2libdir, meson.project_name())

reuse = find_program('reuse', required : false)
cp = find_program('cp')
lv2_validate = find_program('lv2_validate', native : true, required : false)
sord_validate = find_program('sord_validate', native : true, required : false)
lv2lint = find_program('lv2lint', required : false)

add_project_arguments('-D_GNU_SOURCE', language : 'c')

conf_data = configuration_data()
cc = meson.get_compiler('c')

m_dep = cc.find_library('m')
lv2_dep = dependency('lv2', version : '>=1.14.0')

lib_deps = []
lib_deps += m_dep
lib_deps += lv2_dep

lib_incs = []
lib_incs += include_directories('')

lib_srcs = []
lib_srcs += join_paths('src', 'timely.c')

timely_lv2 = declare_dependency(
  include_directories : lib_incs,
  dependencies : lib_deps,
  sources : lib_srcs)

test_args = []
test_args += '-fvisibility=hidden'
test_args += '-ffast-math'

test_deps = []
test_deps += timely_lv2

test_srcs = []
test_srcs += join_paths('test', 'timely.c')

install = not meson.is_subproject()

mod = shared_module('timely', test_srcs,
	c_args : test_args,
	name_prefix : '',
	dependencies : test_deps,
	install : install,
	install_dir : inst_dir)

version = get_option('version').split('.')
conf_data.set('MAJOR_VERSION', version[0])
conf_data.set('MINOR_VERSION', version[1])
conf_data.set('MICRO_VERSION', version[2])

suffix = mod.full_path().strip().split('.')[-1]
conf_data.set('MODULE_SUFFIX', '.' + suffix)

manifest_ttl = configure_file(
  input : join_paths('test', 'manifest.ttl.in'),
  output : 'manifest.ttl',
	configuration : conf_data,
	install : install,
	install_dir : inst_dir)

dsp_ttl = configure_file(
	input : join_paths('test', 'timely.ttl'),
	output : 'timely.ttl',
	copy : true,
	install : install,
	install_dir : inst_dir)

if lv2_validate.found() and sord_validate.found()
	test('LV2 validate', lv2_validate,
		args : [manifest_ttl, dsp_ttl])
endif

if lv2lint.found()
	test('LV2 lint', lv2lint,
    args : ['-M', 'pack',
      '-E', 'warn',
      '-I', meson.current_build_dir(),
			'http://open-music-kontrollers.ch/lv2/timely#test'])
endif

if reuse.found()
  test('REUSE', reuse, args : [
    '--root', meson.current_source_dir(),
    'lint'
  ])
endif
