# SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
# SPDX-License-Identifier: CC0-1.0

project('ser_atom.lv2', 'c', default_options : [
	'buildtype=release',
	'warning_level=3',
	'werror=false',
	'b_lto=true',
	'c_std=c11'])

reuse = find_program('reuse', required : false)

add_project_arguments('-D_GNU_SOURCE', language : 'c')

lv2_dep = dependency('lv2')

lib_deps = []
lib_deps += lv2_dep

lib_incs = []
lib_incs += include_directories('')

lib_srcs = []
lib_srcs += join_paths('src', 'ser_atom.c')

ser_atom_lv2 = declare_dependency(
  include_directories : lib_incs,
  dependencies : lib_deps,
  sources : lib_srcs)

test_deps = []
test_deps += ser_atom_lv2

test_srcs = []
test_srcs += join_paths('test', 'test.c')

test_ser_atom = executable('Test', test_srcs,
	dependencies : test_deps,
	install : false)

test('Test', test_ser_atom)

if reuse.found()
  test('REUSE', reuse, args : [
    '--root', meson.current_source_dir(),
    'lint'
  ])
endif
