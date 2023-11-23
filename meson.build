project(
  'fftpack',
  'fortran',
  version: '4.0.0',
  meson_version: '>=0.57',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
    'c_std=c11',
  ],
)

install = not (meson.is_subproject() and get_option('default_library') == 'static')

# Collect source of the project
srcs = []
subdir('src')

# General configuration information
inc_dirs = []
lib_deps = []
subdir('configure')

# Library target
fftpack_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
fftpack_dep = declare_dependency(
  link_with: fftpack_lib,
  dependencies: lib_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)

# add the testsuite
subdir('test')

if install
  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('configure'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    fftpack_lib,
    description: 'fftpack',
    subdirs: ['', module_id],
  )
endif
