#
# Meson Build file for the Bump project
# Panos Asproulis <p.asproulis@icloud.com>
#
project(
    'Bump',
    'cpp',
    version: '1.0.0',
    default_options: ['cpp_std=c++20', 'buildtype=debugoptimized'],
    license: 'MIT',
    meson_version: '>=0.50.0',
)
#
# Source files
#
src = ['src/driver.cpp', 'src/Bump.cpp', 'src/ShallowWater.cpp']
#
# Include directories
#
hdrs = include_directories('src')
#
# Documentation generation
#
doxygen = find_program('doxygen', required: false)
if doxygen.found()
    message('Doxygen found')
    run_target('docs', command: [doxygen, meson.source_root() + '/Doxyfile'])
else
    warning('Documentation disabled without doxygen')
endif
#
# Executable target
#
executable('bump', src, include_directories: hdrs)