-- Global rules
add_rules("mode.debug", "mode.release")

-- Package Requirements
add_requires("openmp")

-- Target Definition
target("bump")
    set_kind("binary")
    set_languages("c++20")
    set_optimize("fastest")
    add_files("src/*.cpp")
    add_packages("openmp")
