-- Global rules
add_rules("mode.debug", "mode.release")

-- Target Definition
target("bump")
    set_kind("binary")
    set_languages("c++20")
    set_optimize("fastest")
    add_files("src/*.cpp")
