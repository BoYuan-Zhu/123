# CMake generated Testfile for 
# Source directory: /home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem_Square
# Build directory: /home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem_Square/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[UnitTests_Square_Lattice]=] "UnitTests_catch2")
set_tests_properties([=[UnitTests_Square_Lattice]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem_Square/CMakeLists.txt;26;add_test;/home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem_Square/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
