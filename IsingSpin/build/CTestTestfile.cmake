# CMake generated Testfile for 
# Source directory: /home/peggy_piggy_zhu/SJTU_CP2024/IsingSpin
# Build directory: /home/peggy_piggy_zhu/SJTU_CP2024/IsingSpin/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTests_IsingSpin "UnitTests_catch2")
set_tests_properties(UnitTests_IsingSpin PROPERTIES  _BACKTRACE_TRIPLES "/home/peggy_piggy_zhu/SJTU_CP2024/IsingSpin/CMakeLists.txt;21;add_test;/home/peggy_piggy_zhu/SJTU_CP2024/IsingSpin/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
