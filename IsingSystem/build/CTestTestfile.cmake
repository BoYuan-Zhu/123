# CMake generated Testfile for 
# Source directory: /home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem
# Build directory: /home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest_IsingSystem "UnitTests_catch2")
set_tests_properties(UnitTest_IsingSystem PROPERTIES  _BACKTRACE_TRIPLES "/home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem/CMakeLists.txt;21;add_test;/home/peggy_piggy_zhu/SJTU_CP2024/IsingSystem/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
