# CMake generated Testfile for 
# Source directory: C:/dev/plate-tectonics/test
# Build directory: C:/dev/plate-tectonics/build-png/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(PlateTectonicsTests "C:/dev/plate-tectonics/build-png/test/Debug/PlateTectonicsTests.exe")
  set_tests_properties(PlateTectonicsTests PROPERTIES  WORKING_DIRECTORY "C:/dev/plate-tectonics/build-png" _BACKTRACE_TRIPLES "C:/dev/plate-tectonics/test/CMakeLists.txt;21;add_test;C:/dev/plate-tectonics/test/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(PlateTectonicsTests "C:/dev/plate-tectonics/build-png/test/Release/PlateTectonicsTests.exe")
  set_tests_properties(PlateTectonicsTests PROPERTIES  WORKING_DIRECTORY "C:/dev/plate-tectonics/build-png" _BACKTRACE_TRIPLES "C:/dev/plate-tectonics/test/CMakeLists.txt;21;add_test;C:/dev/plate-tectonics/test/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(PlateTectonicsTests "C:/dev/plate-tectonics/build-png/test/MinSizeRel/PlateTectonicsTests.exe")
  set_tests_properties(PlateTectonicsTests PROPERTIES  WORKING_DIRECTORY "C:/dev/plate-tectonics/build-png" _BACKTRACE_TRIPLES "C:/dev/plate-tectonics/test/CMakeLists.txt;21;add_test;C:/dev/plate-tectonics/test/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(PlateTectonicsTests "C:/dev/plate-tectonics/build-png/test/RelWithDebInfo/PlateTectonicsTests.exe")
  set_tests_properties(PlateTectonicsTests PROPERTIES  WORKING_DIRECTORY "C:/dev/plate-tectonics/build-png" _BACKTRACE_TRIPLES "C:/dev/plate-tectonics/test/CMakeLists.txt;21;add_test;C:/dev/plate-tectonics/test/CMakeLists.txt;0;")
else()
  add_test(PlateTectonicsTests NOT_AVAILABLE)
endif()
subdirs("../_deps/googletest-build")
