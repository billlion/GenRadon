cmake_minimum_required(VERSION 3.16)
project(GenRadon VERSION 0.1.0 LANGUAGES CXX)
include(GNUInstallDirs)

add_library(genradon INTERFACE)
add_library(GenRadon::genradon ALIAS genradon)
target_compile_features(genradon INTERFACE cxx_std_17)
target_include_directories(genradon INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

option(GENRADON_BUILD_EXAMPLES "Build small synthetic acquisitions" ON)
include(CTest)
if(BUILD_TESTING)
    add_executable(genradon_tests tests/test_core.cpp)
    target_link_libraries(genradon_tests PRIVATE GenRadon::genradon)
    add_test(NAME core COMMAND genradon_tests)
endif()
if(GENRADON_BUILD_EXAMPLES)
    add_executable(genradon_toy examples/toy_transforms.cpp)
    target_link_libraries(genradon_toy PRIVATE GenRadon::genradon)
    if(BUILD_TESTING)
        add_test(NAME toy_midpoint COMMAND genradon_toy --smoke)
        add_test(NAME toy_affine COMMAND genradon_toy --smoke --affine)
    endif()
endif()

option(GENRADON_TEST_CODEGEN "Generate and test the SymPy example (requires Python and SymPy)" OFF)
if(GENRADON_TEST_CODEGEN)
    if(NOT BUILD_TESTING)
        message(FATAL_ERROR "GENRADON_TEST_CODEGEN requires BUILD_TESTING")
    endif()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generated/spheroid3.hpp
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_spheroid.py
                --output ${CMAKE_CURRENT_BINARY_DIR}/generated/spheroid3.hpp
        DEPENDS tools/generate_spheroid.py VERBATIM)
    add_executable(genradon_codegen_test tests/test_generated.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/generated/spheroid3.hpp)
    target_include_directories(genradon_codegen_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
    target_link_libraries(genradon_codegen_test PRIVATE GenRadon::genradon)
    add_test(NAME generated_spheroid COMMAND genradon_codegen_test)
endif()

include(CMakePackageConfigHelpers)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS genradon EXPORT GenRadonTargets)
install(EXPORT GenRadonTargets NAMESPACE GenRadon:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GenRadon)
configure_package_config_file(cmake/GenRadonConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GenRadonConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GenRadon)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/GenRadonConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GenRadonConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/GenRadonConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GenRadon)
