cmake_minimum_required(VERSION 3.15)
project(AIDAProject CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/bin)

add_compile_definitions(_USE_MATH_DEFINES WITH_ALGEBRAICMATRIX)

# -------------------------------------------------------------
# INCLUDE DIRECTORIES
# -------------------------------------------------------------
set(DACE_INCLUDE     /usr/local/include/dace)
set(CSPICE_INCLUDE   /usr/local/include/cspice)
set(ASTRO_INCLUDE    /usr/local/include/astro)
set(ATMOS_INCLUDE    /usr/local/include/atmos)
set(DYNORB_INCLUDE   /usr/local/include/dynorb)
set(GECO_INCLUDE     /usr/local/include/geco)
set(ODEPROP_INCLUDE  /usr/local/include/odeprop)
set(SADA_INCLUDE     /usr/local/include/sada)
set(JSONCPP_INCLUDE  /usr/include/jsoncpp)
set(NLOHMANN_INCLUDE /usr/include/nlohmann)

# -------------------------------------------------------------
# FIND LIBRARIES
# -------------------------------------------------------------
set(LIB_HINTS /usr/local/lib /usr/lib /usr/lib/x86_64-linux-gnu)

find_library(DACE_LIB    NAMES dace    HINTS ${LIB_HINTS})
find_library(ASTRO_LIB   NAMES astro   HINTS ${LIB_HINTS})
find_library(ATMOS_LIB   NAMES atmos   HINTS ${LIB_HINTS})
find_library(GECO_LIB    NAMES geco    HINTS ${LIB_HINTS})
find_library(DYNORB_LIB  NAMES dynorb  HINTS ${LIB_HINTS})
find_library(ODEPROP_LIB NAMES odeprop HINTS ${LIB_HINTS})
find_library(SADA_LIB    NAMES sada    HINTS ${LIB_HINTS})
find_library(JSONCPP_LIB NAMES jsoncpp HINTS ${LIB_HINTS})

# CSPICE static libs
set(CSPICE_LIB     /usr/local/lib/cspice.a)
set(CSUPPORT_LIB   /usr/local/lib/csupport.a)

# -------------------------------------------------------------
# COLLECT ALL NON-CSPICE LIBS FIRST
# -------------------------------------------------------------
set(BASE_LIBS m)

if (DACE_LIB)     
	list(APPEND BASE_LIBS ${DACE_LIB})     
endif()
if (ASTRO_LIB)    
	list(APPEND BASE_LIBS ${ASTRO_LIB})    
endif()
if (ATMOS_LIB)    
	list(APPEND BASE_LIBS ${ATMOS_LIB})    
endif()
if (DYNORB_LIB)   
	list(APPEND BASE_LIBS ${DYNORB_LIB})   
endif()
if (GECO_LIB)     
	list(APPEND BASE_LIBS ${GECO_LIB})     
endif()
if (ODEPROP_LIB)  
	list(APPEND BASE_LIBS ${ODEPROP_LIB})  
endif()
if (SADA_LIB)     
	list(APPEND BASE_LIBS ${SADA_LIB})     
endif()
if (JSONCPP_LIB)  
	list(APPEND BASE_LIBS ${JSONCPP_LIB})  
endif()

# -------------------------------------------------------------
# HELPERS — configure target
# -------------------------------------------------------------
function(configure_target tgt)
    target_include_directories(${tgt} BEFORE PRIVATE
        ${DACE_INCLUDE}
        ${CSPICE_INCLUDE}
        ${ASTRO_INCLUDE}
        ${ATMOS_INCLUDE}
        ${DYNORB_INCLUDE}
        ${GECO_INCLUDE}
        ${ODEPROP_INCLUDE}
        ${SADA_INCLUDE}
        ${JSONCPP_INCLUDE}
        ${NLOHMANN_INCLUDE}
        /usr/local/include
        /usr/include   
	)

    # Link non-CSPICE libs FIRST
    target_link_libraries(${tgt} PRIVATE ${BASE_LIBS})

    # Then CSPICE last (static library ordering requirement)
    target_link_libraries(${tgt} PRIVATE
        ${CSUPPORT_LIB}
        ${CSPICE_LIB}
    )

    # Set RPATH
    set_target_properties(${tgt} PROPERTIES
        BUILD_RPATH "/usr/local/lib"
        INSTALL_RPATH "/usr/local/lib"
    )
endfunction()

# -------------------------------------------------------------
# EXECUTABLES
# -------------------------------------------------------------
add_executable(propAida src/propAida.cpp src/prop_utils.cpp)
configure_target(propAida)

add_executable(propAidaDebris src/propAidaDebris.cpp src/prop_utils.cpp)
configure_target(propAidaDebris)

add_executable(nliGmm src/nliGmm.cpp src/prop_utils.cpp)
configure_target(nliGmm)

add_executable(findRealTca src/findRealTca.cpp src/prop_utils.cpp)
configure_target(findRealTca)

add_executable(pocMaps src/pocMaps.cpp src/prop_utils.cpp)
configure_target(pocMaps)

add_executable(findGeo src/findGeo.cpp src/prop_utils.cpp)
configure_target(findGeo)