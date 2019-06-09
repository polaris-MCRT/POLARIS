set(CMAKE_CXX_STANDARD 98)
set(CXX_STANDARD_REQUIRED ON)

if(BUILD_STATIC_LIBS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98 -O2 -g -fopenmp -MMD -MP -MF -static-libgcc -static-libstdc++ -static")
else(BUILD_STATIC_LIBS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98 -O2 -g -fopenmp -MMD -MP -MF")
endif(BUILD_STATIC_LIBS)

SET(SRC_FILES CommandParser.cpp Cylindrical.cpp Dust.cpp 
    Faddeeva.cc Grid.cpp MathFunctions.cpp GasSpecies.cpp  
    OcTree.cpp Pipeline.cpp RadiativeTransfer.cpp Source.cpp 
    Spherical.cpp Voronoi.cpp main.cpp Synchrotron.cpp
)

add_executable (polaris ${SRC_FILES})

target_link_libraries(polaris ${CCFITS_LIBRARY})
target_link_libraries(polaris ${CFITSIO_LIBRARY})

install(TARGETS polaris
        DESTINATION ${CMAKE_CURRENT_LIST_DIR}/../bin/)
