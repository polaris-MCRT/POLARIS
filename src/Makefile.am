bin_PROGRAMS = polaris
CLEANFILES = $(bin_PROGRAMS)

polaris_CXXFLAGS = -fopenmp -MMD -MP -MF

polaris_SOURCES = CommandParser.cpp CommandParser.h Cylindrical.cpp Cylindrical.h Detector.h Dust.cpp Dust.h Faddeeva.cc Faddeeva.hhh Grid.cpp Grid.h MathFunctions.cpp MathFunctions.h Matrix2D.h GasSpecies.cpp GasSpecies.h OcTree.cpp OcTree.h Pipeline.cpp Pipeline.h RadiativeTransfer.cpp RadiativeTransfer.h Source.cpp Source.h Spherical.cpp Spherical.h Stokes.h Vector.h Voronoi.cpp Voronoi.h chelper.h main.cpp typedefs.h Raytracing.h Synchrotron.cpp Synchrotron.h

clean::
	rm -rf polaris
	rm -rf *.o
