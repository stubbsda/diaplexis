include ../../Makefile.config

HEADERS    = event.h simplex.h complex.h spacetime.h sheet.h

OBJECTS = spacetime.o spacetime_disk.o spacetime_visualization.o spacetime_optimization.o spacetime_geometry.o \
spacetime_hyphansis.o spacetime_initialization.o spacetime_explication.o spacetime_implication.o event.o complex.o \
spacetime_topology.o complex_topology.o complex_dynamics.o simplex.o sheet.o

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

diaplexis: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o libdiaplexis.so $(OBJECTS) $(LIBS)

spacetime.o: spacetime.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime.cxx

spacetime_disk.o: spacetime_disk.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_disk.cxx

spacetime_optimization.o: spacetime_optimization.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_optimization.cxx

spacetime_explication.o: spacetime_explication.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_explication.cxx

spacetime_implication.o: spacetime_implication.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_implication.cxx

spacetime_hyphansis.o: spacetime_hyphansis.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_hyphansis.cxx

spacetime_topology.o: spacetime_topology.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_topology.cxx

spacetime_geometry.o: spacetime_geometry.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_geometry.cxx

spacetime_visualization.o: spacetime_visualization.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_visualization.cxx

spacetime_initialization.o: spacetime_initialization.cxx $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_initialization.cxx

complex.o: complex.cxx event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex.cxx

complex_dynamics.o: complex_dynamics.cxx event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex_dynamics.cxx

complex_topology.o: complex_topology.cxx event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex_topology.cxx

sheet.o: sheet.cxx sheet.h
	$(CXX) $(CXX_FLAGS) -c sheet.cxx

simplex.o: simplex.cxx simplex.h
	$(CXX) $(CXX_FLAGS) -c simplex.cxx

event.o: event.cxx event.h
	$(CXX) $(CXX_FLAGS) -c event.cxx

clean:
	rm -f *~
	rm -f $(OBJECTS)
	rm -f libdiaplexis.so










