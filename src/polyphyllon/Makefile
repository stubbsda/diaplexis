include ../../Makefile.config

HEADERS    = event.h simplex.h complex.h spacetime.h sheet.h

SOURCES = spacetime.cxx spacetime_disk.cxx spacetime_initialization.cxx spacetime_hyphansis.cxx \
spacetime_geometry.cxx spacetime_topology.cxx spacetime_explication.cxx spacetime_implication.cxx \
spacetime_optimization.cxx spacetime_visualization.cxx

OBJECTS = event_wrapper.o simplex.o complex_wrapper.o spacetime_wrapper.o sheet.o

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

diaplexis: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o libdiaplexis.so $(OBJECTS) $(LIBS)

spacetime_wrapper.o: spacetime_wrapper.cxx $(SOURCES) $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_wrapper.cxx

complex_wrapper.o: complex_wrapper.cxx complex.cxx complex_topology.cxx complex_dynamics.cxx event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex_wrapper.cxx

sheet.o: sheet.cxx sheet.h
	$(CXX) $(CXX_FLAGS) -c sheet.cxx

simplex.o: simplex.cxx simplex.h
	$(CXX) $(CXX_FLAGS) -c simplex.cxx

event_wrapper.o: event_wrapper.cxx event.cxx event.h
	$(CXX) $(CXX_FLAGS) -c event_wrapper.cxx

clean:
	rm -f *~
	rm -f $(OBJECTS)
	rm -f libdiaplexis.so










