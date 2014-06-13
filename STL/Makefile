include ../Makefile.def

OBJECTS = sheet.o spacetime.o spacetime_disk.o spacetime_visualize.o spacetime_dynamics.o spacetime_topology.o \
spacetime_geometry.o spacetime_hyphansis.o vertex.o simplex.o

TARGET = libdiaplexis-stl.so

diaplexis: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

spacetime.o: spacetime.cpp $(HEADERS) 
	$(CXX) $(CXX_FLAGS) -c spacetime.cpp

spacetime_disk.o: spacetime_disk.cpp $(HEADERS) 
	$(CXX) $(CXX_FLAGS) -c spacetime_disk.cpp

spacetime_dynamics.o: spacetime_dynamics.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_dynamics.cpp

spacetime_visualize.o: spacetime_visualize.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_visualize.cpp

spacetime_topology.o: spacetime_topology.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_topology.cpp

spacetime_hyphansis.o: spacetime_hyphansis.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_hyphansis.cpp

spacetime_geometry.o: spacetime_geometry.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_geometry.cpp

sheet.o: sheet.cpp sheet.h 
	$(CXX) $(CXX_FLAGS) -c sheet.cpp

simplex.o: simplex.cpp simplex.h
	$(CXX) $(CXX_FLAGS) -c simplex.cpp

vertex.o: vertex.cpp vertex.h
	$(CXX) $(CXX_FLAGS) -c vertex.cpp

clean:
	rm -f $(OBJECTS)
	rm -f $(TARGET) 
	rm -f *~










