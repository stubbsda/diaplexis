CXX_FLAGS += -I$(BIBLIOTHEK)/include $(OPENMP) -fPIC #-DDISCRETE
LD_FLAGS  += -L$(BIBLIOTHEK)/lib $(OPENMP) -shared

#CXX_FLAGS += $(CXX_OPT)
CXX_FLAGS += $(DEBUG) 

#LD_FLAGS += $(CXX_OPT)
LD_FLAGS  += $(DEBUG)

UNAME = $(shell uname)

ifeq ($(UNAME),Darwin)
  LD_FLAGS += -install_name $(BIBLIOTHEK)/lib/libdiaplexis.so
endif

HEADERS    = event.h simplex.h spacetime.h

LIBS       = -lsynarmosma $(LAPACK) $(BOOST_TIMER) $(BOOST_FILESYSTEM) $(BOOST_SYSTEM) -lpugixml -lntl -lgmp -lm

OBJECTS = spacetime.o spacetime_disk.o spacetime_visualize.o spacetime_dynamics.o spacetime_topology.o \
spacetime_geometry.o spacetime_hyphansis.o event.o simplex.o

install: diaplexis 
	mkdir -p $(BIBLIOTHEK)/lib
	install -p libdiaplexis.so $(BIBLIOTHEK)/lib/
	mkdir -p $(BIBLIOTHEK)/include/diaplexis
	install -p -m 444 *.h $(BIBLIOTHEK)/include/diaplexis/

diaplexis: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o libdiaplexis.so $(OBJECTS) $(LIBS)

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

simplex.o: simplex.cpp simplex.h
	$(CXX) $(CXX_FLAGS) -c simplex.cpp

event.o: event.cpp event.h
	$(CXX) $(CXX_FLAGS) -c event.cpp

clean:
	rm -f *~
	rm -f $(OBJECTS)
	rm -f libdiaplexis.so
	rm -f $(BIBLIOTHEK)/lib/libdiaplexis.so
	rm -rf $(BIBLIOTHEK)/include/diaplexis









