# Which C++ compiler to use, default is Gnu g++
CXX         = g++
# Where is the Synarmosma library located?
SYNARMOSMA  = $(HOME)/fabrica/local
# General flags for C++ compilation, assumes compiler with OpenMP support
CXX_FLAGS   = -Wall -std=c++14 -fopenmp -march=native -I$(SYNARMOSMA)/include -fPIC
# For debug mode
DEBUG       = -g -DDEBUG
# Optimization values for the g++ on an x86-64 platform
OPT         = -O3 -fstrict-aliasing -ftree-vectorize -funroll-loops 
# Choose either the debug or optimized version of the library
#CXX_FLAGS  += $(OPT)
CXX_FLAGS  += $(DEBUG) 
# Verbose output 
CXX_FLAGS  += -DVERBOSE 
# If the Synarmosma library was compiled with the DISCRETE flag, so must 
# the Diaplexis library
#CXX_FLAGS += -DDISCRETE
# Standard linking flags, shouldn't need to be changed
LD_FLAGS    = -L$(SYNARMOSMA)/lib -Wall -shared -fopenmp
# Debug or optimized version for linking...
#LD_FLAGS   += $(OPT)
LD_FLAGS   += $(DEBUG)
# Name of the LAPACK and BLAS libraries on this platform
LAPACK      = -llapack -lblas
# The Boost library on this platform
BOOST       = -lboost_timer -lboost_filesystem -lboost_system
# All the necessary libraries: Synarmosma, LAPACK, Boost, PugiXML and NTL
LIBS        = -lsynarmosma $(LAPACK) $(BOOST) -lpugixml -lntl -lgmp -lm
# Where to install the Diaplexis library, requires root privileges
INSTALL_DIR = $(HOME)/fabrica/local
# End of user-modifiable parameters!

HEADERS    = event.h simplex.h complex.h spacetime.h

OBJECTS = spacetime.o spacetime_disk.o spacetime_visualization.o spacetime_optimization.o spacetime_geometry.o \
spacetime_hyphansis.o spacetime_initialization.o spacetime_explication.o spacetime_implication.o event.o complex.o \
spacetime_topology.o complex_topology.o complex_dynamics.o simplex.o

install: diaplexis 
	mkdir -p $(INSTALL_DIR)/lib
	install -p libdiaplexis.so $(INSTALL_DIR)/lib/
	mkdir -p $(INSTALL_DIR)/include/diaplexis
	install -p -m 444 *.h $(INSTALL_DIR)/include/diaplexis/

diaplexis: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o libdiaplexis.so $(OBJECTS) $(LIBS)

spacetime.o: spacetime.cpp $(HEADERS) 
	$(CXX) $(CXX_FLAGS) -c spacetime.cpp

spacetime_disk.o: spacetime_disk.cpp $(HEADERS) 
	$(CXX) $(CXX_FLAGS) -c spacetime_disk.cpp

spacetime_optimization.o: spacetime_optimization.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_optimization.cpp

spacetime_explication.o: spacetime_explication.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_explication.cpp

spacetime_implication.o: spacetime_implication.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_implication.cpp

spacetime_topology.o: spacetime_topology.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_topology.cpp

spacetime_visualization.o: spacetime_visualization.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_visualization.cpp

spacetime_initialization.o: spacetime_initialization.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_initialization.cpp

spacetime_hyphansis.o: spacetime_hyphansis.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_hyphansis.cpp

spacetime_geometry.o: spacetime_geometry.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c spacetime_geometry.cpp

complex.o: complex.cpp event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex.cpp

complex_dynamics.o: complex_dynamics.cpp event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex_dynamics.cpp

complex_topology.o: complex_topology.cpp event.h simplex.h complex.h
	$(CXX) $(CXX_FLAGS) -c complex_topology.cpp

simplex.o: simplex.cpp simplex.h
	$(CXX) $(CXX_FLAGS) -c simplex.cpp

event.o: event.cpp event.h
	$(CXX) $(CXX_FLAGS) -c event.cpp

clean:
	rm -f *~
	rm -f $(OBJECTS)
	rm -f libdiaplexis.so
	rm -f $(INSTALL_DIR)/lib/libdiaplexis.so
	rm -rf $(INSTALL_DIR)/include/diaplexis









