CXX = g++

CXX_FLAGS = -std=c++14 -Wall -march=native

LD_FLAGS = -Wall -fopenmp

INSTALL_DIR = $(HOME)/fabrica/local

CXX_FLAGS += -I$(INSTALL_DIR)/include

LD_FLAGS += -L$(INSTALL_DIR)/lib

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

LAPACK = -llapack -lblas

BOOST = -lboost_timer -lboost_filesystem -lboost_system

LIBS   = -ldiaplexis -lsynarmosma $(LAPACK) $(BOOST) -lpugixml -lntl -lgmp -lm

euplecton: euplecton.o 
	$(CXX) $(LD_FLAGS) -o euplecton euplecton.o $(LIBS)

euplecton.o: euplecton.cpp 
	$(CXX) $(CXX_FLAGS) -c euplecton.cpp

clean:
	rm -f euplecton euplecton.o 
	rm -f *~










