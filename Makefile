CXX = g++

CXX_FLAGS = -Wall -std=c++14 -Wall -march=native -fopenmp

INSTALL_DIR = $(HOME)/fabrica/local

CXX_FLAGS += -I$(INSTALL_DIR)/include -L$(INSTALL_DIR)/lib

CXX_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

LIBS   = -ldiaplexis -lsynarmosma

euplecton: euplecton.cpp
	$(CXX) $(CXX_FLAGS) -o euplecton euplecton.cpp $(LIBS)

clean:
	rm -f euplecton
	rm -f *~










