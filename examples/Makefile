include ../Makefile.config

CXX_FLAGS += -I$(INSTALL_DIR)/include -L$(INSTALL_DIR)/lib -Wl,-rpath $(INSTALL_DIR)/lib

ifeq ($(TYPE),MONOPHYLLON)
  CXX_FLAGS += -DMONOPHYLLON
  LIBS = -ldiaplexis-monophyllon
else 
  CXX_FLAGS += -DPOLYPHYLLON
  LIBS = -ldiaplexis-polyphyllon
endif

LIBS += -lsynarmosma

euplecton:
	$(CXX) $(CXX_FLAGS) -o euplecton euplecton.cxx $(LIBS)

clean:
	rm -f euplecton
	rm -f *~









