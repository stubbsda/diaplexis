HEADERS    = event.h simplex.h sheet.h spacetime.h

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

LIBS       = -lsynarmosma $(LAPACK) $(BOOST_TIMER) $(BOOST_SYSTEM) -lpugixml -lntl -lm
