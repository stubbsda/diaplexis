HEADERS    = event.h simplex.h sheet.h spacetime.h

CXX_FLAGS += -I$(SYNARMOSMA)/include $(OPENMP) -fPIC
LD_FLAGS  += -L$(SYNARMOSMA)/lib $(OPENMP) -shared

#CXX_FLAGS += $(CXX_OPT)
CXX_FLAGS += $(DEBUG) 

#LD_FLAGS += $(CXX_OPT)
LD_FLAGS  += $(DEBUG)

LIBS       = -lsynarmosma $(LAPACK) $(BOOST) -lpugixml -lntl -lm

