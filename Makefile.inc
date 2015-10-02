HEADERS    = event.h simplex.h sheet.h spacetime.h

DEBUG      = -g

OPT        = $(CXX_OPT)

CXX_FLAGS += -I$(SYNARMOSMA)/include -fPIC -DVERBOSE
LD_FLAGS  += -L$(SYNARMOSMA)/lib -shared

CXX_FLAGS += $(DEBUG) 
LD_FLAGS  += $(DEBUG)

LIBS       = -lsynarmosma $(LAPACK) $(BOOST) -lpugixml -lntl -lm

