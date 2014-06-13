HEADERS    = vertex.h simplex.h sheet.h spacetime.h

DEBUG      = -g

OPT        = $(CXX_OPT)

CXX_FLAGS += -I$(SYNARMOSMA)/include -Wall -fPIC -DVERBOSE
LD_FLAGS  += -L$(SYNARMOSMA)/lib -Wall -shared

CXX_FLAGS += $(DEBUG) 
LD_FLAGS  += $(DEBUG)

LIBS       = -lsynarmosma $(LAPACK) -lboost_system -lpugixml -lntl -lm

