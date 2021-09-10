# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../amrex

# libamrex location -- assumes it is $AMREX_HOME/dist => use ./configure
# --prefix=./dist
INCLUDE_LOCATIONS += $(AMREX_HOME)/dist/include
LIBRARY_LOCATIONS += $(AMREX_HOME)/dist/lib

# Localized compiler settings matching the one in AMReX
include $(AMREX_HOME)/Tools/GNUMake/Make.local


src := ../../..


include $(src)/src_hydro/Make.package
INCLUDE_LOCATIONS += $(src)/src_hydro/

include $(src)/src_gmres/Make.package
INCLUDE_LOCATIONS += $(src)/src_gmres/

include $(src)/src_common/Make.package
INCLUDE_LOCATIONS += $(src)/src_common/

include $(src)/src_rng/Make.package
INCLUDE_LOCATIONS += $(src)/src_rng/


CXXFLAGS := -fPIC -g -fopenmp
FCFLAGS  := -cpp
LDFLAGS  := -shared
OBJECTS  := $(CEXE_sources:.cpp=.o)
OBJECTS  += $(F90EXE_sources:.F90=.o)


LIBRARIES := $(addprefix -L, $(LIBRARY_LOCATIONS)) -lmpi -lamrex
INCLUDES  := $(addprefix -I, $(INCLUDE_LOCATIONS))


TARGET  = libhydro.so

.PHONY: all

all: $(TARGET)


%.o: $(src)/src_common/%.cpp
	$(CXX) $(CXXFLAGS) -c $(INCLUDES) -o $@ $<

%.o: $(src)/src_hydro/%.cpp common_namelist.o
	$(CXX) $(CXXFLAGS) -c $(INCLUDES) -o $@ $<

%.o: $(src)/src_gmres/%.cpp
	$(CXX) $(CXXFLAGS) -c $(INCLUDES) -o $@ $<

%.o: $(src)/src_rng/%.cpp
	$(CXX) $(CXXFLAGS) -c $(INCLUDES) -o $@ $<

# $(TARGET): $(OBJECTS)
# 	$(CXX) $(LIBRARIES) -o $(TARGET) -shared -Wl,--whole-archive $(AMREX_HOME)/dist/lib/libamrex.a -Wl,--no-whole-archive $(OBJECTS)
$(TARGET): $(OBJECTS)
	$(CXX) -fopenmp $(LIBRARIES) -o $(TARGET) $(OBJECTS)



clean:
	rm *.o
	rm *.mod
