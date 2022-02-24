# Tunable parameters (currently set for Cray YMP)
#
# LDR           Program to load the objects into an executable
# LDFLAGS       For definition of lmset and csu on the Cray 2
#               Flags to the loader
# LIBS          A list of libraries to use on the Cray 2
#               Note that the C library (-lc) must ALWAYS be the last library
#               in the load step on both the XMP and the Cray 2
# COMMAND       Name of the executable
#
#F77 = ifort                 # name of f77 compiler
#C = icc
F77 = gfortran                 # name of f77 compiler
CC = gcc
FFLAGS = -c  -O2          # -O2 optimization
#FFLAGS =   -O2          # -O2 optimization
#FFLAGS = -c  -g          # -O2 optimization
CCFLAGS = -c  -O2          # -O2 optimization
#CCFLAGS = -c  -g          # -O2 optimization
LDR =	gfortran               # name of linker
LDRFLAGS=            -g 
LIBS = 	   -lm # 
#
# List of objects to build multidimensional VH-1 hydro code:

VHOBJS = driver.o stats.o flux.o
SOURCE = driver.f flux.f
OBJS = stats.o



# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build the executable
#
.f.o:
	$(F77) $(FFLAGS) $<

.c.o:
	$(CC) $(CCFLAGS) -DBL_FORT_USE_UNDERSCORE $<

ppmcns_new:	$(VHOBJS)
	$(LDR) $(LDRFLAGS) -o llns_1d.ex $(VHOBJS) $(LIBS)


#
#
clean:
	rm -f *.o *.l

clobber:	clean
	rm -f vhone


# Dependencies for the object files
stats.o:    stats.c    
driver.o:   driver.f
flux.o:   flux.f





