###Set up some paths
TOP = $(PWD)
###################make file for DPMJET with pythia################
EXE = BeAGLE
MAIN = user-eA3.0-5
#MAIN  = $(TOP)main

###################compiler - options#######################
#32 bit
#CXX      = g77
#CXXFLAGS = -c -g -Wall -m32 -fno-inline -fno-automatic
#LD = $(CXX) -g -m32

#64 bit
CXX      = gfortran
#gcc 4.6.3 is used because of the fluka version.
#The 64 bit fluka lib can be compiled with gcc versions
#higher than 4.6, so our eic installed 64bit fluka
#is compiled with this version.
#To be compatible with that compiled fluka lib
#this dpmjet code has to be done with the same gcc version.
#When everything has been done, the gcc 4.6.3 lib has to
#be included in LD_LIBRARY_PATH. On eic, add
#/opt/gcc/4.6.3/lib64/ to to your LD_LIBRARY_PATH in 
#your .cshrc
#
#CXXFLAGS = -c -g -Wall -m64 -fno-inline -fno-automatic -O -fcheck=all -fbacktrace
CXXFLAGS = -c -g -Wall -m64 -fno-inline -fno-automatic
LD = $(CXX) -g -m64 

##################directories###############################
#PYTHIA   = PYTHIA-6.4.13
PYTHIA   = PYTHIA-6.4.28
PROGRAM  = dpmjet
PYQM = PyQM
#RADGEN = radgen
RADGEN = radgen-6.4.28

#to use the eic updated PYTHIA source
#uncomment the following two lines and comment out
#the pythia and radgen directory in current source list
#PYTHIASRC = /afs/rhic/eic/PACKAGES/PYTHIA-PP/PYTHIA-6.4.28/
#RADGENSRC = linkToRadgen/

#src, obj directory definition for
#pythia, pyqm and dpmjet radgen
PYTHIASRC = $(TOP)/$(PYTHIA)/
PYTHIAOBJ = $(TOP)/obj/$(PYTHIA)/
PYQMSRC = $(TOP)/$(PYQM)/
PYQMOBJ = $(TOP)/obj/$(PYQM)/
PROGRAMSRC = $(TOP)/src/
PROGRAMOBJ = $(TOP)/obj/$(PROGRAM)/
RADGENSRC = $(TOP)/$(RADGEN)/
RADGENOBJ = $(TOP)/obj/$(RADGEN)/

MAINOBJ = $(TOP)/obj/$(MAIN)
MAINSRC = $(TOP)/$(MAIN)

#####################libraries######################################
#change the following library path according to the setting in 
#your environment
#FLUKA = /afs/rhic/eic/PACKAGES/fluka-32/
#LIB1 = -L/cern/pro/lib -lmathlib -lkernlib -lpacklib_noshift -ldl -lm 
#LIB2 = -L$(FLUKA) -lflukahp
#LIB3 = -L/afs/rhic.bnl.gov/eic/lib32 -lLHAPDF 

# BNL:
#FLUKA = /afs/rhic/eic/PACKAGES/fluka-64/
#LIB1 = -L/cern64/pro/lib -lmathlib -lkernlib -lpacklib_noshift -ldl -lm 
#LIB2 = -L$(FLUKA) -lflukahp
#LIB3 = -L/afs/rhic/eic/lib -lLHAPDF 
#
# JLAB:
FLUKA = /u/group/ldgeom/PACKAGES/fluka-64/
LIB1 = -L/u/site/cernlib/x86_64_rhel6/2005/lib -lmathlib -lkernlib -lpacklib -ldl -lm 
LIB2 = -L$(FLUKA) -lflukahp
LIB3 = -L/u/group/ldgeom/PACKAGES/LHAPDF-5.9.1-64BIT/lib -lLHAPDF 

#####################include directories###########################
INCLUDE = $(TOP)/include
#change the following fluka path according to the setting in 
#your environment
FLUINC = $(FLUKA)flukapro

##This fpe.o module can not be used when PYTHIA runs
##since some pythia routines like pysppa which orders a
##initial shower will introduce some +-INF value but not
##allowed by this routine.
#TRAP = fpe.o

######################dependence rules#####################################
#pythia library
#get the stripped off file name list by patsubst
pythia_objects_ := $(patsubst %.f,%.o,$(notdir $(wildcard $(PYTHIASRC)*.f)))
#get the object file list which prefixed with the obj path directory
pythia_objects := $(patsubst %,$(PYTHIAOBJ)/%,$(pythia_objects_))

obj/libpythia6.a : $(pythia_objects)
	cd $(TOP)/obj; ar rcf libpythia6.a $(pythia_objects); ranlib libpythia6.a

#build the compiling rule for the source files
$(PYTHIAOBJ)/%.o : $(PYTHIASRC)%.f
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDE) $< -o $@

#radgen library
radgen_objects_ := $(patsubst %.f,%.o,$(notdir $(wildcard $(RADGENSRC)*.f)))
radgen_objects := $(patsubst %,$(RADGENOBJ)/%,$(radgen_objects_))

obj/libradgen.a : $(radgen_objects)
	cd $(TOP)/obj; ar rcf libradgen.a $(radgen_objects); ranlib libradgen.a

$(RADGENOBJ)/%.o : $(RADGENSRC)/%.f
	$(CXX)  $(CXXFLAGS) -I$(INCLUDE) $< -o $@

#pyqm library
pyqm_objects_ := $(patsubst %.f,%.o,$(notdir $(wildcard $(PYQMSRC)*.f)))
pyqm_objects := $(patsubst %,$(PYQMOBJ)/%,$(pyqm_objects_))

obj/libpyqm.a : $(pyqm_objects)
	cd $(TOP)/obj; ar rcf libpyqm.a $(pyqm_objects); ranlib libpyqm.a

$(PYQMOBJ)/%.o : $(PYQMSRC)%.f
	$(CXX)  $(CXXFLAGS) -I$(INCLUDE) $< -o $@

#dpmjet library
dpmjet_objects_ := $(patsubst %.f,%.o,$(notdir $(wildcard $(PROGRAMSRC)*.f)))
dpmjet_objects := $(patsubst %,$(PROGRAMOBJ)/%,$(dpmjet_objects_))

obj/libdpmjet.a : $(dpmjet_objects)
	cd $(TOP)/obj; ar rcf libdpmjet.a $(dpmjet_objects); ranlib libdpmjet.a

$(PROGRAMOBJ)/%.o : $(PROGRAMSRC)%.f
	$(CXX)  $(CXXFLAGS) -I$(INCLUDE) -I$(FLUINC) $< -o $@

#main program object
$(MAINOBJ).o: $(MAINSRC).f
	cd $(TOP); $(CXX) $(CXXFLAGS) -I$(INCLUDE) -I$(FLUINC) $(MAIN).f  -o $@

all: dir_maker	$(EXE)

#generate subdirectories in obj/
dir_maker:
	mkdir -p $(TOP)/obj $(PYTHIAOBJ) $(PYQMOBJ) $(PROGRAMOBJ) $(RADGENOBJ)

#The linking sequence must be dpmjet, pyqm, pythia6, radge
#Since there is a dependence relationship
#dpmjet needs routines from pythia6, pythia6 needs routines from radge
#The fundamental one always comes later!!!!
$(EXE): $(MAINOBJ).o obj/libdpmjet.a obj/libpythia6.a obj/libradgen.a obj/libpyqm.a
	cd $(TOP); $(LD) -I$(INCLUDE) -I$(FLUINC) -o $(EXE) $(MAINOBJ).o \
						-Lobj/ -ldpmjet -lpyqm -lpythia6 -lradgen \
							$(LIB1) $(LIB2) $(LIB3) 
#							-L/opt/gcc/4.6.3/lib64/
#						-Lobj/ -ldpmjet -L/afs/rhic/eic/PACKAGES/PYTHIA-RAD-CORR-32BIT/PYTHIA-6.4.28/ -lpythia6 -Lobj/ -lradgen -lpyqm \

#########################clean up utility##############################
clean:   
	rm -f $(EXE) core.* $(MAINOBJ).o $(PYTHIAOBJ)/*.o $(PROGRAMOBJ)/*.o $(PYQMOBJ)/*.o $(RADGENOBJ)/*.o obj/*.o obj/*.a

