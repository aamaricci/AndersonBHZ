#TO BE CHANGED BY USER:
#$ DRIVER NAME without .f90 extension
#$ COMPILER: supported compilers are ifort, gnu >v4.7 or use mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default is $HOME/.bin)

EXE=Abhz_2d
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

# LIBRARIES TO BE INCLUDED
#$ LIB_ED: either use *edlat* or *dmft_ed* (until we fix the naming conventions)
#$ LIB_SS: specify slave spins library if any
#$ LIB_TB: specify custom tight-binding library if any
#LIB_ED=edipack2
#LIB_SS=slave_spins
#LIB_TB=honeytools


#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################

ifdef LIB_ED
GLOB_INC+=$(shell pkg-config --cflags ${LIB_ED})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_ED})
endif

ifdef LIB_SS
GLOB_INC+=$(shell pkg-config --cflags ${LIB_SS})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_SS})
endif

ifdef LIB_TB
GLOB_INC+=$(shell pkg-config --cflags ${LIB_TB})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_TB})
endif

GLOB_INC+=$(shell pkg-config --cflags dmft_tools scifor)
GLOB_LIB+=$(shell pkg-config --libs   dmft_tools scifor)


ifeq ($(PLAT),intel)
FFLAG=-O3 -ftz	
DFLAG=-p -O0 -g
AFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPSERIAL =-fpp -D_
FPPMPI =-fpp -D_	
endif

ifeq ($(PLAT),gnu)
FFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
DFLAG = -w -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
AFLAG = -w -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime
FPPSERIAL= -cpp -D_
FPPMPI= -cpp -D_SCALAPACK
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif

OBJS = COMMON.o LCM_SQUARE.o


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

define colorecho	
	@tput setaf $2
	@tput bold
	@echo $1
	@tput sgr0
endef



all: FLAG:=${FFLAG} ${FPPSERIAL}
all: ${OBJS}
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) ${OBJS} $(FLAG)    $(EXE).f90 -o    $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(FC) ${OBJS} $(FLAG) mf_$(EXE).f90 -o $(DIREXE)/mf_$(EXE) ${GLOB_INC} ${GLOB_LIB}	
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)


mpi: FLAG:=${FFLAG} ${FPPMPI}
mpi: ${OBJS}
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) ${OBJS} $(FLAG)    $(EXE).f90 -o    $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(FC) ${OBJS} $(FLAG) mf_$(EXE).f90 -o $(DIREXE)/mf_$(EXE) ${GLOB_INC} ${GLOB_LIB}		
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)


debug: FLAG:=${DFLAG} ${FPPSERIAL}
debug: ${OBJS}
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) ${OBJS} $(FLAG)    $(EXE).f90 -o    $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(FC) ${OBJS} $(FLAG) mf_$(EXE).f90 -o $(DIREXE)/mf_$(EXE) ${GLOB_INC} ${GLOB_LIB}		
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)


aggressive: FLAG:=${AFLAG} ${FPPSERIAL}
aggressive: ${OBJS}
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) ${OBJS} $(FLAG)    $(EXE).f90 -o    $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(FC) ${OBJS} $(FLAG) mf_$(EXE).f90 -o $(DIREXE)/mf_$(EXE) ${GLOB_INC} ${GLOB_LIB}		
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)


clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/$(EXE) $(DIREXE)/mf_$(EXE)

.f90.o:	
	$(FC) $(FLAG) -c $< ${GLOB_INC}

version:
	@echo $(VER)


#########################################################################
