EXE=Abhz_2d
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

# LIBRARIES TO BE INCLUDED
LIB_ED=edipack2


#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################

ifdef LIB_ED
GLOB_INC+=$(shell pkg-config --cflags ${LIB_ED})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_ED})
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
FFLAG = -ffree-line-length-none -w  -fallow-argument-mismatch -O3   -funroll-loops
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
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



all: FLAG:=${FFLAG} ${FPPMPI}
all: ${OBJS}
all: compile


debug: FLAG:=${DFLAG} ${FPPMPI}
debug: ${OBJS}
debug: compile





compile: tb mf dmft data
	@echo ""
	@echo "Done"
	@echo ""

tb:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG)    $(EXE).f90 -o    $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)

mf:
	@echo ""
	$(call colorecho,"compiling mf_$(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG) mf_$(EXE).f90 -o $(DIREXE)/mf_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created mf_$(EXE) in  $(DIREXE)", 1)

dmft:
	@echo ""
	$(call colorecho,"compiling dmft_$(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG) dmft_$(EXE).f90 -o $(DIREXE)/dmft_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created dmft_$(EXE) in  $(DIREXE)", 1)

data:FLAG:=${FFLAG} ${FPPMPI}
data:
	@echo ""
	$(call colorecho,"compiling data_analysis_$(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG) get_data_$(EXE).f90 -o $(DIREXE)/get_data_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(FC) ${OBJS} $(FLAG) get_gf_$(EXE).f90 -o $(DIREXE)/get_gf_$(EXE) ${GLOB_INC} ${GLOB_LIB}	
	$(call colorecho,"created data_analysis_$(EXE) in  $(DIREXE)", 1)

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/*$(EXE)

.f90.o:	
	$(FC) $(FLAG) -c $< ${GLOB_INC}

version:
	@echo $(VER)


#########################################################################
