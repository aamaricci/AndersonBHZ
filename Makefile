EXE=Abhz_2d
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
foo_ := $(shell pkg-config --exists scifor && echo yes || echo no)
ifneq ($(foo_),yes)
$(error ERROR: scifor is not installed or loaded)
endif
foo_ := $(shell pkg-config --exists dmft_tools && echo yes || echo no)
ifneq ($(foo_),yes)
$(error ERROR: dmft_tools is not installed or loaded)
endif
foo_ := $(shell pkg-config --exists edipack && echo yes || echo no)
ifneq ($(foo_),yes)
$(error ERROR: edipack is not installed or loaded)
endif

GLOB_INC+=$(shell pkg-config --cflags edipack dmft_tools scifor)
GLOB_LIB+=$(shell pkg-config --libs   edipack dmft_tools scifor)


ifeq ($(PLAT),intel)
FFLAG=-O3 -ftz
DFLAG=-p -O0 -g
AFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
endif

ifeq ($(PLAT),gnu)
FFLAG = -ffree-line-length-none -w  -fallow-argument-mismatch -O2
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
DFLAG = -w -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
AFLAG = -w -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif

OBJS = COMMON.o LCM_SQUARE.o

F2PY=f2py

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

define colorecho
	@tput setaf $2
	@tput bold
	@echo $1
	@tput sgr0
endef



all: FLAG:=${FFLAG}
all: ${OBJS}
all: compile


debug: FLAG:=${DFLAG}
debug: ${OBJS}
debug: compile





compile: tb mf dmft
	@echo ""
	@echo "Done"
	@echo ""

tb:FLAG:=${FFLAG}
tb: ${OBJS}
tb:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG)     $(EXE).f90 -o    $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)

mf:FLAG:=${FFLAG}
mf: ${OBJS}
mf:
	@echo ""
	$(call colorecho,"compiling mf_$(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG) mf_$(EXE).f90 -o $(DIREXE)/mf_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created mf_$(EXE) in  $(DIREXE)", 1)

dmft:FLAG:=${FFLAG}
dmft: ${OBJS}
dmft:
	@echo ""
	$(call colorecho,"compiling dmft_$(EXE).f90 ", 3)
	$(FC) ${OBJS} $(FLAG) dmft_$(EXE).f90 -o $(DIREXE)/dmft_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created dmft_$(EXE) in  $(DIREXE)", 1)





data:FLAG:=${FFLAG}
data: ${OBJS}
	@echo ""
	$(call colorecho,"compiling data_analysis codes", 3)
	$(FC) ${OBJS} $(FLAG) post-process/get_data_$(EXE).f90 -o $(DIREXE)/get_data_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(FC) ${OBJS} $(FLAG) post-process/get_gf_$(EXE).f90 -o $(DIREXE)/get_gf_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created data analysis codes in  $(DIREXE)", 1)

clean:
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
#	@rm -fv  $(DIREXE)/*$(EXE)

.f90.o:
	$(FC) $(FLAG) -c $< ${GLOB_INC}

version:
	@echo $(VER)


#########################################################################
