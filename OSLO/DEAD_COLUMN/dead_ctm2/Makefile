# $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/Makefile,v 1.1 2003/04/15 14:41:35 alfgr Exp $

# Purpose: Makefile for aerosol module aer

# Usage:
# Debugging             : cd ~/aer;make OPTS=D aer;cd - 
# Double precision reals: cd ~/aer;make OPTS=D PRC=D aer;cd -
# Single precision reals: cd ~/aer;make OPTS=D PRC=S aer;cd -
# Turn-on sandblasting  : cd ~/aer;make USR_TKN=-DAlG01 aer;cd -

# Known problems: 

# Use separate directory for aer objects because they are often 
# compiled with -r8 -i4 and thus incompatible with, e.g., map
ifndef PVM_ARCH
 PVM_ARCH := $(shell pvmgetarch)
endif
MY_OBJ_DIR := ./obj/${PVM_ARCH}
MY_BLD_DIR := .
MY_DAT_DIR := .
MY_DOC_DIR := .
MY_DPN_DIR := ${MY_OBJ_DIR}
MY_INC_DIR := .
MY_SRC_DIR := .
ifndef MY_BIN_DIR
 MY_BIN_DIR := ./bin/${PVM_ARCH}
endif
ifndef MY_LIB_DIR
 MY_LIB_DIR := ../lib
endif

# Primary tokens which determine build options
# Specify non-default when invoking make, e.g. make DODS=Y
ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # endif DBG
ifndef DPN_SFX
 DPN_SFX := d
endif # endif DPN_SFX
ifndef MDL_SFX
 MDL_SFX := mod
endif # endif MDL_SFX
ifndef NETCDF_INC
 NETCDF_INC := /usr/local/include 
endif
ifndef NETCDF_LIB
 NETCDF_LIB := /usr/local/lib
endif
ifndef ${OMP} # OpenMP
 OMP := Y
endif # endif OMP
ifndef OPTS
 OPTS := X
endif
ifndef ${PRC}
# Precision token D=Double (default) S=Single
 PRC := S
endif # endif PRC
ifndef STC
# Created statically linked executable
 STC := N
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef USR_TKN
 USR_TKN := 
endif # endif USR_TKN
ifndef VERSION
 VERSION := $(shell cat ${MY_DOC_DIR}/VERSION)
endif

# Derived-tokens based on primary tokens
# These tokens should not be altered by hand
ifeq (${PRC},D)
 PRC_TKN := -DPRC_DBL
else
 PRC_TKN := -DPRC_FLT
endif # endif PRC

# NB: CPP_TKN is recursively expanded variable, define only when components are ready
CPP_TKN = ${USR_TKN} -DDST -DBXM -DCCM -D${PVM_ARCH} ${PRC_TKN}

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n

# Directories to search for source files
MDL_PTH := ${MY_SRC_DIR}
# Find all Fortran source files in directory
FIND_FNC = $(wildcard ${dir}/*.F ${dir}/*.F90)
# Assemble source files from all directories
SRC_LST = $(foreach dir, ${MDL_PTH},$(FIND_FNC))
# Source file names with directories removed
MDL_SRC := $(notdir $(SRC_LST))
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .${DPN_SFX}, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.${DPN_SFX} file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := aer # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All executables in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := # Include targets
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := ${MDL_INC_TRG} ${MDL_INC_SYM_LNK} # All header files in MY_INC_DIR
MDL_INC := $(addprefix ${MY_INC_DIR}/,${MDL_INC_STB}) # dst_cln removes these files

# C and C++ rules
${MY_OBJ_DIR}/%.o : %.c
	${CC} ${CPPFLAGS} ${CFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)
${MY_OBJ_DIR}/%.o : %.cc
	${C++} ${CPPFLAGS} ${C++FLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)

# Fortran rules
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
%.${MDL_SFX} : %.F90
	-ln -s -f ${NETCDF_INC}/netcdf.${MDL_SFX} .
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
%.f : %.F
	${CPP} -P ${CPPFLAGS} -o $@ $<

# Dependency rules
# Automatically generate a dependency file for each source file
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.${DPN_SFX}
DPN_GNR := makdep
${MY_DPN_DIR}/%.${DPN_SFX} : %.c
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.cc
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.F
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.F90
	@echo "Building dependency file $@"
	${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@

# The first LDFLAGS is for typical C programs with netCDF, math, and networking
# The second LDFLAGS enables C/Fortran linking

# Works on AIX and AIX46K
ifneq (${null},$(findstring AIX,${PVM_ARCH}))
ifeq (${OMP},Y)
 CC := xlc_r
 FC := xlf95_r
# -qsmp=noauto : Turn on SMP/OMP code generation but do no automatic parallelization
# -qsmp=omp : Use industry standard OMP without IBM extensions
 OMP_FLG := -qsmp=noauto
else
 CC := xlc
 CPP_DFN += -U_OPENMP
 FC := xlf95
endif # endif OMP
C++ := xlC
CFLAGS := -qmaxmem=-1
CPP := /lib/cpp -P
# AIX xlf requires wrapping all CPP tokens of type -DFOO and -UFOO in -WF flag
CPP_DFN := $(subst $(space)$(space),$(space),${CPP_TKN})
CPP_DFN_F := -WF,$(subst $(space),$(comma),$(CPP_DFN))
CPP_DFN_F := $(subst $(comma)$(comma),$(comma),$(CPP_DFN_F))
CPPINC := ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FFLAGS := ${OMP_FLG} -qmaxmem=-1 -qarch=auto -qspillsize=2500
FIXEDFLAGS := -qsuffix=f=f:cpp=F -qfixed=132
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
# MATCH adds to LDFLAGS: -bmaxdata:2147483648 -bmaxstack:268435456
LDFLAGS := ${OMP_FLG} -L${NETCDF_LIB} -lnetcdf -bnoquiet
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
 FFLAGS += -qrealsize=8 -qintsize=4 -qdpc=e
else
 FFLAGS += -qrealsize=4 -qintsize=4
endif
ifeq (${OPTS},O)
 CFLAGS +=  -O2
 FFLAGS += -O3 -qstrict -Q
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qsigtrap -qextchk
endif
ifeq (${OPTS},X)
 CFLAGS += -g
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qsigtrap -qextchk -C
endif
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # endif OMP
C++ := cxx -std strict_ansi
CC := cc
CFLAGS := ${OMP_FLG}
CPP := cpp
CPPFLAGS := -DOSF1 ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f95
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS := ${OMP_FLG} -L${NETCDF_LIB} -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif # PRC != D
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS += -O3 -inline speed
endif
ifeq (${OPTS},D)
 CFLAGS += -g -check_bounds -check -check_omp
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS := -g -check bounds -check omp_bindings -check overflow -check underflow
endif
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif ALPHA

ifeq (${PVM_ARCH},CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90
LD := ld
LDFLAGS := -L${NETCDF_LIB} -lnetcdf
LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O
# 1998/09/09: Added -F -O inline1 since BEE uses them on MATCH
 FFLAGS = -N 132 -F -O inline1
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS = -g -N 132 -F -O inline1
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS = -g -N 132 -e i -F -O inline1
endif
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
# 1999/05/24 Added -C for cpp to retain comments (else cpp thinks Fortran string concatenator // is a C++ comment and removes it)
${MY_OBJ_DIR}/%.o : %.F
	${CPP} -C -P ${CPPFLAGS} $< > $(patsubst %.F,%.f,$(notdir $<))
	${FC} -c ${FFLAGS} $(patsubst %.F,%.f,$(notdir $<)) 
#	${FC} -e p -c ${FFLAGS} $(notdir $<) 
#	${FC} -e p -Wp"${CPPFLAGS}" -c ${FFLAGS} $(notdir $<) 
	-mv -f $(notdir $@) ${MY_OBJ_DIR}
	rm -f $(patsubst %.F,%.f,$(notdir $<)) 
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} $<
	mv -f $(notdir $@) ${MY_OBJ_DIR}
endif
# endif CRAY

# Works on LINUX and LINUXALPHA
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
C++ := g++
CC := gcc
CPP := ${CC}
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
ifeq (${PVM_ARCH},LINUX)
 FC := lf95
# FC := pgf90
# FC := ifc
else # LINUXALPHA
 FC := g77
endif
LD := ld
LDFLAGS := -mp -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
# GNU Compiler Collection
ifeq (${CC},gcc)
 CFLAGS := -Wall
ifeq (${OPTS},O)
 CFLAGS += -O
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS +=
endif
ifeq (${OPTS},X)
 CFLAGS += -g -O
endif
 C++FLAGS := ${CFLAGS}
endif
# Intel (Kai) C Compiler
ifeq (${CC},icc)
 CFLAGS := 
 OMP_FLG_C := -openmp
ifeq (${OPTS},O)
 CFLAGS += -O
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS +=
endif
ifeq (${OPTS},X)
 CFLAGS += -g
endif
 C++FLAGS := ${CFLAGS}
endif
# Intel (Kai) Fortran Compiler
ifeq (${FC},ifc)
# ifc produces "dictionary" *.d file
# ifc produces "program unit catalog" *.pc file
# -fpp2 necessary, but not sufficient, for  OpenMP
# -w95 supresses warnings for non-standard fortran (which includes OpenMP)
 FFLAGS := -extend_source -implicitnone -vms -fpp2 -w95
# -lVaxlib needed for iargc_, getarg_
 LDFLAGS += -lPEPCF90
 MDL_SFX := d
 OMP_FLG_F := -openmp
ifeq (${PRC},D)
 FFLAGS += -i4 -r8 -doubletemps
else
 FFLAGS += -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -O2
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
# -e95 issues warnings for non-standard fortran (which includes OpenMP)
 FFLAGS += -g -C -e95
endif
endif # end Intel (Kai) Fortran Compiler
# Lahey-Fujitsu Fortran Compiler
ifeq (${FC},lf95)
 FFLAGS := -DLinux -DHIDE_SHR_MSG --mod .:${NETCDF_INC} --wide
 OMP_FLG_F := --openmp
ifeq (${PRC},D)
 FFLAGS += --dbl
else
 FFLAGS += 
endif
ifeq ($(STC),Y)
 LDFLAGS := --staticlink ${LDFLAGS}
endif # endif STC
ifeq (${OPTS},O)
# --tpp: Pentium Pro/Pentium II code generation
# --[n]sav: [Do not] save local variables
 FFLAGS += -O3 --tpp --nchk --ntrace --nsav --nwarn
endif
ifeq (${OPTS},D)
 FFLAGS += -g --chk a,e,s --pca --trace --trap
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
 FFLAGS += -g --chk a,e,s,u,x --chkglobal --pca --trace --trap
endif
endif # end Lahey-Fujitsu Fortran Compiler
# Portland Group Fortran Compiler
ifeq (${FC},pgf90)
# 20020403: Omit -ktrap=fp so MATCH works in single precision
# FFLAGS := -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign -ktrap=fp
 FFLAGS := -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign
 OMP_FLG_F := -mp
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -fast
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
 FFLAGS += -g -Mbounds
endif
endif # end Portland Group Fortran Compiler
# G77 Fortran compiler
ifeq (${FC},g77)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore
ifeq (${OPTS},O)
 FFLAGS += -O
endif
ifeq (${OPTS},D)
 FFLAGS += -g -fdebug-kludge
endif
ifeq (${OPTS},R)
 FFLAGS += -fdebug-kludge
endif
ifeq (${OPTS},X)
 FFLAGS += -g -O -fdebug-kludge -fbounds-check
endif
endif # endif G77 Fortran compiler
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_F}
 LDFLAGS := ${OMP_FLG_F} ${LDFLAGS}
endif # endif OMP
endif
# endif LINUX

ifeq (${UNAMES},IRIX64)
C++ := CC
CC := cc
CFLAGS := -64 -mips4
CPP := /usr/lib/cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90 -cpp
FFLAGS := -64 -mips4 -extend_source -mp -mpio
LD := ld
LDFLAGS := -64 -mp -mips4 -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${PRC},D)
 FFLAGS += -r8 -i4 
else
 FFLAGS += -r4 -i4
endif
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS += -O2
endif
ifeq (${OPTS},R)
 CFLAGS += 
 FFLAGS += 
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g -trapuv
 FFLAGS += -g -check_bounds -ansi -uvar
endif
endif
# endif IRIX64=SGI64,SGIMP64

ifeq (${UNAMES},SunOS) 
C++ := g++
CC := gcc -ansi
CFLAGS := -Wall
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90
FFLAGS := -xs -stackvar -e
LD := ld
LDFLAGS := -L/opt/SUNWspro/lib -lf77compat -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS += -fast -inline=%auto
 LDFLAGS += -fast # 20010420 Yes, -fast tells Sun linker to add extra libraries
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g 
 FFLAGS += -g
ifeq (${FC},f77)
# NB: 19980601 -C (range-checking) is not supported by Sun f90
 FFLAGS += -C
endif
endif
ifeq (${OMP},Y)
 FFLAGS += -openmp # Equivalent to -mp=openmp -explicitpar -stackvar -D_OPENMP
 LDFLAGS += -openmp # 
endif
ifeq (${PRC},D)
# Do not add -r8const flag unless all constants are given precision in the code
# FFLAGS += -r8const
# Using xtypemap appears to cause aer to fail
 FFLAGS += -xtypemap=real:64,double:64,integer:32
ifeq (${FC},f77)
 FFLAGS += -r8 -i4
endif
else
 FFLAGS += -xtypemap=real:32,double:64,integer:32
ifeq (${FC},f77)
 FFLAGS += -r4 -i4
endif
endif
# Sun compilers seem unable to find netcdf.${MDL_SFX}
${MY_OBJ_DIR}/%.o : %.F90
	-ln -s -f ${NETCDF_INC}/netcdf.${MDL_SFX} .
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif SunOS=SUN4SOL2,SUNMP

ifeq (${DBG},Y)
# Enable Debugging token 
 CPPFLAGS += -DDST_DBG
endif # endif DBG

# Disable OpenMP on platforms that automatically support it
ifeq (${OMP},N)
 FFLAGS := $(filter-out -omp,${FFLAGS})
 FIXEDFLAGS := $(filter-out -omp,${FIXEDFLAGS})
 CFLAGS := $(filter-out -mpio,${CFLAGS})
 FFLAGS := $(filter-out -mpio,${FFLAGS})
 LDFLAGS := $(filter-out -mpio,${LDFLAGS})
 CFLAGS := $(filter-out -mp,${CFLAGS})
 FFLAGS := $(filter-out -mp,${FFLAGS})
 LDFLAGS := $(filter-out -mp,${LDFLAGS})
 ifneq (${UNAMES},SunOS)
  CPPFLAGS += -U_OPENMP
 endif # endif SunOS
endif # endif OMP

ifeq ($(STC),Y)
# Create statically linked executable
 LDFLAGS := -Wl,-Bstatic ${LDFLAGS}
endif # endif STC

# Define any remaining variables

# Default targets
all: dir ${MDL_BIN_TRG}
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .${DPN_SFX} .F .F90 .f90 .${MDL_SFX} .o
# Delete targets which were not successfully made
.DELETE_ON_ERROR:
# Target directories which may not exist
dir: bin_dir obj_dir 
bin_dir:
	-mkdir -p ${MY_BIN_DIR}
obj_dir:
	-mkdir -p ${MY_OBJ_DIR}

# Targets in bin
aer : ${MDL_OBJ}
	${FC} -o ${MY_BIN_DIR}/$@ ${MDL_OBJ} ${LDFLAGS}
	chmod 755 ${MY_BIN_DIR}/$@

MAKDEP_OBJ :=	${MY_OBJ_DIR}/makdep.o
# Use system getopt() rather than BSD getopt() to eliminate unnecessary dependencies
CPPFLAGS:= -I/usr/include ${CPPFLAGS}
makdep:	${MY_BIN_DIR}/makdep
${MY_BIN_DIR}/makdep:	$(MAKDEP_OBJ) 
	${CC} -o $@ $(MAKDEP_OBJ)
	chmod 755 $@

bin : ${MDL_BIN_TRG}
	rm -f ${MDL_BIN}
binclean: bin_cln
bin_cln: 
	rm -f ${MDL_BIN}
# Targets in bld
buildclean: bld_cln
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
# Targets in dat
data:	dat
dat_cln :
	-cd ${MY_DAT_DIR}; rm -f aer.nc dst_mss_bdg.nc
# Targets in dpn
depend: dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in inc
include: inc
inc : 
	cp -f -p ${MDL_INC_TRG} ${MY_INC_DIR}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
library: lib
lib :
	-cd ${HOME}/f; gmake lib
libclean: lib_cln
lib_cln :
	-cd ${MY_LIB_DIR};
# Targets in mod
mod : ${MDL_OBJ}
	mv -f $(wildcard ${MY_OBJ_DIR}/*.mod) .
mod_cln :
	rm -f $(wildcard *.mod) $(wildcard ${MY_OBJ_DIR}/*.mod)
# Targets in obj
object: obj
obj : ${MDL_OBJ}
objclean: obj_cln
obj_cln :
	rm -f ${MDL_OBJ}

# Housekeeping
clean : cln
cln : dat_cln dpn_cln lib_cln mod_cln obj_cln
debug: dbg
dbg:
	@printf "C++ = ${C++}\n"
	@printf "C++FLAGS = ${C++FLAGS}\n"
	@printf "CC = ${CC}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CPP_TKN = ${CPP_TKN}\n"
	@printf "DBG = ${DBG}\n"
	@printf "FC = ${FC}\n"
	@printf "FFLAGS = ${FFLAGS}\n"
	@printf "FIXEDFLAGS = ${FIXEDFLAGS}\n"
	@printf "FREEFLAGS = ${FREEFLAGS}\n"
	@printf "GSL_LIB = ${GSL_LIB}\n"
	@printf "HOST = ${HOST}\n"
	@printf "HOSTNAME = ${HOSTNAME}\n"
	@printf "LDFLAGS = ${LDFLAGS}\n"
	@printf "MDL_DPN = ${MDL_DPN}\n"
	@printf "MDL_OBJ = ${MDL_OBJ}\n"
	@printf "MDL_PTH = ${MDL_PTH}\n"
	@printf "MDL_SRC = ${MDL_SRC}\n"
	@printf "MY_BIN_DIR = ${MY_BIN_DIR}\n"
	@printf "MY_BLD_DIR = ${MY_BLD_DIR}\n"
	@printf "MY_DAT_DIR = ${MY_DAT_DIR}\n"
	@printf "MY_DOC_DIR = ${MY_DOC_DIR}\n"
	@printf "MY_DPN_DIR_RX = ${MY_DPN_DIR_RX}\n"
	@printf "MY_ES_DIR = ${MY_ES_DIR}\n"
	@printf "MY_INC_DIR = ${MY_INC_DIR}\n"
	@printf "MY_LIB_DIR = ${MY_LIB_DIR}\n"
	@printf "MY_SHR_DIR = ${MY_SHR_DIR}\n"
	@printf "MY_OBJ_DIR = ${MY_OBJ_DIR}\n"
	@printf "MY_OBJ_DIR_RX = ${MY_OBJ_DIR_RX}\n"
	@printf "NETCDF_INC = ${NETCDF_INC}\n"
	@printf "NETCDF_LIB = ${NETCDF_LIB}\n"
	@printf "OMP = ${OMP}\n"
	@printf "OPTS = ${OPTS}\n"
	@printf "PRC = ${PRC}\n"
	@printf "PVM_ARCH = ${PVM_ARCH}\n"
	@printf "UNAMES = ${UNAMES}\n"
	@printf "VPATH = ${VPATH}\n"
distclean : dst_cln
dst_cln : bld_cln cln lib_cln
DST_MFL := ChangeLog FAQ index.shtml INSTALL MANIFEST Makefile README TODO VERSION # Distribution metafiles
DOT_H_FL := dst.h misc.h params.h pmgrid.h preproc.h # Distribution .h files
dst : ${DST_MFL} ${SRC_LST}
	tar cvzf ${DATA}/tmp/dead-${VERSION}.tar.gz ${DST_MFL} ${DOT_H_FL} ${SRC_LST}
	scp -p ${DST_MFL} ${HOME}/c/makdep.c ${DATA}/tmp/dead-${VERSION}.tar.gz dust.ess.uci.edu:/var/www/html/dead
	-ssh dust.ess.uci.edu "cd /var/www/html/dead;ln -s -f dead-${VERSION}.tar.gz dead.tar.gz"
	scp -p ${DATA}/tmp/dead-${VERSION}.tar.gz dust.ess.uci.edu:/var/ftp/pub/zender/dead
	-ssh dust.ess.uci.edu "cd /var/ftp/pub/zender/dead;ln -s -f dead-${VERSION}.tar.gz dead.tar.gz"
pch : 
	-mkdir -p ${MY_SRC_DIR}/glue/match/readers;cp -f ${HOME}/match_dst/readers/mss.c ${MY_SRC_DIR}/glue/match/readers
	-mkdir -p ${MY_SRC_DIR}/glue/match/src;cp -f ${HOME}/match_dst/src/main.F ${HOME}/match_dst/src/optDepth.F ${HOME}/match_dst/src/physlic.F ${HOME}/match_dst/src/histout.F ${MY_SRC_DIR}/glue/match/src
	-mkdir -p ${MY_SRC_DIR}/glue/ccm3/src/physics;cp -f ${HOME}/ccm_dst/src/physics/aphys.F ${MY_SRC_DIR}/glue/ccm3/src/physics
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.F90 *.f90 *.hh *.h Makefile
uninstall : bin_cln dst_cln inc_cln lib_cln

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean dir distclean dst_cln dpn_cln tags uninstall
ifeq (${null},$(findstring $(MAKECMDGOALS),${GOALS_WHICH_DELETE_DEPENDENCY_FILES}))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq (${INCLUDE_DPN},TRUE)
# Following incorporates dependency files into Makefile rules
-include ${MDL_DPN}
endif
