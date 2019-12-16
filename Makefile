##################################################################
# Oslo CTM3
# Makefile for at least one platform.
##################################################################
#
#     Does not use a dependency generator, but lists specific
#     dependencies at the bottom of the file.
#
#     Tested and works on UiO linux cluster Abel.
#
#     Ole Amund Sovde, September 2008 - November 2015
##################################################################

#Get machine info
UNAMES   = $(shell uname -s)
UNAMEM   = $(shell uname -m)
UNAMEN   = $(shell uname -n)
MYNAME   = $(shell whoami)
AMDPROC  = $(shell grep -i AMD /proc/cpuinfo)
INTELPROC= $(shell grep -i Intel /proc/cpuinfo)
#Set file names
MAKEFILE = Makefile
MAIN     = osloctm3
PROGNAME = Oslo CTM3
COMPDIR  = $(shell pwd)

##################################################################
#--- SET USER CHOICES (NO SPACES BEFORE OR AFTER Y/N SIGNS) ------
#-----------------------------------------------------------------
#Compiler options O=Optimize, D=Debug, A=Agressive optimization
OPTS :=A
#Native (metdata) horizontal resolution (HNATIVET42 / HNATIVET159 / ...)
HNATIVE :=HNATIVET159
#Native (metdata) vertical resolution (L40 / L60)
VNATIVE :=VNATIVEL60
#
#CTM run resolution (degrade horizontal, collapse vertical)
#Horizontal resolution (HORIGINAL/HTWO/HFOUR)
HWINDOW :=HTWO
#Collapse layer 1:3 and 4:5?
COLLAPSE :=N
#-----------------------------------------------------------------
#Compile with Oslo chemistry/physics
OSLOCHEM :=Y
#Compile with Oslo tropospheric chemistry
TROPCHEM :=Y
#Compile with Oslo stratospheric chemistry
STRATCHEM :=Y
#Sulfur scheme
SULPHUR :=Y
#BC/OC
BCOC :=N
#Nitrate (requires SALT and SULPHUR)
NITRATE :=Y
#Sea salt
SEASALT :=Y
#Dust
DUST :=N
#Secondary organic aerosols
SOA :=N
#E90 (MUST be Y to calculate STE fluxes!)
E90 :=Y
#LINOZ (NOT set up to replace stratospheric chemistry, only included
#       to calculate STE flux as UCI does it.
LINOZ :=N
#-----------------------------------------------------------------
#M7 (not implemented)
M7 :=N
#-----------------------------------------------------------------
#Lightning scaling factors: To produce lightning scaling factors,
#turn off all the above (except OC) and set LIT:=Y (sets NPAR=1).
#You should also remove all transport calls from pmain.f90, e.g.
#in a file pmain_lit.f90, and update Makefile accordingly.
#Also set getFlashFiles = .true. in lightning.f90.
LIT :=N
#-----------------------------------------------------------------
#Emission & dry deposition treatment (U=UCI, O=OSLO)
#UCI: Separate processes, OSLO: Production and loss terms in chemistry
#Note that SALT and DUST use OSLO approach anyway.
EMISDEP_TREATMENT :=O
#-----------------------------------------------------------------
#Fortran compiler (ifort, pgf90)
FC=ifort
#FC=pgf90
#-----------------------------------------------------------------
#Object files are stored in their source directories
#--- END USER CHOICES --------------------------------------------
##################################################################


##################################################################
#Normally you don't need to touch anything below this line.
##################################################################


#-----------------------------------------------------------------
#--- INITIALIZATIONS ---------------------------------------------
#-----------------------------------------------------------------
null=
space=$(null) $(null)
kolon=$(null):$(null)
#Directory where object files are put during compilation
MY_OBJ_DIR:=bin
#Directory where .mod files can be put during compilation
MY_MOD_DIR:=mod
#Initialize includes
INCLUDES :=-I./
#Initialize LDFLAGS (Libraries needed when linking object files)
LDFLAGS :=
#-----------------------------------------------------------------


#-----------------------------------------------------------------
#--- COMPILER OPTIONS --------------------------------------------
#-----------------------------------------------------------------
# Sets the right compiler and compiler options for the right
# architecture. In general these are machine dependent things, but
# you may want to change the name of the compiler for the given
# architecture.
#-----------------------------------------------------------------

#--- LINUX MACHINES ----------------------------------------------
ifeq ($(UNAMES),Linux)
# First specific options for 64/32 bits
ifneq ($(null),$(findstring _64,$(UNAMEM)))
# === 64 bit === x86_64 ===
ifneq ($(null),$(findstring uci.edu,$(UNAMEN)))
#UCI
NETCDF_LIB=netCDF/lib
NETCDF_LIB=netCDF/include
NETCDF_LLIBS=-lnetcdf
endif #UCI netCDF
ifneq ($(null),$(findstring .local,$(UNAMEN)))
#UiO Abel nodes
NETCDF_INC=$(NETCDF_ROOT)/include
NETCDF_LIB=$(NETCDF_ROOT)/lib
NETCDF_LLIBS=-lnetcdff -lnetcdf
endif
ifneq ($(null),$(findstring login-,$(UNAMEN)))
#UiO SAGA nodes
NETCDF_INC=$(NETCDF_ROOT)/include
NETCDF_LIB=$(NETCDF_ROOT)/lib
NETCDF_LLIBS=-lnetcdff -lnetcdf
endif
ifneq ($(null),$(findstring .hpc.uio.,$(UNAMEN)))
#UiO enso/nao
NETCDF_INC=$(NETCDF_ROOT)/include
NETCDF_LIB=$(NETCDF_ROOT)/lib
NETCDF_LLIBS=-lnetcdff -lnetcdf
endif
endif
# === END 64 bit === x86_64 ===
# === 32 bit not available ===


# *** IFORT ******************************************************
ifeq ($(FC),ifort)
CC=cc
FIXEDFLAGS = -fixed
FREEFLAGS  = -free
#
# intel version 10-15
#STDOPTS= -auto -fpp -mcmodel=large -shared-intel -mp1 -module $(MY_MOD_DIR)
#FCOPT  = -inline-forceinline -ip -O2 -openmp -ftz -fno-alias -fomit-frame-pointer
#FCAGR  = -inline-forceinline -ip -O3 -openmp -ftz -fno-alias -fomit-frame-pointer
#DEBUG  = -check bounds -check uninit -traceback -openmp -g
#
# intel version >= 17
STDOPTS= -auto -fpp -mcmodel=large -shared-intel -mp1 -module $(MY_MOD_DIR)
FCOPT  = -inline-forceinline -ip -O2 -qopenmp -ftz -fno-alias -fomit-frame-pointer
FCAGR  = -inline-forceinline -ip -O3 -qopenmp -ftz -fno-alias -fomit-frame-pointer
DEBUG  = -check bounds -check uninit -traceback -qopenmp -g
#
#Special x86_64 options:
ifneq ($(null),$(findstring _64,$(UNAMEM)))
#intel processor?
ifneq ($(null),$(INTELPROC))
FCOPT +=
FCAGR +=
endif  #check on intel processor
#AMD processor?
ifneq ($(null),$(AMDPROC))
FCOPT+=
FCAGR+=
endif #check on AMD processor
endif #check on Special x86_64 options
endif #check on ifort compiler
# *** END IFORT **************************************************


# *** PORTLAND (PGI) *********************************************
ifeq ($(FC),pgf90)
#INCLUDES += -I/usr/lib64
#LDFLAGS += -L/usr/lib64 -lnuma
CC=pgcc
FIXEDFLAGS = -Mfixed
FREEFLAGS = -Mfree
STDOPTS = -mcmodel=medium -Mdepchk 
FCOPT  = -O2 -Mprefetch -Minline=reshape
FCAGR  = -O3 -mp -Mprefetch -Minline -Munroll
#         -Mvect=sse does not work for some loops
DEBUG  = -Mbounds -g
#Check if you are on titan==> some more options apply
ifneq ($(null),$(findstring _64,$(UNAMEM)))
FCOPT +=
FCAGR += 
endif  #check on titan
endif  #check on using pgf90 compiler
# *** END PORTLAND (PGI) *****************************************


endif  #check on linux
#--- END LINUX MACHINES ------------------------------------------

#--- COMPILER FLAGS ----------------------------------------------
# Use Optimization or Debug 
ifeq ($(OPTS),A)
OPTFLAGS=$(FCAGR)
endif
ifeq ($(OPTS),O)
OPTFLAGS=$(FCOPT)
endif
ifeq ($(OPTS),D)
OPTFLAGS=$(DEBUG)
endif

#INCLUDE LIBRARIES (e.g. netCDF)
INCLUDES += -I$(NETCDF_INC)
LDFLAGS += -L$(NETCDF_LIB) $(NETCDF_LLIBS)


# Token for collapsing layers
ifeq ($(COLLAPSE),Y)
DEGRADE_TKN := -DCOLLAPSE
endif
#--- END COMPILER FLAGS ------------------------------------------

#-----------------------------------------------------------------
#--- END COMPILER OPTIONS ----------------------------------------
#-----------------------------------------------------------------


##################################################################


#-----------------------------------------------------------------
#--- SOURCE FILES ------------------------------------------------
#-----------------------------------------------------------------
# Initialize list of all directories in use
ALL_DIRS:=
#------------------ CORE -----------------------------------------
# List of source files. Some have changed from UCI to CTM3 (_oc-files)
CORE_SRC := \
	cmn_precision.f90 \
	cmn_size.F90 \
	cmn_ctm.f90 \
	cmn_met.f90 \
	cmn_chem.f90 \
	cmn_diag.f90 \
	cmn_fjx.f90 \
	cmn_sfc.f90 \
	cmn_parameters.f90 \
	utilities.f90 \
	averages.f90 \
	budgets.f90 \
	cloudjx.f90 \
	p-cloud2.f \
	convection.f90 \
	p-dyn0.f \
	p-dyn2.f \
	fastjx.f90 \
	initialize.f90 \
	regridding.f90 \
	grid.f90 \
	p-linoz.f \
	lightning.f90 \
	pmain.f90 \
	omp.f90 \
	pbl_mixing.f90 \
	p-phot_oc.f \
	scavenging_drydep_uci.f90 \
	scavenging_largescale_uci.f90 \
	p-series.f \
	source_uci.f90 \
	p-vect3.f \
	steflux.f90 \
	spectral_routines.f \
	stt_save_load.f90 \
	metdata_ecmwf.f90

#------------------ END CORE -------------------------------------


#------------------ OSLO CHEMISTRY -------------------------------
# The use of DUMMY routines is essential to Oslo CTM3. Learn how
# to use them.
#-----------------------------------------------------------------
ifeq ($(OSLOCHEM),Y)
# Initialize source code list
OSLO_SRC :=
# Oslo chemistry source path
OSLO_PATH:= OSLO
# Add to ALL_DIRS
ALL_DIRS += $(OSLO_PATH)
# Search for includes in this directory
INCLUDES += -I$(OSLO_PATH)
# Token for compiling params.h
OSLO_TKN := -DOSLOCHEM

# Path for DUMMIES
OSLODUM_PATH := DUMMIES
# Add to ALL_DIRS
ALL_DIRS += $(OSLO_PATH)/$(OSLODUM_PATH)

# CTM3 core -----------------------------------------
# First list utilities used by other routines.
# List of source files. The order is arbitrary, just remember to
# specify dependencies at the bottom of the Makefile.
OSLO_SRC := \
	cmn_oslo.f90 \
	cnv_oslo.f90 \
	dateconv.f90 \
	atom.f90 \
	caribic2.f90 \
	diagnostics_general.f90 \
	diagnostics_scavenging.f90 \
	drydeposition_oslo.f90 \
	gmdump3hrs.f90 \
	emissions_aircraft.f90 \
	emissions_megan.f90 \
	emissions_ocean.f90 \
	emissions_oslo.f90 \
	emissions_volcanoes.f90 \
	emisutils_oslo.f90 \
	fallingaerosols.f90 \
	input_oslo.f90 \
	main_oslo.f90 \
	ncutils.f90 \
	physics_oslo.f90 \
	qssa_integrator.f90 \
	ch4routines.f90 \
	satelliteprofiles_mls.f90 \
	seasaltprod.f90 \
	strat_aerosols.f90 \
	strat_h2o.f90 \
	strat_loss.f90 \
	troccinox_ban.f90 \
	troccinox_fal.f90 \
	troccinox_geo.f90 \
	utilities_oslo.f90 \
	verticalprofiles_stations2.f90 \
	aerosols2fastjx.f90 \
	hippo.f90

# Add packages / dummies ----------------------------

# Sulphur -------------------------------------------
ifeq ($(SULPHUR),Y)
OSLO_SRC += sulphur_oslo.f90
OSLO_TKN += -DSULPHUR
else
OSLO_SRC += \
	$(OSLODUM_PATH)/sulphur_oslo.f90
endif


# Emission & drydeposition treatment ---------------
# Check only for O; assume U if not.
ifeq ($(EMISDEP_TREATMENT),O)
OSLO_TKN += -DEMISDEPCHEM
OSLO_SRC += emisdep4chem_oslo.f90
else
OSLO_SRC += $(OSLODUM_PATH)/emisdep4chem_oslo.f90
endif


# Tropospheric chemistry ----------------------------
ifeq ($(TROPCHEM),Y)
OSLO_SRC += \
	tropchem_oslo.f90 \
	chem_oslo_rates.f90 \
	pchemc_ij.f90
OSLO_TKN += -DTROPCHEM
else
OSLO_SRC += \
	$(OSLODUM_PATH)/chem_oslo_rates.f90 \
	$(OSLODUM_PATH)/tropchem_oslo.f90
endif

# Stratospheric chemistry ---------------------------
ifeq ($(STRATCHEM),Y)
OSLO_SRC += \
	psc_microphysics.f90 \
	stratchem_oslo.f90 \
	$(OSLODUM_PATH)/strat_o3noy_clim.f90 \
	pchemc_str_ij.f90
OSLO_TKN += -DSTRATCHEM
else
OSLO_SRC += \
	$(OSLODUM_PATH)/stratchem_oslo.f90 \
	strat_o3noy_clim.f90 \
	$(OSLODUM_PATH)/psc_microphysics.f90
endif


# M7 ------------------------------------------------
ifeq ($(M7),Y)
M7_TKN := -DM7
endif


# BCOC ----------------------------------------------
ifeq ($(BCOC),Y)
OSLO_SRC += bcoc_oslo.f90
OSLO_TKN += -DBCOC
else
OSLO_SRC += $(OSLODUM_PATH)/bcoc_oslo.f90
endif

# SOA ----------------------------------------------
ifeq ($(SOA),Y)
OSLO_SRC += soa_oslo.f90
OSLO_TKN += -DSOA
else
OSLO_SRC += $(OSLODUM_PATH)/soa_oslo.f90
endif


# SEA SALT ------------------------------------------
ifeq ($(SEASALT),Y)
OSLO_SRC += seasalt.f90
OSLO_TKN += -DSEASALT
else
OSLO_SRC += $(OSLODUM_PATH)/seasalt.f90
endif


# Nitrate -------------------------------------------
ifeq ($(NITRATE),Y)
OSLO_SRC += \
	eqsam_v03d.f90 \
	nitrate.f90
OSLO_TKN += -DNITRATE
else
OSLO_SRC += \
	$(OSLODUM_PATH)/nitrate.f90
endif


# DUST ----------------------------------------------
# Number of DUST species: default=8, M7=2
ifeq ($(DUST),Y)
OSLO_SRC += dust_oslo.f90
OSLO_TKN += -DDUST
# Alfaro & Gomes Sandblasting: add -DAlG01 to OSLO_TKN
# Directory for DUST/DEAD
DUST_PATH := $(OSLO_PATH)/DEAD_COLUMN
# Includes for separate DUST compilation; may add -DDST_DBG
DUST_INC := -I$(DUST_PATH) -I$(NETCDF_INC) -I$(NETCDF_LIB) -D$(HNATIVE) \
	-D$(HWINDOW) -D$(VNATIVE) -DDST $(DEGRADE_TKN) $(M7_TKN)
# DUST source file list
DUST_LIST = \
	dead_precision.F90 pmgrid.F dstgrd.F90 dstdbg.F90 \
	dbg_mdl.F90 sng_mdl.F90 xtr_mdl.F90 utl_mdl.F90 \
	dstctl.F90 dstcst.F90 vec_mdl.F90 dpsdryutl.F90 \
	erf_mdl.F90 psdlgn.F90 dstaer.F90 nf90_utl.F90 \
	dstpsd.F90 dstodx.F90 dstchm.F90 dstbdg.F90 \
	dstmssutl.F90 dstblm.F90 dstnm.F90 dead_history.f90 \
	dstlsm.F90 blmutl.F90 dstscv.F90 dstsltsbl.F90 \
	dstcmnini.F90 dsttvbds.F90 dstsfc.F90 dsttibds.F90 \
	dstdpsdry.F90 dstmblutl.F90 gmm_mdl.F90 wbl_mdl.F90 \
	dstmbl.F90 dead_inirun.f90
#	 aer.F90 qneg.F aernvr.F90 dstdpswet.F90 dstrad.F90 phyzlic.F90
# Expand source list with DUST_PATH
DUST_SRC := $(addprefix $(DUST_PATH)/,$(DUST_LIST))
# Dust objects
DUST_OBJ= $(addprefix $(MY_OBJ_DIR)/,$(addsuffix .o,$(basename $(DUST_SRC))))
# Dust paths
DUSTOBJDIRS:= $(addprefix $(MY_OBJ_DIR)/, $(DUST_PATH))
# Add to ALL_DIRS
ALL_DIRS += $(DUST_PATH)
else
# Dust dummy
OSLO_SRC += $(OSLODUM_PATH)/dust_oslo.f90
endif


# E90 ------------------------------------------------------------
ifeq ($(E90),Y)
OSLO_TKN += -DE90
endif
# LINOZ ----------------------------------------------------------
ifeq ($(LINOZ),Y)
OSLO_TKN += -DLINOZ
endif
# LIGTHNING generation / read-through ----------------------------
ifeq ($(LIT),Y)
OSLO_TKN += -DLITGEN
endif


endif # ifeq ($(OSLOCHEM),Y)
#------------------ END OSLO CHEMISTRY ---------------------------



#------------------ Update ALL_OBJ & VPATH -----------------------
# All source files
ALL_SRC:= \
	$(CORE_SRC) \
	$(DUST_SRC) \
	$(addprefix $(OSLO_PATH)/,$(OSLO_SRC))

# All object files
ALL_OBJ= $(addprefix $(MY_OBJ_DIR)/,$(addsuffix .o,$(basename $(ALL_SRC))))

# All object directories (excluding DUST directory due to compiling rules)
OBJDIRS:= $(MY_OBJ_DIR) $(addprefix $(MY_OBJ_DIR)/, $(ALL_DIRS)) $(MY_MOD_DIR)

# VPATH is used by make (directories to be searched for cmn-files etc)
VPATH:= .$(subst $(space),$(kolon),$(ALL_DIRS))

#-----------------------------------------------------------------
#--- END SOURCE FILES --------------------------------------------
#-----------------------------------------------------------------



#-----------------------------------------------------------------
# COMPILE RULES
#-----------------------------------------------------------------
# For printing title and end note
TITLE=TITLE_COMPILE
DONE=DONE_COMPILE
# Compile .f
$(MY_OBJ_DIR)/%.o: %.f | $(OBJDIRS)
	$(FC) $(FIXEDFLAGS) $(STDOPTS) $(INCLUDES) $(OPTFLAGS) -o $@ -c $<
# Compile .f90
$(MY_OBJ_DIR)/%.o: %.f90 | $(OBJDIRS)
	$(FC) $(FREEFLAGS) $(STDOPTS) $(INCLUDES) $(OPTFLAGS) -o $@ -c $<
# Compile .F90
$(MY_OBJ_DIR)/%.o: %.F90 | $(OBJDIRS)
	$(FC) $(FREEFLAGS) $(STDOPTS) $(INCLUDES) $(OPTFLAGS) $(OSLO_TKN) -D$(HNATIVE) -D$(HWINDOW) -D$(VNATIVE) $(DEGRADE_TKN)  -o $@  -c $<
$(MY_OBJ_DIR)/%.o: %.F | $(OBJDIRS)
	$(FC) $(FIXEDFLAGS) $(STDOPTS) $(INCLUDES) $(OPTFLAGS) $(OSLO_TKN) -D$(HNATIVE) -D$(HWINDOW) -D$(VNATIVE) $(DEGRADE_TKN)  -o $@  -c $<
#--- END of Oslo CTM3 compilation rules --------------------------



#-----------------------------------------------------------------
# COMPILE ALL AND LINK MAIN (do not edit)
#-----------------------------------------------------------------
# Final info to screen
$(DONE): $(MAIN)
	@echo $(PROGNAME) compiled as $(MAIN)
	@echo

# To make sure object dirs are generated (works so far...)
$(MY_OBJ_DIR)/$(OSLO_PATH)/$(OSLODUM_PATH): | $(MY_OBJ_DIR)/$(OSLO_PATH)
$(MY_OBJ_DIR)/$(DUST_PATH): | $(MY_OBJ_DIR)/$(OSLO_PATH)
# Set up  object directories
$(OBJDIRS): | $(TITLE)
	mkdir -p $@

# I need the tokens to work on cmn_size.F90, but have not found out how yet,
# so until that is taken care of, I chose to recompile all if Makefile have changed.
$(ALL_OBJ): $(MAKEFILE)


# Link main
$(MAIN): $(ALL_OBJ)
	@echo
	@echo ALL SHOULD BE OK :: LINKING MAIN
	${FC} -o $@ ${STDOPTS} ${OPTFLAGS} $(INCLUDES) $(ALL_OBJ) ${LDFLAGS}
	@echo

# Print title
$(TITLE):
	@echo
	@echo Compiling Oslo CTM3
	@echo

#-----------------------------------------------------------------
# DEPENDENCIES (edit when adding new files)
#-----------------------------------------------------------------
# Each file may depend on other files. The files you need to list are
# common files and parameter files and also module files.
# Regular subroutine files, i.e. not modules, used by another file
# does not have to be listed.
#
# Files are listed alphabetically for
#   CORE
#   OSLO
#   DUST
#-----------------------------------------------------------------


# === CORE dependencies ==========================================
# cmn_precision.o dependencies
# NONE
#
# cmn_ctm.o dependencies
$(filter %cmn_ctm.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ))
# cmn_met.o dependencies
$(filter %cmn_met.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ))
# cmn_diag.o dependencies
$(filter %cmn_diag.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ))
# cmn_chem.o dependencies
$(filter %cmn_chem.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ))
# cmn_fjx.o dependencies
$(filter %cmn_fjx.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ))
# cmn_sfc.o dependencies
$(filter %cmn_sfc.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ))
# cmn_parameters.o dependencies
$(filter %cmn_parameters.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ))

# averages.o dependencies
$(filter %averages.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# budgets.o dependencies
$(filter %budgets.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))
# cloudjx.o dependencies
$(filter %cloudjx.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ))
# p-cloud2.o dependencies
$(filter %p-cloud2.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cloudjx.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# convection.o dependencies
$(filter %convection.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %cnv_oslo.o, $(ALL_OBJ))
# fastjx.o dependencies
$(filter %fastjx.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# grid.o dependencies
$(filter %grid.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))
# initialize.o dependencies
$(filter %initialize.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cloudjx.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %stt_save_load.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %input_oslo.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ))
# lightning.o dependencies
$(filter %lightning.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# metdata_ecmwf.o dependencies
$(filter %metdata_ecmwf.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cloudjx.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %spectral_routines.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ))
# omp.o dependencies
$(filter %omp.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ))

# p-dyn0.o dependencies
$(filter %p-dyn0.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
# p-dyn2.o dependencies
$(filter %p-dyn2.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %p-vect3.o, $(ALL_OBJ))
# p-linoz.o dependencies
$(filter %p-linoz.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# p-series.o dependencies
$(filter %p-series.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ))
# p-phot_oc.o dependencies
$(filter %p-phot_oc.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %aerosols2fastjx.o, $(ALL_OBJ))
# p-setc_oc.o dependencies
#$(filter %p-setc_oc.o, $(ALL_OBJ)): \
#	$(filter %cmn_precision.o, $(ALL_OBJ)) \
#	$(filter %cmn_size.o, $(ALL_OBJ)) \
#	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
#	$(filter %cmn_chem.o, $(ALL_OBJ)) \
#	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
#	$(filter %cmn_met.o, $(ALL_OBJ)) \
#	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
#	$(filter %utilities.o, $(ALL_OBJ))
# p-vect3.o dependencies
$(filter %p-vect3.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))

# pbl_mixing.o dependencies
$(filter %pbl_mixing.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# regridding.o dependencies
$(filter %regridding.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# scavenging_drydep_uci.o dependencies
$(filter %scavenging_drydep_uci.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %metdata_ecmwf.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# scavenging_largescale_uci.o dependencies
$(filter %scavenging_largescale_uci.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %metdata_ecmwf.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))
# source_uci.o dependencies
$(filter %source_uci.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_aircraft.o, $(ALL_OBJ)) \
	$(filter %emissions_megan.o, $(ALL_OBJ)) \
	$(filter %emissions_ocean.o, $(ALL_OBJ)) \
	$(filter %sulphur_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_volcanoes.o, $(ALL_OBJ))
# spectral_routines.o dependencies
$(filter %spectral_routines.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ))
# steflux.o dependencies
$(filter %steflux.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))
# utilities.o dependencies
$(filter %utilities.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))
# stt_save_load.o dependencies
$(filter %stt_save_load.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \

# pmain.o dependencies
$(filter %pmain.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %averages.o, $(ALL_OBJ)) \
	$(filter %budgets.o, $(ALL_OBJ)) \
	$(filter %cloudjx.o, $(ALL_OBJ)) \
	$(filter %convection.o, $(ALL_OBJ)) \
	$(filter %fastjx.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %initialize.o, $(ALL_OBJ)) \
	$(filter %lightning.o, $(ALL_OBJ)) \
	$(filter %metdata_ecmwf.o, $(ALL_OBJ)) \
	$(filter %omp.o, $(ALL_OBJ)) \
	$(filter %p-linoz.o, $(ALL_OBJ)) \
	$(filter %pbl_mixing.o, $(ALL_OBJ)) \
	$(filter %scavenging_drydep_uci.o, $(ALL_OBJ)) \
	$(filter %scavenging_largescale_uci.o, $(ALL_OBJ)) \
	$(filter %source_uci.o, $(ALL_OBJ)) \
	$(filter %steflux.o, $(ALL_OBJ)) \
	$(filter %stt_save_load.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_aircraft.o, $(ALL_OBJ)) \
	$(filter %emissions_ocean.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %diagnostics_general.o, $(ALL_OBJ)) \
	$(filter %diagnostics_scavenging.o, $(ALL_OBJ)) \
	$(filter %drydeposition_oslo.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_oslo.o, $(ALL_OBJ)) \
	$(filter %fallingaerosols.o, $(ALL_OBJ)) \
	$(filter %gmdump3hrs.o, $(ALL_OBJ)) \
	$(filter %input_oslo.o, $(ALL_OBJ)) \
	$(filter %main_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %ch4routines.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %soa_oslo.o, $(ALL_OBJ)) \
	$(filter %stratchem_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))






# === OSLO dependencies ==========================================
# cmn_oslo.o dependencies
$(filter %cmn_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ))
# eqsam_v03d.o dependencies
$(filter %eqsam_v03d.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ))
# aerosols2fastjx.o dependencies
$(filter %aerosols2fastjx.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))
# atom.o dependencies
$(filter %atom.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %dateconv.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# bcoc_oslo.o dependencies
$(filter %bcoc_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %qssa_integrator.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))
# caribic2.o dependencies
$(filter %caribic2.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %dateconv.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# ch4routines.o dependencies
$(filter %ch4routines.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %emisutils_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))
# chem_oslo_rates.o dependencies
$(filter %chem_oslo_rates.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ)) \
	$(filter %psc_microphysics.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %utilities_troe.o, $(ALL_OBJ))

# cnv_oslo.o dependencies
$(filter %cnv_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %scavenging_largescale_uci.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))

# diagnostics_gemeral.o dependencies
$(filter %diagnostics_general.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %atom.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %caribic2.o, $(ALL_OBJ)) \
	$(filter %ch4routines.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %diagnostics_scavenging.o, $(ALL_OBJ)) \
	$(filter %emissions_megan.o, $(ALL_OBJ)) \
	$(filter %hippo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %satelliteprofiles_mls.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ)) \
	$(filter %troccinox_fal.o, $(ALL_OBJ)) \
	$(filter %troccinox_geo.o, $(ALL_OBJ)) \
	$(filter %troccinox_ban.o, $(ALL_OBJ)) \
	$(filter %verticalprofiles_stations2.o, $(ALL_OBJ))
# diagnostics_scavenging.o dependencies
$(filter %diagnostics_scavenging.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %soa_oslo.o, $(ALL_OBJ))
# drydeposition_oslo.o dependencies
$(filter %drydeposition_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %ch4routines.o, $(ALL_OBJ)) \
	$(filter %soa_oslo.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))
# dust_oslo.o dependencies
$(filter %dust_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dsttvbds.o, $(ALL_OBJ)) \
	$(filter %dstdpsdry.o, $(ALL_OBJ)) \
	$(filter %dstsfc.o, $(ALL_OBJ)) \
	$(filter %dstmbl.o, $(ALL_OBJ)) \
	$(filter %dead_inirun.o, $(ALL_OBJ))
# emisdep4chem_oslo.o dependencies
$(filter %emisdep4chem_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_aircraft.o, $(ALL_OBJ)) \
	$(filter %emissions_megan.o, $(ALL_OBJ)) \
	$(filter %emissions_ocean.o, $(ALL_OBJ)) \
	$(filter %sulphur_oslo.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_volcanoes.o, $(ALL_OBJ))
# emissions_aircraft.o dependencies
$(filter %emissions_aircraft.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))

# emissions_megan.o dependencies
$(filter %emissions_megan.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))

# emissions_ocean.o dependencies
$(filter %emissions_ocean.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %seasaltprod.o, $(ALL_OBJ))
# emissions_oslo.o dependencies
$(filter %emissions_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_aircraft.o, $(ALL_OBJ)) \
	$(filter %emissions_megan.o, $(ALL_OBJ)) \
	$(filter %emissions_ocean.o, $(ALL_OBJ)) \
	$(filter %emisutils_oslo.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %sulphur_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_volcanoes.o, $(ALL_OBJ))
# emissions_volcanoes.o dependencies
$(filter %emissions_volcanoes.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))
# emisutils_oslo.o dependencies
$(filter %emisutils_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_megan.o, $(ALL_OBJ)) \
	$(filter %emissions_ocean.o, $(ALL_OBJ)) \
	$(filter %emissions_volcanoes.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %sulphur_oslo.o, $(ALL_OBJ))
# fallingaerosols.o dependencies
$(filter %fallingaerosols.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ))
# gmdump3hrs.o dependencies
$(filter %gmdump3hrs.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))
# hippo.o dependencies
$(filter %hippo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %dateconv.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# input_oslo.o dependencies
$(filter %input_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %fastjx.o, $(ALL_OBJ)) \
	$(filter %omp.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %scavenging_largescale_uci.o, $(ALL_OBJ)) \
	$(filter %scavenging_drydep_uci.o, $(ALL_OBJ)) \
	$(filter %steflux.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %aerosols2fastjx.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %diagnostics_general.o, $(ALL_OBJ)) \
	$(filter %diagnostics_scavenging.o, $(ALL_OBJ)) \
	$(filter %drydeposition_oslo.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_oslo.o, $(ALL_OBJ)) \
	$(filter %emissions_ocean.o, $(ALL_OBJ)) \
	$(filter %fallingaerosols.o, $(ALL_OBJ)) \
	$(filter %ch4routines.o, $(ALL_OBJ)) \
	$(filter %psc_microphysics.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %nitrate.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %soa_oslo.o, $(ALL_OBJ)) \
	$(filter %sulphur_oslo.o, $(ALL_OBJ))

# main_oslo.o dependencies
$(filter %main_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %averages.o, $(ALL_OBJ)) \
	$(filter %cloudjx.o, $(ALL_OBJ)) \
	$(filter %p-chemflux.o, $(ALL_OBJ)) \
	$(filter %p-linoz.o, $(ALL_OBJ)) \
	$(filter %bcoc_oslo.o, $(ALL_OBJ)) \
	$(filter %diagnostics_general.o, $(ALL_OBJ)) \
	$(filter %drydeposition_oslo.o, $(ALL_OBJ)) \
	$(filter %dust_oslo.o, $(ALL_OBJ)) \
	$(filter %emisdep4chem_oslo.o, $(ALL_OBJ)) \
	$(filter %emisutils_oslo.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %aerosols2fastjx.o, $(ALL_OBJ)) \
	$(filter %psc_microphysics.o, $(ALL_OBJ)) \
	$(filter %nitrate.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %ch4routines.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %strat_aerosols.o, $(ALL_OBJ)) \
	$(filter %stratchem_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ)) \
	$(filter %strat_loss.o, $(ALL_OBJ)) \
	$(filter %strat_o3noy_clim.o, $(ALL_OBJ)) \
	$(filter %tropchem_oslo.o, $(ALL_OBJ))
# psc_microphysics.o dependencies
$(filter %psc_microphysics.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))
# ncutils.o dependencies
$(filter %ncutils.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ))
# nitrate.o dependencies
$(filter %nitrate.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %seasalt.o, $(ALL_OBJ)) \
	$(filter %eqsam_v03d.o, $(ALL_OBJ))
# physics_oslo.o dependencies
$(filter %physics_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %grid.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))
# qssa_integrator.o dependencies
$(filter %qssa_integrator.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ))
# satelliteprofiles_mls.o dependencies
$(filter %satelliteprofiles_mls.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# seasalt.o dependencies
$(filter %seasalt.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %seasaltprod.o, $(ALL_OBJ))

# seasaltprod.o dependencies
$(filter %seasaltprod.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ))
# soa_oslo.o dependencies
$(filter %soa_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ))
# strat_aerosols.o dependencies
$(filter %strat_aerosols.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
# strat_loss.o dependencies
$(filter %strat_loss.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))
# sulphur_oslo.o dependencies
$(filter %sulphur_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))

# stratchem_oslo.o dependencies
$(filter %stratchem_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %diagnostics_general.o, $(ALL_OBJ)) \
	$(filter %emisdep4chem_oslo.o, $(ALL_OBJ)) \
	$(filter %psc_microphysics.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %soa_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ)) \
	$(filter %strat_o3noy_clim.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ)) \
	$(filter %pchemc_str_ij.o, $(ALL_OBJ))

# strat_h2o.o dependencies
$(filter %strat_h2o.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ))

# strat_o3noy_clim.o dependencies
$(filter %strat_o3noy_clim.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %regridding.o, $(ALL_OBJ)) \
	$(filter %utilities.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %ncutils.o, $(ALL_OBJ))

# troccinox_ban.o dependencies
$(filter %troccinox_ban.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %dateconv.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# troccinox_fal.o dependencies
$(filter %troccinox_fal.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %dateconv.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# troccinox_geo.o dependencies
$(filter %troccinox_geo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %dateconv.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# tropchem_oslo.o dependencies
$(filter %tropchem_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_fjx.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_sfc.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %diagnostics_general.o, $(ALL_OBJ)) \
	$(filter %diagnostics_scavenging.o, $(ALL_OBJ)) \
	$(filter %emisdep4chem_oslo.o, $(ALL_OBJ)) \
	$(filter %soa_oslo.o, $(ALL_OBJ)) \
	$(filter %sulphur_oslo.o, $(ALL_OBJ)) \
	$(filter %utilities_oslo.o, $(ALL_OBJ)) \
	$(filter %pchemc_ij.o, $(ALL_OBJ))

# utilities_oslo.o dependencies
$(filter %utilities_oslo.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_chem.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_diag.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))
# verticalprofiles_stations2.o dependencies
$(filter %verticalprofiles_stations2.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ)) \
	$(filter %cmn_met.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ)) \
	$(filter %physics_oslo.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))

# pchemc_ij.o dependencies
$(filter %pchemc_ij.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %qssa_integrator.o, $(ALL_OBJ)) \
	$(filter %strat_h2o.o, $(ALL_OBJ))
# pchemc_str_ij.o dependencies
$(filter %pchemc_str_ij.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ)) \
	$(filter %chem_oslo_rates.o, $(ALL_OBJ)) \
	$(filter %qssa_integrator.o, $(ALL_OBJ))




# === DUST dependencies ==========================================

# blmutl.o dependencies
$(filter %blmutl.o, $(ALL_OBJ)): \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstlsm.o, $(ALL_OBJ)) \
	$(filter %dstblm.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dbg_mdl.o dependencies
#	None
# dead_inirun.o dependencies
$(filter %dead_inirun.o, $(ALL_OBJ)): \
	$(filter %dstcmnini.o, $(ALL_OBJ)) \
	$(filter %dstbdg.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstctl.o, $(ALL_OBJ)) \
	$(filter %dstdbg.o, $(ALL_OBJ)) \
	$(filter %dstnm.o, $(ALL_OBJ)) \
	$(filter %dstpsd.o, $(ALL_OBJ)) \
	$(filter %dsttibds.o, $(ALL_OBJ)) \
	$(filter %dsttvbds.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstaer.o dependencies
$(filter %dstaer.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ))
# dstbdg.o dependencies
$(filter %dstbdg.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dbg_mdl.o, $(ALL_OBJ)) \
	$(filter %nf90_utl.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ)) \
	$(filter %vec_mdl.o, $(ALL_OBJ)) \
	$(filter %dbg_mdl.o, $(ALL_OBJ))
# dstblm.o dependencies
$(filter %dstblm.o, $(ALL_OBJ)): \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstchm.o dependencies
$(filter %dstchm.o, $(ALL_OBJ)): \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstcmnini.o dependencies
$(filter %dstcmnini.o, $(ALL_OBJ)): \
	$(filter %dstblm.o, $(ALL_OBJ)) \
	$(filter %dstchm.o, $(ALL_OBJ)) \
	$(filter %dstdbg.o, $(ALL_OBJ)) \
	$(filter %dstodx.o, $(ALL_OBJ)) \
	$(filter %dstscv.o, $(ALL_OBJ)) \
	$(filter %dstsltsbl.o, $(ALL_OBJ))
# dstcst.o dependencies
$(filter %dstcst.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %cmn_parameters.o, $(ALL_OBJ))
# dstctl.o dependencies
$(filter %dstctl.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ))
# dstdbg.o dependencies
$(filter %dstdbg.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstdpsdry.o dependencies
$(filter %dstdpsdry.o, $(ALL_OBJ)): \
	$(filter %blmutl.o, $(ALL_OBJ)) \
	$(filter %dpsdryutl.o, $(ALL_OBJ)) \
	$(filter %dstaer.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstdbg.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dstmssutl.o, $(ALL_OBJ)) \
	$(filter %dstsfc.o, $(ALL_OBJ)) \
	$(filter %dead_history.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dpsdryutl.o dependencies
$(filter %dpsdryutl.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ))
# dstgrd.o dependencies
$(filter %dstgrd.o, $(ALL_OBJ)): \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstlsm.o dependencies
$(filter %dstlsm.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstmbl.o dependencies
$(filter %dstmbl.o, $(ALL_OBJ)): \
	$(filter %blmutl.o, $(ALL_OBJ)) \
	$(filter %dstaer.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstdbg.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dstmblutl.o, $(ALL_OBJ)) \
	$(filter %dstmssutl.o, $(ALL_OBJ)) \
	$(filter %dstnm.o, $(ALL_OBJ)) \
	$(filter %dstsfc.o, $(ALL_OBJ)) \
	$(filter %dstsltsbl.o, $(ALL_OBJ)) \
	$(filter %dsttvbds.o, $(ALL_OBJ)) \
	$(filter %dead_history.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ)) \
	$(filter %wbl_mdl.o, $(ALL_OBJ))
# dstmblutl.o dependencies
$(filter %dstmblutl.o, $(ALL_OBJ)): \
	$(filter %blmutl.o, $(ALL_OBJ)) \
	$(filter %dbg_mdl.o, $(ALL_OBJ)) \
	$(filter %dstblm.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dstlsm.o, $(ALL_OBJ)) \
	$(filter %dstpsd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstmssutl.o dependencies
$(filter %dstmssutl.o, $(ALL_OBJ)): \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %vec_mdl.o, $(ALL_OBJ))
# dstnm.o dependencies
$(filter %dstnm.o, $(ALL_OBJ)): \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ))
# dstodx.o dependencies
$(filter %dstodx.o, $(ALL_OBJ)): \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dstpsd.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %vec_mdl.o, $(ALL_OBJ)) \
	$(filter %xtr_mdl.o, $(ALL_OBJ))
# dstpsd.o dependencies
$(filter %dstpsd.o, $(ALL_OBJ)): \
	$(filter %dbg_mdl.o, $(ALL_OBJ)) \
	$(filter %dstaer.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstctl.o, $(ALL_OBJ)) \
	$(filter %dstdryutl.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %nf90_utl.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %psdlgn.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ)) \
	$(filter %vec_mdl.o, $(ALL_OBJ)) \
	$(filter %cmn_oslo.o, $(ALL_OBJ))
# dstsfc.o dependencies
$(filter %dstsfc.o, $(ALL_OBJ)): \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dstscv.o dependencies
$(filter %dstscv.o, $(ALL_OBJ)): \
	$(filter %dstaer.o, $(ALL_OBJ)) \
	$(filter %dstcst.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %vec_mdl.o, $(ALL_OBJ)) \
	$(filter %xtr_mdl.o, $(ALL_OBJ))
# dstsltsbl.o dependencies
$(filter %dstsltsbl.o, $(ALL_OBJ)): \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dsttibds.o dependencies
$(filter %dsttibds.o, $(ALL_OBJ)): \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dstgrd.o, $(ALL_OBJ)) \
	$(filter %dstsfc.o, $(ALL_OBJ)) \
	$(filter %nf90_utl.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dsttvbds.o dependencies
$(filter %dsttvbds.o, $(ALL_OBJ)): \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dstdbg.o, $(ALL_OBJ)) \
	$(filter %nf90_utl.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %utl_mdl.o, $(ALL_OBJ))
# erf_mdl.o dependencies
$(filter %erf_mdl.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ))
# gmm_mdl.o dependencies
$(filter %gmm_mdl.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ))
# dead_history.o dependencies
$(filter %dead_history.o, $(ALL_OBJ)): \
	$(filter %pmgrid.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# nf90_utl.o dependencies
$(filter %nf90_utl.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %dbg_mdl.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ))
# pmgrid.o dependencies
$(filter %pmgrid.o, $(ALL_OBJ)): \
	$(filter %cmn_size.o, $(ALL_OBJ)) \
	$(filter %cmn_ctm.o, $(ALL_OBJ))
# dead_precision.o dependencies
$(filter %dead_precision.o, $(ALL_OBJ)): \
	$(filter %cmn_precision.o, $(ALL_OBJ))

# psdlgn.o dependencies
$(filter %psdlgn.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %erf_mdl.o, $(ALL_OBJ))
# sng_mdl.o dependencies
$(filter %sng_mdl.o, $(ALL_OBJ)): \
	$(filter %dbg_mdl.o, $(ALL_OBJ))
# utl_mdl.o dependencies
$(filter %utl_mdl.o, $(ALL_OBJ)): \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %dbg_mdl.o, $(ALL_OBJ))
# vec_mdl.o dependencies
$(filter %vec_mdl.o, $(ALL_OBJ)): \
	$(filter %dbg_mdl.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ)) \
	$(filter %utl_mdl.o, $(ALL_OBJ)) \
	$(filter %xtr_mdl.o, $(ALL_OBJ))
# wbl_mdl.o dependencies
$(filter %wbl_mdl.o, $(ALL_OBJ)): \
	$(filter %gmm_mdl.o, $(ALL_OBJ)) \
	$(filter %dead_precision.o, $(ALL_OBJ))
# xtr_mdl.o dependencies
$(filter %xtr_mdl.o, $(ALL_OBJ)): \
	$(filter %dbg_mdl.o, $(ALL_OBJ)) \
	$(filter %sng_mdl.o, $(ALL_OBJ))




#-----------------------------------------------------------------
# CLEAN
#-----------------------------------------------------------------
.PHONY: clean
clean :
	@echo Cleaning make of Oslo CTM3
	rm -f fort.* core *.mod  $(MAIN) $(TITLE) $(DONE)
	rm -f cmn_size_mod.f90
	rm -rf $(MY_OBJ_DIR) $(MY_MOD_DIR)


#-----------------------------------------------------------------
# CHECK
#-----------------------------------------------------------------
.PHONY: check
check:
	@echo "MACHINE:    " $(UNAMEM)
	@echo "NODE NAME:  " $(UNAMEN)
	@echo "KERNEL NAME:" $(UNAMES)
	@echo "COMPILED BY:" $(MYNAME)
	@echo "MAKEFILE:   " $(MAKEFILE)
	@echo
	@echo "INCLUDES: " $(INCLUDES)
	@echo "NETCDF_INC: " $(NETCDF_INC)
	@echo "NETCDF_LIB: " $(NETCDF_LIB)
	@echo "NETCDF_LLIBS: " $(NETCDF_LLIBS)
	@echo "VPATH: " $(VPATH)
	@echo
	@echo "ALL_DIRS: " $(ALL_DIRS)
	@echo
	@echo "ALL_SRC: " $(ALL_SRC)
	@echo
	@echo "ALL_OBJ: " $(ALL_OBJ)
	@echo
	@echo "MY_OBJ_DIR: " $(MY_OBJ_DIR)
	@echo "OBJDIRS: " $(OBJDIRS)
	@echo
	@echo "DUST_SRC: " $(DUST_SRC)
	@echo "DUST_OBJ: " $(DUST_OBJ)
	@echo "DUSTOBJDIRS: " $(DUSTOBJDIRS)
	@echo "HNATIVE: " $(HNATIVE)
	@echo "HWINDOW: " $(HWINDOW)
	@echo "VNATIVE: " $(VNATIVE)

##################################################################
# END Oslo CTM3 Makefile
##################################################################
