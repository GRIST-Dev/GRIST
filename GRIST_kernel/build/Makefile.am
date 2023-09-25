
##################################################################
## PLEASE DO NOT MODIFY BELOW UNLESS YOU KNOW WHAT YOU ARE DOING##
##################################################################

UTOMAKE_OPTIONS = foreign

NETCDF   =${NETCDF_PATH}
PNETCDF  =${PNETCDF_PATH}
LAPACK   =${LAPACK_PATH}
METIS    =${METIS_LIB_PATH}
NCVERSION=${netcdf_version}
LIBTORCH =${LIBTORCH_PATH}

if IS_FTORCH

AM_FFLAGS =  -I$(NETCDF)/include/ -I$(PNETCDF)/include/ -I$(LAPACK) -DUSE_HALO2 -I${LIBTORCH}/include/ftorch/
if NETCDF_V3
AM_LDFLAGS= -L$(PNETCDF)/lib/ -L$(NETCDF)/lib/ -L${LAPACK} -llapack -lrefblas -L${METIS} -lnetcdf -L${LIBTORCH}/lib64/ -lftorch
endif
if NETCDF_V4
AM_LDFLAGS= -L$(PNETCDF)/lib/ -L$(NETCDF)/lib/ -L${LAPACK} -llapack -lrefblas -L${METIS} -lnetcdf -lnetcdff -L${LIBTORCH}/lib64/ -lftorch
endif

else

AM_FFLAGS =  -I$(NETCDF)/include/ -I$(PNETCDF)/include/ -I$(LAPACK) -DUSE_HALO2
if NETCDF_V3
AM_LDFLAGS= -L$(PNETCDF)/lib/ -L$(NETCDF)/lib/ -L${LAPACK} -llapack -lrefblas -L${METIS} -lnetcdf
endif
if NETCDF_V4
AM_LDFLAGS= -L$(PNETCDF)/lib/ -L$(NETCDF)/lib/ -L${LAPACK} -llapack -lrefblas -L${METIS} -lnetcdf -lnetcdff
endif

endif

AM_CPPFLAGS= -O3 -DBOOST_LOG_DYN_LINK -fopenmp -fPIC --std=c++0x -Werror=return-type -fno-trapping-math

FC = ${Fortran_Compiler}
CC = ${FC}
CXX = ${CXX_Compiler}

VPATH=$(shell cat Filepath)
BUILD_DIR = ./
SOURCES = $(shell cat Srcfiles)
ALLOBJS = $(addprefix $(BUILD_DIR),$(addsuffix .o,$(basename $(SOURCES))))
bin_PROGRAMS_DESTDIR  = ${prefix}
oplevel3= -O3 -DBOOST_LOG_DYN_LINK -fopenmp -fno-trapping-math

##################################################################
## If new Entry is needed, please add following the template    ##
##################################################################

if BUILD_GRIST_SWM
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_swe.o,$(ALLOBJS))

bin_PROGRAMS=ParGRIST_swm.exe partition.exe
ParGRIST_swm_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
ParGRIST_swm_exe_LDADD = $(OBJS)  $(AM_LDFLAGS)  -lmetis -lmpi -lstdc++ -lpnetcdf
partition_exe_LDADD= $(PART_OBJS) $(AM_LDFLAGS)  -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel = -O1 
endif

if BUILD_GRIST_GCM
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_atmos.o,$(ALLOBJS))

bin_PROGRAMS=ParGRIST_gcm.exe partition.exe
ParGRIST_gcm_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
ParGRIST_gcm_exe_LDADD = $(OBJS)  $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf				 
partition_exe_LDADD= $(PART_OBJS) $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel = -O1 -DCDATE
endif

if BUILD_GRIST_AMIPC
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_atmos.o,$(ALLOBJS))

bin_PROGRAMS=ParGRIST_amipc.exe partition.exe
ParGRIST_amipc_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
ParGRIST_amipc_exe_LDADD = $(OBJS)  $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf				 
partition_exe_LDADD= $(PART_OBJS) $(AM_LDFLAGS)   -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel = -O1 -DSPIO -DCMAPI -DCDATE -DAMIPC_PHYSICS -DUSE_NOAHMP -DOCNABD -DCAM3OCNABD
endif

if BUILD_GRIST_AMIPW
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_atmos.o,$(ALLOBJS))

bin_PROGRAMS=ParGRIST_amipw.exe partition.exe
ParGRIST_amipw_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
ParGRIST_amipw_exe_LDADD = $(OBJS)  $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf				 
partition_exe_LDADD= $(PART_OBJS) $(AM_LDFLAGS)   -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel = -O1 -DRRTMG_V381 -DSCALE_TBCU -DSPIO -DUSE_LEAP_YEAR \
          -convert big_endian -r8 -DAMIPW_PHYSICS -DAMIPW_CLIMATE -DUSE_NOAHMP -DCDATE
endif

if BUILD_GRIST_LAM_AMIPW
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_atmos.o,$(ALLOBJS))

bin_PROGRAMS=ParGRIST_lam_amipw.exe partition.exe
ParGRIST_lam_amipw_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
ParGRIST_lam_amipw_exe_LDADD = $(OBJS)  $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf				 
partition_exe_LDADD= $(PART_OBJS) $(AM_LDFLAGS)   -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel =-DLAM_DOMAIN -O1 -DRRTMG_V381 -DSCALE_TBCU -DSPIO -DUSE_LEAP_YEAR \
          -convert big_endian -r8 -DAMIPW_PHYSICS -DAMIPW_CLIMATE -DUSE_NOAHMP -DCDATE
endif

if BUILD_GRIST_SCM_PHYSC
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_scm.o,$(ALLOBJS))

bin_PROGRAMS=GRIST_scm_physc.exe partition.exe
GRIST_scm_physc_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
GRIST_scm_physc_exe_LDADD = $(OBJS)  $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf				 
partition_exe_LDADD=    $(PART_OBJS) $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel =-O1 -DSPIO -DCDATE -DAMIPC_PHYSICS -DSCM_PHYSICS #-DCMAPI -DUSE_NOAHMP -DOCNABD -DCAM3OCNABD
endif

if BUILD_GRIST_SCM_PHYSW
OBJS = $(filter-out ./grist_partition.o,$(ALLOBJS))
PART_OBJS = $(filter-out ./grist_scm.o,$(ALLOBJS))

bin_PROGRAMS=GRIST_scm_physw.exe partition.exe
GRIST_scm_physw_exe_SOURCES = $(OBJS)
partition_exe_SOURCES =  $(PART_OBJS)
GRIST_scm_physw_exe_LDADD = $(OBJS)  $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf				 
partition_exe_LDADD=    $(PART_OBJS) $(AM_LDFLAGS) -lmetis -lmpi -lstdc++ -lpnetcdf
oplevel =-O1 -DRRTMG_V381 -DSCALE_TBCU -DSPIO -DUSE_LEAP_YEAR \
          -convert big_endian -r8 -DAMIPW_PHYSICS -DAMIPW_CLIMATE -DSCM_PHYSICS -DCDATE
endif

##################################################################
## PLEASE DO NOT MODIFY BELOW UNLESS YOU KNOW WHAT YOU ARE DOING##
##################################################################

.cpp.o:
	$(CXX) $(AM_CPPFLAGS) -c $< -o $@

grist_mpi.o: grist_mpi.F90
	$(FC) ${oplevel3} -c $< -o $@
grist_map.o: grist_map.F90
	$(FC) ${oplevel3} -c $< -o $@
grist_metis.o: grist_metis.F90
	$(FC) ${oplevel3} -c $< -o $@
grist_set.o: grist_set.F90
	$(FC) ${oplevel3} -c $< -o $@
grist_lib.o: grist_lib.F90
	$(FC) ${oplevel3} -c $< -o $@
grist_utils.o: grist_utils.F90
	$(FC) ${oplevel3} -c $< -o $@

.F.o:
	$(FC) ${oplevel} $(AM_FFLAGS) -c $< -o $@
.f.o:
	$(FC) ${oplevel} $(AM_FFLAGS) -c $< -o $@
.F90.o:
	$(FC) ${oplevel} $(AM_FFLAGS) -c $< -o $@
.f90.o:
	$(FC) ${oplevel} $(AM_FFLAGS) -c $< -o $@

install-exec-hook:
	@$(NORMAL_INSTALL)
	test -z "$(DESTDIR)$(bin_PROGRAMS_DESTDIR)" || $(MKDIR_P) "$(DESTDIR)$(bin_PROGRAMS_DESTDIR)"
	$(INSTALL_PROGRAM) $(bin_PROGRAMS) "$(DESTDIR)$(bin_PROGRAMS_DESTDIR)"
	@echo "========================================================================"        
	@echo "--->    Executables successfully built $(prefix)/$(bin_PROGRAMS)   <---"

include Depends

install-exec-hook:
	@$(NORMAL_INSTALL)
	test -z "$(DESTDIR)$(bin_PROGRAMS_DESTDIR)" || $(MKDIR_P) "$(DESTDIR)$(bin_PROGRAMS_DESTDIR)"
	$(INSTALL_PROGRAM) $(bin_PROGRAMS) "$(DESTDIR)$(bin_PROGRAMS_DESTDIR)"
	@echo "========================================================================"        
	@echo "--->    Executables successfully built $(prefix)/$(bin_PROGRAMS)   <---"

include Depends

RM  = rm
.PHONY : clean
clean:
	-$(RM) -f $(ALLOBJS)  *.mod *exe $(prefix)/*exe  #$(EXEDIR)/$(EXENAME) #$(EXEDIR)/$(PARTEXE)
