#/bin/bash -f
model=$1
PREFIX=/fs2/home/zhangyi/zhangyi/ai-phys/GRIST_kernel-main/bin/
export NETCDF_PATH=/fs2/software/netcdf/4.8.0-icc19.0-openmpi/
export PNETCDF_PATH=/fs2/software/pnetcdf/1.12.2-icc19.0-openmpi/
export LAPACK_PATH=/fs2/home/zhangyi/softwares/lapack-3.8.0/
export METIS_LIB_PATH=/fs2/home/zhangyi/softwares/metis-5.1.0/lib/
export LIBTORCH_PATH=/fs2/home/zhangyi/wangym/model/test/libtorch/FTorch-intel

export Fortran_Compiler="mpifort"
export CXX_Compiler=mpicxx
export netcdf_version="4"  #this is to distinguish netcdf code version(before 4.1 set 3; otherwise set 4)
export use_ftorch="no"    #

##################################################################
## PLEASE DO NOT MODIFY BELOW UNLESS YOU KNOW WHAT YOU ARE DOING##
##################################################################

if [[ ${model} != "amipw" && \
      ${model} != "amipc" && \
      ${model} != "gcm"   && \
      ${model} != "scm_physc" && \
      ${model} != "scm_physw" ]];then
echo "A proper model prompt for GRIST is: amipw, amipc, gcm, scm_physc, scm_physw"
echo "Your input : "${model}" is not one of them"
exit
fi

if [[ ${model} = "amipc" ]] ;then
cp ./path/Filepath.GCM_AMIPC Filepath
fi
if [[ ${model} = "amipw" || ${model} = "lam_amipw" ]] ;then
cp ./path/Filepath.GCM_AMIPW Filepath
fi
if [[ ${model} = "gcm" ]] ;then
cp ./path/Filepath.GCM   Filepath
fi
if [[ ${model} = "swm" ]] ;then
cp ./path/Filepath.SWM   Filepath
fi
if [[ ${model} = "scm_physc" ]] ;then
cp ./path/Filepath.SCM_PhysC   Filepath
fi
if [[ ${model} = "scm_physw" ]] ;then
cp ./path/Filepath.SCM_PhysW   Filepath
fi

make clean
make distclean
./mkSrcfiles
./mkDepends Filepath Srcfiles >Depends

aclocal
autoconf
autoheader
automake --add-missing

if [[ ${model} = "swm" ]] ;then
./configure --with-grist_swm --prefix=${PREFIX}
fi

if [[ ${model} = "gcm" ]] ;then
./configure --with-grist_gcm --prefix=${PREFIX}
fi

if [[ ${model} = "amipw" ]] ;then
./configure --with-grist_amipw --prefix=${PREFIX}
fi

if [[ ${model} = "lam_amipw" ]] ;then
./configure --with-grist_lam_amipw --prefix=${PREFIX}
fi

if [[ ${model} = "amipc" ]] ;then
./configure --with-grist_amipc --prefix=${PREFIX}
fi

if [[ ${model} = "scm_physc" ]] ;then
./configure --with-grist_scm_physc --prefix=${PREFIX}
fi

if [[ ${model} = "scm_physw" ]] ;then
./configure --with-grist_scm_physw --prefix=${PREFIX}
fi

make -j8

make install-exec-hook
