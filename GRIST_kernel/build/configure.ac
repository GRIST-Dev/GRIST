#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([GRIST_MODEL], [A23.10.20], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([../src/infrastructure/utils/grist_list.F90])
AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT([/usr/local])

# Checks for programs.
AC_PROG_FC
#AC_ARG_ENABLE

# Checks for libraries.
AC_PROG_RANLIB
# Check for ftorch
AM_CONDITIONAL([IS_FTORCH],[test "x$use_ftorch" = "xyes"])
# Checks for netcdf version
AM_CONDITIONAL([NETCDF_V3],[test "x$netcdf_version" = "x3"])
AM_CONDITIONAL([NETCDF_V4],[test "x$netcdf_version" = "x4"])

# config options---control compile and dep 定义WITH_变量
AC_ARG_WITH([grist_swm], 
	[AS_HELP_STRING([--with-grist_swm], [Build with grist_swm options])])
AM_CONDITIONAL([BUILD_GRIST_SWM],[test "x$with_grist_swm" = "xyes"])

AC_ARG_WITH([grist_gcm], 
	[AS_HELP_STRING([--with-grist_gcm], [Build with grist_gcm options])])
AM_CONDITIONAL([BUILD_GRIST_GCM],[test "x$with_grist_gcm" = "xyes"])

AC_ARG_WITH([grist_amipc], 
	[AS_HELP_STRING([--with-grist_amipc], [Build with grist_amipc options])])
AM_CONDITIONAL([BUILD_GRIST_AMIPC],[test "x$with_grist_amipc" = "xyes"])

AC_ARG_WITH([grist_amipw], 
	[AS_HELP_STRING([--with-grist_amipw], [Build with grist_amipw options])])
AM_CONDITIONAL([BUILD_GRIST_AMIPW],[test "x$with_grist_amipw" = "xyes"])

AC_ARG_WITH([grist_lam_amipw], 
	[AS_HELP_STRING([--with-grist_lam_amipw], [Build with grist_lam_amipw options])])
AM_CONDITIONAL([BUILD_GRIST_LAM_AMIPW],[test "x$with_grist_lam_amipw" = "xyes"])

AC_ARG_WITH([grist_scm_physc], 
	[AS_HELP_STRING([--with-grist_scm_physc], [Build with grist_scm_physc options])])
AM_CONDITIONAL([BUILD_GRIST_SCM_PHYSC],[test "x$with_grist_scm_physc" = "xyes"])

AC_ARG_WITH([grist_scm_physw], 
	[AS_HELP_STRING([--with-grist_scm_physw], [Build with grist_scm_physw options])])
AM_CONDITIONAL([BUILD_GRIST_SCM_PHYSW],[test "x$with_grist_scm_physw" = "xyes"])

# 根据配置参数进行相应的设置
if test "x$with-grist_swm" = "xyes";then
	AC_DEFINE([BUILD_GRIST_SWM],[],[Build grist_swm])
	GRIST_SWM_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_SWM_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_SWM_OUTPUT_DIR/configure])
fi

if test "x$with-grist_gcm" = "xyes";then
	AC_DEFINE([BUILD_GRIST_GCM],[],[Build grist_gcm])
	GRIST_GCM_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_GCM_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_GCM_OUTPUT_DIR/configure])
fi

if test "x$with-grist_amipc" = "xyes";then
	AC_DEFINE([BUILD_GRIST_AMIPC],[],[Build grist_amipc])
	GRIST_AMIPC_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_AMIPC_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_AMIPC_OUTPUT_DIR/configure])
fi

if test "x$with-grist_amipw" = "xyes";then
	AC_DEFINE([BUILD_GRIST_AMIPW],[],[Build grist_amipw])
	GRIST_AMIPW_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_AMIPW_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_AMIPW_OUTPUT_DIR/configure])
fi

if test "x$with-grist_lam_amipw" = "xyes";then
	AC_DEFINE([BUILD_GRIST_LAM_AMIPW],[],[Build grist_lam_amipw])
	GRIST_LAM_AMIPW_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_LAM_AMIPW_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_LAM_AMIPW_OUTPUT_DIR/configure])
fi

if test "x$with-grist_scm_physc" = "xyes";then
	AC_DEFINE([BUILD_GRIST_SCM_PHYSC],[],[Build grist_scm_physc])
	GRIST_SCM_PHYSC_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_SCM_PHYSC_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_SCM_PHYSC_OUTPUT_DIR/configure])
fi

if test "x$with-grist_scm_physw" = "xyes";then
	AC_DEFINE([BUILD_GRIST_SCM_PHYSW],[],[Build grist_scm_physw])
	GRIST_SCM_PHYSW_OUTPUT_DIR="../bin/"
	AC_CONFIG_FILES([$GRIST_SCM_PHYSW_OUTPUT_DIR/Makefile])
	AC_OUTPUT([$GRIST_SCM_PHYSW_OUTPUT_DIR/configure])
fi

#AC_SUBST([FFLAGS_GRIST_GCM])
# Checks for header files.
# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT()
