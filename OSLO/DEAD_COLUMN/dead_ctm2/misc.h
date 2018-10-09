/* $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/misc.h,v 1.1 2003/04/15 14:41:36 alfgr Exp $ */

#ifndef MISC_SET
#define MISC_SET
#define NCPREC NF_FLOAT
#define PVP 
#ifdef CRAY
#define REALTYPE MPI_REAL
#undef  SHELL_MSS
#undef  FORTFFT
#else /* not CRAY */
#define REALTYPE MPI_DOUBLE_PRECISION
#define SHELL_MSS
#define FORTFFT
#endif /* not CRAY */
#undef  COUP_SOM
#undef  COUP_CSM
#undef  SPMD
#if ( ! defined AIX ) && ( ! defined CRAY ) && ( ! defined DEC ) && ( ! defined LINUX ) && ( ! defined OSF1 ) && ( ! defined SGI ) && ( ! defined SUN ) && ( ! defined T3D )
You must define one of AIX, CRAY, DEC, LINUX, OSF1, SGI, SUN, or T3D
#endif /* not CCM_ARCH */ 
#endif /* not MISC_SET */


