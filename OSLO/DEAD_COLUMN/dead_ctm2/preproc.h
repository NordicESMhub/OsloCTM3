/* $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/preproc.h,v 1.1 2003/04/15 14:41:36 alfgr Exp $ */

#ifndef PREPROC_SET
#define PREPROC_SET
#define COUP_CCM

#ifdef T5
#define LSMLON 16
#define LSMLAT 8
#endif /* not T5 */

#ifdef T21
#define LSMLON 64
#define LSMLAT 32
#endif /* not T21 */

#ifdef T31
#define LSMLON 96
#define LSMLAT 48
#endif /* not T31 */

#ifdef T42
#define LSMLON 128
#define LSMLAT 64
#endif /* not T42 */

#ifdef T62
#define LSMLON 192
#define LSMLAT 94
#endif /* not T62 */
#endif /* PREPROC_SET */


