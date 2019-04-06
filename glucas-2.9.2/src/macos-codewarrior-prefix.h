/*
  (c) 2001,2002 Klaus Kastens
 
	"Prefix File" for Metrowerks CodeWarrior.
	This file is included in every source file, and contains macro settings
	similar to 'gcc -Dmacro'.
*/
#define Y_MEM_THRESHOLD 8192
#define Y_BLOCKSIZE 4096
#define Y_SHIFT 6
#define Y_LONG_MACROS
#define Y_VECTORIZE
#define Y_VECTORIZE2
#define Y_KILL_BRANCHES
/* use '#define Y_TARGET 0' to disable prefetch hinting */
/* use '#define Y_TARGET 21' to enable prefetch hinting for PowerPC 601 */
/* use '#define Y_TARGET 23' to enable prefetch hinting for PowerPC 604e, 7xx, 74xx */
#define Y_TARGET 23
