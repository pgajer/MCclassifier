
/* All but, _fOpen() routines are modifications of
 *
 * Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 *
 */

#if !defined(IOUTILITIES_H)
#define IOUTILITIES_H

#include "fasta.h"
#include "CUtilities.h"

#define fOpen(x,y)   _fOpen((x), (y), __FILE__, __LINE__)
#define isgap(c) ((c) == '.' || (c) == '-')

char * sre_fgets( FASTArec_t *rp );
FASTArec_t * openFasta(char *seqfile);
void closeFasta(FASTArec_t *rp);
int readFASTArecord( FASTArec_t *rp );
FILE *_fOpen ( const char *filetoopen, const char *format, const char * srcfile, int line );

#endif
