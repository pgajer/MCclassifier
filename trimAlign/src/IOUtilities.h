
/* All but, _fOpen() routines are modifications of
 *
 * Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 *
 */

#if !defined(IOUTILITIES_H)
#define IOUTILITIES_H

#include "msa.h"
#include "CUtilities.h"

#define fOpen(x,y)   _fOpen((x), (y), __FILE__, __LINE__)
#define isgap(c) ((c) == '.' || (c) == '-')

char * sre_fgets( MSArec_t *rp );
MSArec_t * openFasta(char *seqfile);
void closeFasta(MSArec_t *rp);
int readMSArecord( MSArec_t *rp );
FILE *_fOpen ( const char *filetoopen, const char *format, const char * srcfile, int line );

#endif
