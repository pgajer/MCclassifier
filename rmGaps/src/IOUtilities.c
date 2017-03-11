// IOUtilities.c -
//
// Pawel Gajer
// Wednesday, September 26 2012
//

#include "IOUtilities.h"

//----------------------------------------------------------------- openFasta ----
// initializes a pointer to FASTArec_t and reads the first MSA record
FASTArec_t * openFasta(char *seqfile)
{
  FASTArec_t *rp = MallocOrDie(sizeof(FASTArec_t));

  rp->bsize    = 1024;
  rp->buflen   = rp->bsize; // initialize buffer
  rp->buf      = MallocOrDie(sizeof(char) * rp->buflen);
  rp->seqlen   = rp->bsize; // initialize sequence
  rp->seq      = MallocOrDie(sizeof(char) * rp->seqlen);

  rp->fp = fOpen(seqfile, "r");

  if (rp->fp == NULL) { free(rp); return NULL; }

  if ( sre_fgets( rp ) == NULL )
  {
    free(rp);
    return NULL;
  }

  //printf("line %d, buflen: %d\n",0 , rp->buflen);

  return rp;
}

//----------------------------------------------------------------- sre_fgets ----
char * sre_fgets( FASTArec_t *rp )
{
  #define DEBUGS 0

  char **buf = &(rp->buf);
  int *n     = &(rp->buflen);
  FILE *fp   = rp->fp;

  char *s;
  int   len;
  int   pos;

  /* Simple case 1. We're sitting at EOF, or there's an error.
   *                fgets() returns NULL, so we return NULL.
   */
  if (fgets(*buf, *n, fp) == NULL)
  {
    #if DEBUGS
    printf("case 1\n");
    #endif
    return NULL;
  }


  /* Simple case 2. fgets() got a string, and it reached EOF.
   *                return success status, so caller can use
   *                the last line; on the next call we'll
   *                return the 0 for the EOF.
   */
  if (feof(fp))
  {
    #if DEBUGS
    printf("case 2\n");
    #endif
    return *buf;
  }

  /* Simple case 3. We got a complete string, with \n,
   *                and don't need to extend the buffer.
   */
  len = strlen(*buf);
  if ((*buf)[len-1] == '\n')
  {
    (*buf)[len-1] = '\0'; // I don't want \n to be part of the buffer string

    #if DEBUGS
    printf("case 3\n");
    #endif

    return *buf;
  }


  /* The case we're waiting for. We have an incomplete string, and we have to
   * extend the buffer one or more times. Make sure we overwrite the previous
   * fgets's \0 (hence +(n-1) in first step, rather than bsize, and reads of
   * (bsize+1), not bsize).
   */
  pos = (*n)-1;
  #if DEBUGS
  printf("case 4, pos: %d\n", pos);
  #endif
  while (1)
  {
    *n  += rp->bsize;
    *buf = ReallocOrDie(*buf, sizeof(char) * (*n));

    s = *buf + pos; // s starts at the '\0' of the previous buffer

    if (fgets(s, (rp->bsize + 1), fp) == NULL) return *buf;

    len = strlen(s);
    if (s[len-1] == '\n')
    {
      s[len-1] = '\0'; // I don't want \n to be part of the buffer string
      return *buf;
    }

    pos += rp->bsize;
  }
  /*NOTREACHED*/
}

//----------------------------------------------------------------- closeFasta ----
void closeFasta(FASTArec_t *rp)
{
  fclose(rp->fp);

  free(rp->buf);
  free(rp->seq);
  free(rp->name);
  free(rp);
}


//--------------------------------------------- readFASTArecord ----------------------
int readFASTArecord( FASTArec_t *rp )
{
  // rp->seq holds only non-gap characters

  // Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
  if (rp->buf[0] != '>') return 0;

  // Parse out the name: the first non-whitespace token after the >
  char * s  = strtok(rp->buf+1, " \t\n");
  rp->name = MallocOrDie(sizeof(char) * (strlen(s)+1));
  strcpy(rp->name, s);

  // Parse out the sequence and the start and end of non-gap characters
  int n = 0;

  while ( sre_fgets( rp ) != NULL )
  {
    if (rp->buf[0] == '>') break;	// reached the next descline

    for (s = rp->buf; *s != '\0'; s++)
    {
      // here one could check if the buffer consists of characters of nucleotide
      // names + gap symbols but I don't do it b/c I am going to use this on the
      // output of some reliable software that generates correct output

      //printf("*s=%c, gap=%d, isgap(*s)=%d, seenNuc=%d, seqStart=%d\n", *s, gap, isgap(*s), seenNuc, rp->seqStart);

      if ( isgap(*s) ) continue;

      rp->seq[n] = *s;          // store the character, bump length n
      n++;
      if (rp->seqlen == n)	// are we out of room in seq? if so, expand
      {			        // (remember, need space for the final '\0')
	rp->seqlen += rp->bsize;
	rp->seq = ReallocOrDie(rp->seq, sizeof(char) * rp->seqlen);
      }
    }
  }
  rp->seq[n] = '\0';
  rp->seqlen = n;

  return 1;
}


//----------------------------------------------------------------- fOpen ----
//! open a file handle and throw an error message if it cannot be opened
FILE *_fOpen ( const char *file, const char *format, const char * cppfile, int line )
{
    FILE *f;

    errno = 0;

    if ( file == NULL )
        file = "\0"; // fopen() may choke on a NULL file

    if ( ( f = fopen ( file, format ) ) == NULL )
    {
        if ( line == -1 )
            fprintf ( stderr, "Cannot open %s for %s: %s\n\n", file, format, strerror ( errno ) );
        else
            fprintf ( stderr, "%s line:%d  Cannot open %s for %s: %s\n\n",
                      cppfile, line, file, format, strerror ( errno ) );

        exit ( 1 );
    }
    return f;
}
