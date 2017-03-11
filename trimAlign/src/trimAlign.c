/*
  Trim multiple sequence alignment in fasta file format

  author: Pawel Gajer

  date: Sept 26, 2012

*/

#include <getopt.h>
#include <time.h>

#include "CUtilities.h"
#include "IOUtilities.h"

#define DEBUG 0

//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
  fprintf(stderr,
	  "\nUSAGE: %s [Options]\n\n\n"
        "\tOptions:\n"
        "\t-i, --in-file <inFile>\n"
        "\t\tset the input fasta multiple sequence alignment file to <inFile> [mandatory option].\n"
        "\t-o, --out-file <outFile>\n"
        "\t\tset output fasta multiple sequence alignment file to <outFile> [mandatory option]\n"
        "\t-s, --start <startPos>\n"
        "\t\tset the start position of the trimming to <startPos> (0-based) [mandatory option]\n"
        "\t\tThis is the first position of sequence that are going to survive after trimming.\n"
        "\t-e, --end <endPos>\n"
        "\t\tset the end position of the trimming to <endPos> (0-based) [mandatory option]\n"
        "\t\tThis is the last position of sequence that are going to survive after trimming.\n"
        "\t-v, --verbose\n"
        "\t\tverbose mode\n"
        "\t-h, --help\n"
        "\t\tthis message\n\n"

        "\n\tExample: \n"
        "\tTo be executed from the top trimAlign directory\n"
        "\t %s -i data/test1.fa -o data/test1_trim_35_53.fa -s 35 -e 53\n\n"
	  ,s, s);
}

typedef struct
{
  char *inFile;   /// input MSA fasta formated file
  char *outFile;  /// output MSA fasta formated file
  int start;      /// trimming start position (1-based)
  int end;        /// trimming start position (1-based)
  int minLen;     /// minimal sequence length
  int verbose;

} inPar_t;


void printHelp( const char *s )
{
  printUsage(s);
}

void parseArgs( int argc, char ** argv, inPar_t *p );
int trimAlign(inPar_t *inPar);
//timespec diff(timespec start, timespec end);


// ============================== main ======================================
int main(int argc, char **argv)
{
  time_t start = time(NULL);
  //timespec timeStart, timeEnd;
  // CLOCK_REALTIME
  // CLOCK_PROCESS_CPUTIME_ID
  // CLOCK_THREAD_CPUTIME_ID
  //clock_gettime(CLOCK_REALTIME, &timeStart);

  //-- setting up init parameters
  inPar_t *inPar =  (inPar_t *)malloc(sizeof(inPar_t));

  inPar->inFile  = NULL;
  inPar->outFile = NULL;
  inPar->start   = -1;
  inPar->end     = -1;
  inPar->minLen  = -1;
  inPar->verbose = 0;

  //-- parsing input parameters
  parseArgs(argc, argv, inPar);

  if ( !inPar->inFile )
  {
    fprintf(stderr,"Input file is required!\n");
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( !inPar->outFile )
  {
    fprintf(stderr,"Output file has to be specified!\n");
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( !inPar->start )
  {
    fprintf(stderr,"Trimming start position is required!\n");
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( inPar->start < 0 )
  {
    fprintf(stderr,"Trimming start has to be a non-negative integer. Current value of start=%d\n", inPar->start);
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( !inPar->end )
  {
    fprintf(stderr,"Trimming end position is required!\n");
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( inPar->end < 0 )
  {
    fprintf(stderr,"Trimming end has to be a non-negative integer. Current value of end=%d\n", inPar->end);
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( inPar->start >= inPar->end )
  {
    printf("Trim start position has to be less than trim end position.  start=%d\tend=%d\n",inPar->start,inPar->end);
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( !inPar->minLen )
  {
    fprintf(stderr,"Min seq's length is required!\n");
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }
  else if ( inPar->minLen < 0 )
  {
    fprintf(stderr,"Min seq's length needs to be > 0!\n");
    printHelp(argv[0]);
    exit(EXIT_FAILURE);
  }


#if DEBUG
  fprintf(stderr,"inFile:%s\tnoutFile:%s\nstart:%d\nend:%d\nverbose=%d\nminLen:%d\n\n",
	  inPar->inFile, inPar->outFile, inPar->start, inPar->end, inPar->verbose, inPar->minLen);
#endif

  int nSeqs = trimAlign(inPar);

  fprintf(stderr,"\rTrimmed alignment written to %s\n", inPar->outFile);
  fprintf(stderr,"\rProcessed %d sequences\n", nSeqs);

  // report elapsed time
  time_t end = time(NULL);
  long runTime = (long)(end - start);

  //double runTime = difftime(end, start); //(double)(end - start) / CLOCKS_PER_SEC;

  //clock_gettime(CLOCK_REALTIME, &timeEnd);
  //time_t runTime = diff(timeStart,timeEnd).tv_sec; //long runTime.ns = diff(time1,time2).tv_nsec

  if ( runTime > 60 )
  {
    long timeMin = runTime / 60;
    long timeSec = (long)runTime % 60;
    fprintf(stderr, "Elapsed time: %ldmin %ldsec     \n", timeMin, timeSec);
  }
  else
  {
    fprintf(stderr, "Elapsed time: %ldsec           \n", runTime);
  }

  return 0;
}


//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
  int c, errflg = 0;
  optarg = NULL;

  static struct option longOptions[] = {
    {"in-file"     ,required_argument, 0, 'i'},
    {"out-file"    ,required_argument, 0, 'o'},
    {"min-seq-len" ,required_argument, 0, 'l'},
    {"start"       ,required_argument, 0, 's'},
    {"end"         ,required_argument, 0, 'e'},
    {"verbose"     ,no_argument,       0, 0},
    {"help"        ,no_argument,       0, 0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"l:i:o:s:e:vh",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'i':
	p->inFile = strdup(optarg);;
	break;

      case 'o':
	p->outFile = strdup(optarg);;
	break;

      case 's':
	p->start = atoi(optarg);
	break;

      case 'e':
	p->end = atoi(optarg);
	break;

      case 'l':
	p->minLen = atoi(optarg);
	break;

      case 'v':
	p->verbose = 1;
	break;

      case 'h':
	printHelp(argv[0]);
	exit (EXIT_SUCCESS);
	break;

      case 0:
	break;

      default:
	fprintf(stderr,
		"\n=========================================\n"
		" ERROR: Unrecognized option %c\n\n",(char)c);

	int i;
	for ( i=0; i < argc; ++i )
	  fprintf(stderr," %s\n", argv[i]);

	fprintf(stderr,
		"\n==========================================\n\n");
	++errflg;
	break;
    }

  if ( errflg )
  {
    printUsage(argv[0]);
    fprintf(stderr,"Try %s -h for more information\n", argv[0]);
    exit (EXIT_FAILURE);
  }
}


//------------------------------------------------------- trimAlign ----
int trimAlign( inPar_t *p )
{
  MSArec_t *rp = openFasta(p->inFile);
  FILE *out    = fOpen(p->outFile, "w");

  fprintf(stderr,"Trimming %s\n",p->inFile);
  //int nSeqs = numRecordsInFasta( p->inFile );
  //fprintf(stderr,"nSeqs: %d\n", nSeqs);

  int i = 0;
  while ( readMSArecord( rp ) )
  {
    if ( i % 100 == 0 )
      fprintf(stderr,"\rProcessed %d sequences", i);
    //fprintf(stderr,"\rseq: %d [%.2f%%]", i, 100 * (double)i/(double)nSeqs);

    i++;

#if DEBUG
    fprintf(stderr,"seqId: %s\nseq: %s\n", rp->name, rp->seq);
#endif

    if ( rp->seqlen < p->end )
    {
      fprintf(stderr,"Trimming end position is beyond the end of alignment!\n");
      exit(EXIT_FAILURE);
    }

    char *tseq = rp->seq;
    tseq[p->end+1] = '\0';
    tseq += p->start;

    char * s;
    for (s = tseq; *s != '\0'; s++)
    {
      if ( !isgap(*s) )
      {
	rp->ngSeqLen++;
      }
    }

#if DEBUG
    fprintf(stderr,"\nngSeqLen: %d\ttseq: %s\n", rp->ngSeqLen, tseq);
#endif

    if ( rp->ngSeqLen > p->minLen )
    {
      fprintf(out, ">%s\n", rp->name);
      fprintf(out, "%s\n", tseq);
    }
  }

  closeFasta( rp );

  return i;
}
