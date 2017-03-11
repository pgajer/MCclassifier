/*
  Remove gaps from sequences of fasta file

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
	  "\t\tset the input fasta file to <inFile> [mandatory option].\n"
	  "\t-o, --out-file <outFile>\n"
	  "\t\tset output fasta file to <outFile> [mandatory option]\n"
	  "\t\tOutput file contains sequences without gaps"
	  "\t-v, --verbose\n"
	  "\t-q, --quiet\n"
	  "\t\tverbose mode\n"
	  "\t-h, --help\n"
	  "\t\tthis message\n\n"

	  "\n\tExample: \n"
	  "\tTo be executed from the top %s directory\n"
	  "\t %s -i data/test1.fa -o data/test1_gapless.fa\n\n"
	  ,s, s, s);
}

typedef struct
{
  char *inFile;   /// input MSA fasta formated file
  char *outFile;  /// output MSA fasta formated file
  int verbose;
  int quiet;

} inPar_t;


void printHelp( const char *s )
{
  printUsage(s);
}

void parseArgs( int argc, char ** argv, inPar_t *p );
void rmGaps(inPar_t *inPar);



// ============================== main ======================================
int main(int argc, char **argv)
{
  time_t start = time(NULL);

  //-- setting up init parameters
  inPar_t *inPar =  (inPar_t *)malloc(sizeof(inPar_t));

  inPar->inFile  = NULL;
  inPar->outFile = NULL;
  inPar->verbose = 0;
  inPar->quiet = 0;

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

#if DEBUG
  fprintf(stderr,"inFile:%s\tnoutFile:%s\nverbose=%d\n\n",
	  inPar->inFile, inPar->outFile, inPar->verbose);
#endif

  rmGaps(inPar);


  // report elapsed time
  if ( !inPar->quiet )
  {
    time_t end = time(NULL);
    long runTime = (long)(end - start);

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
    {"in-file"  ,required_argument, 0, 'i'},
    {"out-file" ,required_argument, 0, 'o'},
    {"verbose"  ,no_argument,       0, 'v'},
    {"quiet"    ,no_argument,       0, 'q'},
    {"help"     ,no_argument,       0, 0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"i:o:qvh",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'i':
	p->inFile = strdup(optarg);;
	break;

      case 'o':
	p->outFile = strdup(optarg);;
	break;

      case 'v':
	p->verbose = 1;
	break;

      case 'q':
	p->quiet = 1;
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


//------------------------------------------------------- rmGaps ----
void rmGaps( inPar_t *p )
{
  FASTArec_t *rp = openFasta(p->inFile);
  FILE *out    = fOpen(p->outFile, "w");

  if ( !p->quiet )
    fprintf(stderr,"Removing gaps from %s\n",p->inFile);
  //int nSeqs = numRecordsInFasta( p->inFile );
  //fprintf(stderr,"nSeqs: %d\n", nSeqs);

  int i = 0;
  while ( readFASTArecord( rp ) )
  {
    if ( !p->quiet && (i % 100 == 0) )
      fprintf(stderr,"\rProcessed %d sequences", i);
    //fprintf(stderr,"\rseq: %d [%.2f%%]", i, 100 * (double)i/(double)nSeqs);

    i++;

#if DEBUG
    fprintf(stderr,"seqId: %s\nseq: %s\n", rp->name, rp->seq);
#endif

    fprintf(out, ">%s\n", rp->name);
    fprintf(out, "%s\n", rp->seq);
  }

  closeFasta( rp );

  if ( !p->quiet )
  {
    fprintf(stderr,"\rOutput file written to %s\n", p->outFile);
    fprintf(stderr,"\rProcessed %d sequences\n", i);
  }
}
