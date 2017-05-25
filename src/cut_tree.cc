//
// cut_tree.cc
//
// Given a tree and a table of leaf sizes, identify the largest subtrees of size
// (over all their children) less than specified threshold.
//
// Pawel Gajer
// Wednesday, May 24, 2017
//

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <deque>
#include <queue>

using namespace std;

#include "Newick.hh"
#include "IOCppUtilities.hh"
#include "IOCUtilities.h"
#include "CUtilities.h"

#define LINE_LEN  10000

// ===================================================================
//                   subroutine declarations
// ===================================================================

//================================================= inPar_t ====
//! holds input parameters
class inPar_t
{
public:
  inPar_t();
  ~inPar_t();

  char *outDir;        /// input directory for MC model files
  char *treeFile;      /// reference tree file
  char *leafSizeFile;  /// file containg 2 column table with leaf IDs and their sizes
  int cltrMaxSize;     /// upper limit for cluster size
  bool verbose;
  int debug;

  void print();
};


void parseArgs( int argc, char ** argv, inPar_t *p );
void printUsage( const char *s );
void printHelp( const char *s );
void read_char_int_tbl( const char *inputFile, map<string, int> &tbl);
int post_ord_trvl( NewickNode_t *node,
		   map<string, int> &sizeMap,
		   inPar_t *inPar,
		   map<NewickNode_t *, int> &cutMap);


// ===================================================================
//                             main
// ===================================================================

int main(int argc, char **argv)
{
  //-- setting up init parameters
  inPar_t *inPar = new inPar_t();

  //-- parsing input parameters
  parseArgs(argc, argv, inPar);

  if ( !inPar->treeFile )
  {
    fprintf(stderr, "\n\ntERROR in %s at line %d: Missing tree file.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  if ( inPar->outDir )
  {
    string cmd("mkdir -p ");
    cmd += string(inPar->outDir);
    system(cmd.c_str());
  }
  else
  {
    fprintf(stderr, "\n\n\tERROR in %s at line %d: Missing output directory.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  if ( !inPar->leafSizeFile )
  {
    fprintf(stderr, "\n\n\tERROR in %s at line %d: Missing leaf size file.\n\n", __FILE__, __LINE__);
    printHelp(argv[0]);
    exit(1);
  }

  if ( inPar->debug )
    inPar->print();

  fprintf(stderr,"--- Reading tree\n");
  NewickTree_t nt;
  if ( inPar->treeFile ) // load ref tree
  {
    if ( !nt.loadTree(inPar->treeFile) )
    {
      fprintf(stderr,"\n\n\tERROR: Could not load Newick tree from %s\n\n", inPar->treeFile);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr,"\n\n\tERROR: Tree file missing\n\n");
    exit(EXIT_FAILURE);
  }

  int depth = nt.getDepth();
  if ( inPar->verbose )
  {
    fprintf(stderr,"Depth of the tree: %d\n", depth);

    bool withIdx = true;
    nt.printTree( withIdx );
    fprintf(stderr,"\n\n");
  }

  fprintf(stderr,"--- Reading leaf size table\n");
  // ToDo: find out a way to read table without the header. readTable() expects
  // the file to have a header with column names
  double **sizeTbl;
  int nrow, ncol;
  char **rowNames;
  char **colNames;
  readTable( inPar->leafSizeFile, &sizeTbl, &nrow, &ncol, &rowNames, &colNames );

  map<string, int> sizeMap;
  for ( int i = 0; i < nrow; i++ )
    sizeMap[string(rowNames[i])] = sizeTbl[i][0];

  if ( inPar->debug )
  {
    //fprintf(stderr,"\nnrow: %d  ncol: %d\n", nrow, ncol);
    fprintf(stderr,"\nsizeTbl\n");
     for ( int j = 0; j < nrow; ++j )
       fprintf( stderr,"%s\t%d\n", rowNames[j], (int)sizeTbl[j][0] );
    fprintf(stderr,"\n");

    fprintf(stderr,"\nsizeMap\n");
    map<string, int>::iterator it;
    for ( it = sizeMap.begin(); it != sizeMap.end(); it++ )
    {
      fprintf(stderr,"%s\t%d\n", (it->first).c_str(), it->second);
    }
    fprintf(stderr,"\n\n");
  }

  fprintf(stderr,"--- Checking for leaves exceeding cltrMaxSize\n");
  map<string, int> tooBigLeaves;   // those are leaves with size > cltrMaxSize
  map<string, int>::iterator it;
  for ( it = sizeMap.begin(); it != sizeMap.end(); it++ )
  {
    if ( it->second > inPar->cltrMaxSize )
    {
      tooBigLeaves[ it->first ] = it->second;
    }
  }

  if ( !tooBigLeaves.empty() )
  {
    fprintf(stderr,"\n\n\tWARNING: The sizes of the following leaves exceeding cltrMaxSize\n");
    map<string, int>::iterator it;
    for ( it = tooBigLeaves.begin(); it != tooBigLeaves.end(); it++ )
    {
      fprintf(stderr,"%s\t%d\n", (it->first).c_str(), it->second);
    }
    fprintf(stderr,"\n");
  }


  fprintf(stderr,"--- Post-order traversal of the tree with update of internal node sizes\n");
  map<NewickNode_t *, int> cutMap; // keys of this map are nodes corresponding to
				   // maximal size subtrees satisfying size(node)
				   // <= cltrMaxSize
  int totalSize = post_ord_trvl( nt.root(), sizeMap, inPar, cutMap );

  if ( inPar->verbose )
  {
    fprintf(stderr,"\nTotal size: %d\n", totalSize);

    fprintf(stderr,"\nUpdated sizeMap\n");
    map<string, int>::iterator it;
    for ( it = sizeMap.begin(); it != sizeMap.end(); it++ )
    {
      fprintf(stderr,"%s\t%d\n", (it->first).c_str(), it->second);
    }
    fprintf(stderr,"\n");

    fprintf(stderr, "\ncut map\n");
    map<NewickNode_t*, int>::iterator it2;
    for ( it2 = cutMap.begin(); it2 != cutMap.end(); it2++ )
    {
      fprintf(stderr,"%s\t%d\n", (it2->first->label).c_str(), it2->second);
    }
    fprintf(stderr,"\n\n");
  }


  // Generating a nice plot of the tree marking individual clusters

  fprintf(stderr, "--- Generating a cluster table: <leafName> <cluterID= -(node->idx)>\n");

  string outFile = string(inPar->outDir) + string("/") + string("cltr_tbl.txt");
  FILE *out = fOpen( outFile.c_str(), "w" );
  map<NewickNode_t*, int>::iterator it2;
  for ( it2 = cutMap.begin(); it2 != cutMap.end(); it2++ )
  {
    vector<string> leaves;
    nt.leafLabels( it2->first, leaves );

    vector<string>::iterator sItr;
    for ( sItr = leaves.begin(); sItr != leaves.end(); sItr++ )
      fprintf( out,"%s\t%d\n", (*sItr).c_str(), (int)fabs(it2->first->idx));
  }
  fclose( out );

  fprintf( stderr, "\n\n\tCluter table written to %s\n\n", outFile.c_str() );

  return EXIT_SUCCESS;
}


// ===================================================================
//                   subroutine definitions
// ===================================================================

//----------------------------------------------- post_ord_trvl ----
// post order traversal of a tree
int post_ord_trvl( NewickNode_t *node,
		   map<string, int> &sizeMap,
		   inPar_t *inPar,
		   map<NewickNode_t *, int> &cutMap )
{
  int n = 0;
  int numChildren = node->children_m.size();
  if ( numChildren==0 ) // leaf
  {
    n = sizeMap[ node->label ];
  }
  else
  {
    for (int i = 0; i < numChildren; i++)
      n += post_ord_trvl( node->children_m[i], sizeMap, inPar, cutMap );

    if ( node->label.empty() )
    {
      char intStr[16];
      sprintf( intStr, "%d", node->idx );
      node->label = string( intStr ); // "node_" +
    }
    sizeMap[ node->label ] = n;
  }

  if ( n <= inPar->cltrMaxSize )
  {
    cutMap[node] = n;
    for (int i = 0; i < numChildren; i++)
    {
      cutMap.erase( node->children_m[i] );
    }
  }

  return n;
}

//------------------------------------------------------- print ----
void inPar_t::print()
{
  cerr << "\n\nInput parameters" << endl;

  cerr << "outDir: ";
  if ( outDir )
    cerr << outDir << endl;

  cerr << "treeFile: ";
  if ( treeFile )
    cerr << treeFile << endl;

  cerr << "leafSizeFile: ";
  if ( leafSizeFile )
    cerr << leafSizeFile << endl;

  cerr << "Cluster size threshold: " << cltrMaxSize << endl << endl;
}

//------------------------------------------------- constructor ----
inPar_t::inPar_t()
{
  outDir       = NULL;
  treeFile     = NULL;
  leafSizeFile = NULL;
  cltrMaxSize  = 3000;
  verbose      = false;
  debug        = 0;
}

//------------------------------------------------- constructor ----
inPar_t::~inPar_t()
{
  if ( outDir )
    free(outDir);

  if ( treeFile )
    free(treeFile);

  if ( leafSizeFile )
    free(leafSizeFile);
}

//----------------------------------------------------------- printUsage ----
void printUsage( const char *s )
{
  cout << endl

       << "USAGE " << endl
       << endl
       << "  Generating for each non-leaf node of the model tree and all ref sequnces of the\n"
       << "  children nodes a table with each seq's pp w/r each child node.\n"
       << endl
       << "  " << s << " -i <tree file> -s <leaf size file> -o <output dir> [Options]\n"
       << endl
       << "\tOptions:\n"
       << "\t-i <file>  - tree Newick format file\n"
       << "\t-s <file>  - two column table of leaf sizes: <leafName> <size>\n"
       << "\t-o <dir>   - output directory\n"
       << "\t-t <thld>  - cluster size threshold\n"

       << endl
       << "  The following files are generated:\n"
       << "  - ???\n"
       << endl

       << "\n\tExample: \n"

       << "\t" << s << " -v -i ex1.tree -s ex1_leaf_size.txt -o test_cut_tree_dir" << endl << endl;
}


//----------------------------------------------- printHelp ----
void printHelp( const char *s )
{
    printUsage(s);
}


//----------------------------------------------------------- parseArgs ----
//! parse command line arguments
void parseArgs( int argc, char ** argv, inPar_t *p )
{
  int c, errflg = 0;
  optarg = NULL;

  static struct option longOptions[] = {
    {"out-dir"        ,required_argument, 0,          'o'},
    {"tree-file"      ,required_argument, 0,          'i'},
    {"leaf-size-file" ,required_argument, 0,          's'},
    {"cltr-size-thld" ,required_argument, 0,          't'},
    {"help"           ,no_argument,       0,            0},
    {"debug"          ,no_argument, &p->debug,          0},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv,"i:o:s:t:hv",longOptions, NULL)) != -1)
    switch (c)
    {
      case 'i':
	p->treeFile = strdup(optarg);
	break;

      case 's':
	p->leafSizeFile = strdup(optarg);
	break;

      case 'o':
	p->outDir = strdup(optarg);
	break;

      case 't':
	p->cltrMaxSize = atoi(optarg);
	break;

      case 'h':
	printHelp(argv[0]);
	exit (EXIT_SUCCESS);
	break;

      case 'v':
	p->verbose = true;
	break;


      case 0:
	break;

      default:
	cerr << "\n"
	     << "=========================================\n"
	     << " ERROR: Unrecognized option " << (char)c << "\n" << endl;

	for ( int i=0; i < argc; ++i )
	  cerr << argv[i] << " ";
	cerr << endl;

	cerr << "==========================================\n" << endl;
	++errflg;
	break;
    }

  if ( errflg )
  {
    printUsage(argv[0]);
    cerr << "Try '" << argv[0] << " -h' for more information" << endl;
    exit (EXIT_FAILURE);
  }
}
