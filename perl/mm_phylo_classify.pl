#!/usr/bin/env perl

=head1 NAME

  mm_phylo_classify.pl

=head1 DESCRIPTION

  Given data generated by mm_validate_pecan.pl, generate taxonomic assignment to
  query sequences of the M&M project.

  Here are the classification rules.

  -  If the cluster consists of query sequences and reference sequences of one
  species, all query sequences of that cluster are assigned to that species.

  -  If the cluster consists of only query sequences, then

      A)  If the number of non-redundant query sequences is greater than 50, its
  going to be assigned to a new Genus\_sp\_index species, where index will be a
  new index if Genus\_sp\_index type species already exist.
      B)  Otherwise (the number of non-redundant query sequences is less than or
  equal to 50), query sequences are removed from the dataset.

  -  If query sequences are present in two or more clusters, each containing
  only query sequences and the same species, spA, then if the number of
  non-redundant query sequences in these clusters is large enough (???), form
  different versions of spA: spA\_1, ... , spA\_n. Otherwise, assign sequences of
  all these clusters to spA.

  - If query sequences are in a cluster with more than one species, then use the
  majority vote to assign to the query sequences the taxonomy of the most
  abundant species. Change the taxonomy of less abundant (<=10) species to the
  most abundant species.


=head1 SYNOPSIS

  mm_phylo_classify.pl

=head1 OPTIONS

=over

=item B<--val-dir, -d>
  Directory containing validation directories generated by mm_validate_pecan.pl.

=item B<--out-dir, -o>
  Output directory.

=item B<--val-dirs-file, -i>
  File with a list of mm_validate_reports_dir directories to be processed.

=item B<--cltr-min-nr-seqs>
  Minimal number of non-redundant query seq's in a cluster with only query
  sequences or a cluster with more than one species. These are the cases were the
  query sequences will be dropped if their number is below this threshold. Default value: 50.

=item B<--cltr-max-nr-seqs>
  Maximal number of non-redundant query seq's in a cluster to be added to a
  training set. Default: 500.

=item B<--report-spp-in-mult-cltrs>
  Print to stdout cluster info for species present in more than 1 cluster with 10
  or more seq's in 2 or more of them.

=item B<--sp-mult-thld>
    We want at least two cluster to have at least spMultThld seq's of that species.

=item B<--min-sp-size>
  Species in a cluster with multiple species that have $minSpSize or less seq's
  are renamed to the taxonomy of the most abundant species.

=item B<--quiet>
  Do not print progress messages.

=item B<--show-tree>
  Open the pdf file with the tree used to do clustering.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  in ~/projects/M_and_M/new_16S_classification

  mm_phylo_classify.pl --debug -d mm_may12_validate_pecan_dir2 -o mm_may12_phylo_clfy_dir2

  mm_phylo_classify.pl --quiet -o mm_phylo_clfy_dir

  mm_phylo_classify.pl --debug --valid-dirs-file valid_dirs.txt --cltr-min-nr-seqs 50 -o phylo_clfy_dir

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd 'abs_path';
use List::Util qw( sum min max );
use File::Temp qw/ tempfile /;
#use Parallel::ForkManager;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $spMultThld = 10; # we want at least two cluster to have at least spMultThld seq's of that species
my $minSpSize  = 10; # species in a cluster with multiple species that have
		     # $minSpSize or less seq's are renamed to the taxonomy of
		     # the most abundant species
my $minNRQseqs = 50;
my $maxNRQseqs = 500;

GetOptions(
  "val-dirs-file|i=s"        => \my $vDirsFile,
  "val-dir|d=s"              => \my $valDir,
  "out-dir|o=s"              => \my $outDir,
  "cltr-min-nr-seqs"         => \$minNRQseqs,
  "cltr-max-nr-seqs"         => \$maxNRQseqs,
  "report-spp-in-mult-cltrs" => \my $reportSppInMultCltrs,
  "sp-mult-thld"             => \$spMultThld,
  "min-sp-size"              => \$minSpSize,
  "quiet"                    => \my $quiet,
  "igs"                      => \my $igs,
  "show-tree"                => \my $showTree,
  "verbose|v"                => \my $verbose,
  "debug"                    => \my $debug,
  "dry-run"                  => \my $dryRun,
  "help|h!"                  => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $baseDir        = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";
#my $mmDir          = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/";
my $mmDir          = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mm_may12_validate_pecan_dir";
my $spToPhGrFile   = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/MM_may12_spp_phGr_tbl.txt";

my $R              = "R";
my $mothur         = "/Users/pgajer/bin/mothur";
my $usearch6       = "/Users/pgajer/bin/usearch6.0.203_i86osx32";
my $readNewickFile = "/Users/pgajer/.Rlocal/read.newick.R";

if ( defined $igs )
{
  $R               = "/home/pgajer/bin/R";
  $mothur          = "/usr/local/packages/mothur-1.36.1/mothur";
  $usearch6        = "/local/projects/pgajer/bin/usearch6.0.203_i86linux32";
  $readNewickFile  = "/local/projects/pgajer/devel/MCclassifier/perl/read.newick.R";
  $spToPhGrFile    = "/local/scratch/MM/MM_may12_spp_phGr_tbl.txt";
}

if ( !$outDir )
{
  print "\n\n\tERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $debugStr = "";
my $quietStr = "--quiet";
if ($debug)
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ($verbose)
{
  $verboseStr = "--verbose";
}


####################################################################
##                               MAIN
####################################################################

my $startRun = time();
my $endRun = time();
my $runTime = $endRun - $startRun;
my $timeStr;
my $timeMin = int($runTime / 60);
my $timeSec = $runTime % 60;

# Get the number of species found by the classifier in the M&M dataset
#my $spToPhGrFile = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mmDir_May5/sp_to_phGr.txt";
my $wcline = qx/ wc -l $spToPhGrFile /;
$wcline =~ s/^\s+//;
my ($nAllSpp, $str) = split /\s+/, $wcline;

my %spIdx = getSpIdx( $spToPhGrFile ); # genus => number of _sp species within that genus

# my $tmpSpIdxFile = "tmpSpIdxFile.txt";
# writeTbl( \%spIdx, $tmpSpIdxFile );
# exit;

my @vDirs;
if ( $vDirsFile )
{
  @vDirs = read_array( $vDirsFile );
}
elsif ( $valDir )
{
  push @vDirs, $valDir;
  # opendir VDIR, $valDir or die "Error in opening $valDir: $OS_ERROR";
  # while( defined (my $dir = readdir(VDIR)))
  # {
  #   next if $dir =~ /^..?$/; # skip . and ..
  #   next if $dir !~ /^mm_*_dir$/;
  #   push @vDirs, $dir;
  # }
}
else
{
  ##
  ## cd /Users/pgajer/projects/M_and_M/new_16S_classification_data
  ## find . -name "mm_validate_reports_dir_*" -type d -maxdepth 1 -mindepth 1
  ##
  @vDirs = ("mm_validate_reports_dir_20170505_160208",
	    "mm_validate_reports_dir_20170505_163202",
	    "mm_validate_reports_dir_20170505_193808",
	    "mm_validate_reports_dir_20170507_032747",
	    "mm_validate_reports_dir_20170507_074358",
	    "mm_validate_reports_dir_20170507_105808",
	    "mm_validate_reports_dir_20170507_180210",
	    "mm_validate_reports_dir_20170507_181736",
	    "mm_validate_reports_dir_20170507_205318");

  my $dir = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/";
  @vDirs = map { $dir . $_ } @vDirs;
}

if ( $debug )
{
  print_array( \@vDirs, "vDirs" );
  print "\n";
  #exit;
}


if ( ! -e $outDir )
{
  my $cmd = "mkdir -p $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?";# if !$dryRun;
}

my $tmpDir = $outDir . "/temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my %phyloTx;          # phylogeny based taxonomy; seqID => taxonomic assignment to seqID

my %rmSeqIDs;         # seqID => phylo-group of the sequence

my %ipSpp;;           #  sp => phylo-group of sp; table of incompletely processed
		      #  species. The table is defined on species for which vicut run
 		      #  data was not found.

my %gr1sppCltrs;      # sp => phylo-group of sp, for species with cluster(s)
		      # containing more than one species; Generate for them combined
		      # report and maybe also treeDir with links to their tree files

my %processed;        # sp => 1 upon complection of processing the given species

my $nMultCltrTxs = 0; # number multiple-cluster taxons (species present in more
		      # than one cluster containing also query sequences).
my $showCtrs = 0;


my $queryTxFile = $outDir . "/query.tx"; # 2 columns; seqID, species; tab
					 # delimited file of all query sequences
					 # (redundant) except those who were
					 # dropped due to presence in a cluster
					 # with no named species and not being
					 # abundant enough.

my $trainTxFile = $outDir . "/train.tx"; # columns: seqID, species, phyloGroup, multiplicity
                                         # multiplicity = size of the 100%
                                         # identity cluster for query seq's, -1
                                         # for ref seq's.  ref and query
                                         # (non-redundant sequences only) are
                                         # listed here. If more than 500 query
                                         # sequences are present in a given
                                         # cluster classified to a species, only
                                         # the first 500 most abundant are used.

my $chgdRefTxFile = $outDir . "/changed.tx"; # if there are <= 10 seq's of a
					     # species in a cluster with multiple
					     # species, its taxonomy is changed
					     # to that of the most abundant
					     # species.

my $droppedSeqs = $outDir . "/dropped.seqIDs"; # seqIDs of query sequences
					       # dropped from the classification
					       # according to the above rules
my $newSpFile = $outDir . "/new_spp.txt"; #

my @droppedQuerySeqs;

# my $trainTxFile = $outDir . "/train_weighted.tx"; # weighted version of the
# trainin set were non-redundant sequences would be in different aboundances
# based on their presence in the population

open QOUT, ">$queryTxFile" or die "Cannot open $queryTxFile for writing: $OS_ERROR";
open TOUT, ">$trainTxFile" or die "Cannot open $trainTxFile for writing: $OS_ERROR";
open ROUT, ">$chgdRefTxFile" or die "Cannot open $chgdRefTxFile for writing: $OS_ERROR";
open NSPOUT, ">$newSpFile" or die "Cannot open $newSpFile for writing: $OS_ERROR";

my $spCounter = 1;
for my $vDir ( @vDirs )
{
  #print "vDir: $vDir\n";
  opendir VDIR, $vDir or die "Error in opening $vDir: $OS_ERROR";
  while( defined (my $phGrDir = readdir(VDIR)))
  {
    next if $phGrDir =~ /^..?$/;     # skip . and ..
    print "phGrDir: $phGrDir\n";
    if ( $phGrDir =~ /^mm_/ && -d "$vDir/$phGrDir" )
    {
      opendir PGDIR, "$vDir/$phGrDir" or die "Error in opening $vDir/$phGrDir: $OS_ERROR";
      while( defined (my $spDir = readdir(PGDIR)))
      {
	next if $spDir =~ /^..?$/;     # skip . and ..
	print "spDir: $spDir\n";
	if ( $spDir =~ /_dir$/ && -d "$vDir/$phGrDir/$spDir" )
	{
	  #print "\nDiscovered $vDir/$phGrDir/$spDir\n";
	  #print "\nDiscovered $spDir\n";

	  if ( $quiet )
	  {
	    printf( "\r%d [%.1f%%]", $spCounter, 100 * $spCounter / $nAllSpp );
	    $spCounter++;
	  }

	  ## extracting species name
	  my ( $sp ) = ( $spDir =~ /(\w+)_dir$/ );
	  #print "sp: $sp\n";

	  if ( !$sp )
	  {
	    warn "\n\n\tERROR Undefined sp from $spDir";
	    print "\n\n";
	    exit;
	  }
	  my ( $phGr ) = ( $phGrDir =~ /mm_(\w+)_dir/ );
	  #print "phGr: $phGr\n";

	  next if exists $processed{$sp};

	  print "\n-------------------------------------------------------------------------------------\n$sp  ($phGr)\n" if !$quiet;

	  ## There are two vicut directory names possible here
	  my $covSuffix = "_nr_cov80";
	  my $vicutDir80  = "$vDir/$phGrDir/$spDir/$sp" . $covSuffix . "_vicut_dir";
	  $covSuffix = "";
	  my $vicutDir  = "$vDir/$phGrDir/$spDir/$sp" . "_vicut_dir";

	  # Testing which one is present
	  if ( ! -e $vicutDir )
	  {
	    #print "Did not find $vicutDir\n";
	    if ( -e $vicutDir80 )
	    {
	      #print "Found $vicutDir80\n";
	      $covSuffix = "_nr_cov80";
	      $vicutDir = $vicutDir80;
	    }
	    else
	    {
	      #print "WARNING: Did not find neither $covSuffix" . "_vicut_dir nor _vicut_dir for $sp\n" if !$quiet;
	      $ipSpp{$sp} = $phGr;
	      next;
	    }
	  }

	  my $vicutCltrsFile = $vicutDir . "/minNodeCut.cltrs";
	  if ( ! -e $vicutCltrsFile )
	  {
	    #print "WARNING: Did not find $vicutCltrsFile\n" if !$quiet;
	    $ipSpp{$sp} = $phGr;
	    next;
	  }

	  my ($rvCltrTbl, $rvTxTbl, $rvExtTxTbl) = readCltrsTbl($vicutCltrsFile);

	  my %vCltrTbl   = %{$rvCltrTbl};  # seqID => vicut cluster ID
	  my %vTxTbl     = %{$rvTxTbl};    # seqID => taxonomy (NA for query seq's)
	  my %vExtTxTbl  = %{$rvExtTxTbl}; # seqID => taxonomy of seqID if seqID is a phGr ref seq and c<vicut cluster ID of seqID> if seqID is a query seq

	  ## Extended taxonomy table should be alread present
	  my $vExtTxTblFile = "$vDir/$phGrDir/$spDir/$sp" . $covSuffix . "_ext.tx";
	  if ( ! -e $vExtTxTblFile )
	  {
	    #print "WARNING: Did not find $vExtTxTblFile\n" if !$quiet;
	    $ipSpp{$sp} = $phGr;
	    next;
	  }

	  ## vicut-cltr/tx frequency table
	  my %vCltrvTxFreq; # $vCltrvTxFreq{cltr}{tx} = # of seq IDs of taxon tx in the cluster cltr
	  my %txCltrFreq;   # $txCltrFreq{tx}{cltr} = # of seq IDs of taxon tx in the cluster cltr
	  my %vCltrvTxIds;  # $vCltrvTxIds{cltr}{tx} = ref to seqID of the cluster's, cltr, taxon, tx.
	  my %vCltrIds;     # cl => seq IDs within the given cluster

	  for my $id ( keys %vCltrTbl )
	  {
	    $vCltrvTxFreq{$vCltrTbl{$id}}{$vTxTbl{$id}}++;
	    $txCltrFreq{ $vTxTbl{$id} }{ $vCltrTbl{$id} }++;
	    push @{$vCltrvTxIds{$vCltrTbl{$id}}{$vTxTbl{$id}}}, $id;
	    push @{$vCltrIds{$vCltrTbl{$id}}}, $id;
	  }

	  ## Extracting seqIDs of non-redundant query sequences
	  my $nrSeqIDsFile = "$vDir/$phGrDir/$spDir/$sp" . "_nr.seqIDs";
	  my @nrAllSeqIDs = read_NR_array( $nrSeqIDsFile );

	  my @nrSeqIDs;
	  if ( $covSuffix ne "" )
	  {
	    $nrSeqIDsFile = "$vDir/$phGrDir/$spDir/$sp" . $covSuffix . ".seqIDs";
	    @nrSeqIDs = read_array( $nrSeqIDsFile );
	  }
	  else
	  {
	    @nrSeqIDs = @nrAllSeqIDs;
	  }

	  ## Identifing clusters that contain query sequences
	  my @querySeqCltrs;
	  for ( @nrSeqIDs )
	  {
	    if ( exists $vCltrTbl{$_} )
	    {
	      push @querySeqCltrs, $vCltrTbl{$_};
	    }
	    else
	    {
	      warn "\n\n\tERROR: $_ undefined in vCltrTbl: $vicutCltrsFile";
	      print "\n\n";
	      #exit;
	    }
	  }
	  my @queryCltrs = unique(\@querySeqCltrs);

	  print "\nqueryCltrs: @queryCltrs\n" if !$quiet;

	  ## size of each cluster
	  my %vicutCltrSize;
	  for my $cl ( @queryCltrs )
	  {
	    if (exists $vCltrvTxFreq{$cl})
	    {
	      my @txs = keys %{$vCltrvTxFreq{$cl}};
	      my $size = 0;
	      for my $tx (@txs)
	      {
		$size += $vCltrvTxFreq{$cl}{$tx};
	      }
	      $vicutCltrSize{$cl} = $size;
	    }
	    else
	    {
	      warn "\n\nERROR $cl not found in vCltrvTxFreq";
	      print "\n\n";
	      exit;
	    }
	  }

	  my $spClstr2File = "$vDir/$phGrDir/$spDir/$sp" . "_nr.clstr2";
	  #print "--- Parsing clstr2 file\n";
	  my %cTbl = parseClstr2($spClstr2File);

	  ##
	  ## The classification section
	  ##
	  my @na = ("NA");
	  my $winnerTx; # taxonomy of the cluster with the largest number of
			# query sequences. Used to propagate the taxonomy of the
			# non-redundant seq's covering 80% of all seq's to the
	                # remaining 20%.
	  my $maxNumQuerySeqs = 0;
	  for my $cl ( @queryCltrs )
	  {
	    my %txFreq = %{$vCltrvTxFreq{$cl}};
	    my @txs    = keys %txFreq;
	    my $nTxs   = @txs;

	    my @nrQueryIDs = @{ $vCltrvTxIds{$cl}{"NA"} };
	    my $newTx;

	    if ( $nTxs == 1 ) # only NAs
	    {
	      my $tx = shift @txs;
	      if ( $tx ne "NA" )
	      {
		warn "\n\n\tERROR: single taxon cluster (selected from clusters that contain NA) and the 'taxon' is not NA";
		print "\n\n";
		exit 1;
	      }

	      if ( @nrQueryIDs > $minNRQseqs )
	      {
		my ($genus, $s) = split "_", $sp;
		$spIdx{$genus}++;
		$newTx = $genus . "_sp_" . $spIdx{$genus};
		print NSPOUT "Cluster $cl of $sp renamed to $newTx\n";
	      }
	      else
	      {
		for my $refID ( @nrQueryIDs )
		{
		  my @qIDs = @{ $cTbl{$refID} };
		  push @droppedQuerySeqs, @qIDs;
		}
		next;
	      }
	    }
	    elsif ( $nTxs == 2 ) # NAs and a named species
	    {
	      my @nonNAtx = diff( \@txs, \@na );
	      $newTx = shift @nonNAtx;
	    }
	    else
	    {
	      # NAs with more than one species
	      @txs = diff( \@txs, \@na );

	      # sorting txs w/r to their number in the cluster
	      @txs = sort { $txFreq{$b} <=> $txFreq{$a} } @txs;
	      $newTx = shift @txs;

	      ## here <= 10 taxons/species should be renamed to $newTx <= ToDO !!!!!!!
	      for my $tx ( @txs )
	      {
		if ( $tx ne "NA" && $txFreq{$tx} <= $minSpSize )
		{
		  print ROUT "$tx\t$newTx\t$phGr\n";
		}
	      }
	    }

	    for my $refID ( @nrQueryIDs )
	    {
	      print TOUT "$refID\t$newTx\t$phGr\t" . @{ $cTbl{$refID} } . "\n";
	      my @qIDs = @{ $cTbl{$refID} };
	      for ( @qIDs )
	      {
		print QOUT "$_\t$newTx\n";
	      }
	    }

	    if ( @nrQueryIDs > $maxNumQuerySeqs )
	    {
	      $maxNumQuerySeqs = @nrQueryIDs;
	      $winnerTx = $newTx;
	    }

	  } # end of for my $cl ( @queryCltrs )

	  ## Propagating the winner taxonomy to the remaining 20% of sequences if
	  ## we are in the case of 80% coverage situation.
	  if ( $covSuffix ne "")
	  {
	    my @tailIDs = diff( \@nrAllSeqIDs, \@nrSeqIDs );

	    for my $refID ( @tailIDs )
	    {
	      if ( !exists $cTbl{$refID} )
	      {
		warn "\n\n\tERROR: $refID does not exist in cTbl";
		#print_array( );
		exit;
	      }
	      my @qIDs = @{ $cTbl{$refID} };
	      for ( @qIDs )
	      {
		print QOUT "$_\t$winnerTx\n";
	      }
	    }
	  }

	  ## Checking out for species present in more than one cluster
	  if ( $reportSppInMultCltrs )
	  {
	    my @txs = keys %txCltrFreq;
	    @txs = diff( \@txs, \@na );
	    for my $tx ( @txs )
	    {
	      my %cltrMult = %{ $txCltrFreq{$tx} }; # cl => # of seq's of tx in cluster cl
	      my @txCltrs  = keys %cltrMult;
	      my @c = comm( \@txCltrs, \@queryCltrs );
	      if ( @c > 1 )
	      {
		my $nBigCounts = 0;
		for my $cl ( @c )
		{
		  $nBigCounts++ if $cltrMult{$cl} >= $spMultThld;
		}
		if ( $nBigCounts > 1 )
		{
		  print "\n-------------------------------------------------------------------------------------\n";
		  print "$sp  ($phGr)\n";
		  print "$tx found in more than one cluster with query sequences (with >= $spMultThld seq's there)\n";
		  $nMultCltrTxs++;
		  $showCtrs = 1;
		}
	      }
	    }
	  }

	  my @querySortedCltrs = sort { $vicutCltrSize{$b} <=> $vicutCltrSize{$a} } @queryCltrs;
	  for my $cl (@querySortedCltrs)
	  {
	    print "\nCluster $cl (" . $vicutCltrSize{$cl} . ")\n" if !$quiet || $showCtrs;

	    ## Generating a list of species present in $cl sorted by size and with NA
	    ## at the end (igoring the size of NA when sorting
	    my @clTxs = keys %{$vCltrvTxFreq{$cl}};
	    @clTxs = sort { $vCltrvTxFreq{$cl}{$b} <=> $vCltrvTxFreq{$cl}{$a} } @clTxs;
	    ## putting NA at the end
	    @clTxs = diff(\@clTxs, \@na);
	    push @clTxs, "NA";
	    my %clTxsize;
	    for my $tx ( @clTxs )
	    {
	      $clTxsize{$tx} = $vCltrvTxFreq{$cl}{$tx};
	      #print "\t$tx\t" . $vCltrvTxFreq{$cl}{$tx} . "\n";
	    }
	    #print "\n";

	    printFormatedTbl(\%clTxsize, \@clTxs) if !$quiet || $showCtrs;
	  }
	  print "\n" if !$quiet || $showCtrs;
	  $showCtrs = 0;

	  if ( exists $ipSpp{$sp} )
	  {
	    delete $ipSpp{$sp};
	  }

	  $processed{$sp} = 1;
	}
      }
      closedir PGDIR;
    }
  }
  closedir VDIR;
}

close QOUT;
close TOUT;
close ROUT;
close NSPOUT;

##
## Summary stats
##

my $nProcessedSpp = keys %processed;

my @ipSpp = sort { $ipSpp{$a} cmp $ipSpp{$b} } keys %ipSpp;;
my $nIPrSpp = @ipSpp;

print "\n\nNumber of species found by the classifier in the M&M dataset: $nAllSpp\n";
print     "Number processed species:                                     $nProcessedSpp\n";
print     "Number of species for which no vicut data was found:          $nIPrSpp\n";

if ( $reportSppInMultCltrs )
{
  print     "\nNumber of species found in more than one cluster with query sequences  (with >= $spMultThld seq's there): $nMultCltrTxs\n";
}


if ( $nIPrSpp )
{
  print "\nSpecies for which no vicut data was found\n";
  printFormatedTbl( \%ipSpp, \@ipSpp );
  print "\n\n";

  my $ipSppFile = $outDir . "/incompletely_processed_spp.txt";
  open OUT, ">$ipSppFile" or die "Cannot open $ipSppFile for writing: $OS_ERROR\n";
  for ( @ipSpp )
  {
    print OUT "$_\t-1\t" . $ipSpp{$_} . "\n";
  }
  close OUT;
}

if ( $nAllSpp > $nProcessedSpp )
{
  # write to outDir table of missing species
  my %spPhGrTbl = parse_phGr_tbl( $spToPhGrFile );
  my @allSpp = keys %spPhGrTbl;

  my @processedSpp = keys %processed;

  my @d = diff( \@allSpp, \@processedSpp );

  my $outFile = $outDir . "/missing_spp.txt";
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  for ( @d )
  {
    if ( !exists $spPhGrTbl{$_} )
    {
      warn "\n\n\tERROR: $_ does not exist in spPhGrTbl";
      print "\n\n";
      exit;
    }
    print OUT "$_\t-1\t" . $spPhGrTbl{$_} . "\n";
    print "$_\t-1\t" . $spPhGrTbl{$_} . "\n";
  }
  close OUT;
  print "\nTable of missing species written to $outFile\n";
}

if ( @droppedQuerySeqs )
{
  writeArray( \@droppedQuerySeqs, $droppedSeqs );
}

print "Output written to $outDir\n\n";

## print "\n\n\tMissing species table written to $ipSppFile\n\n";

####################################################################
##                               SUBS
####################################################################

sub parse_phGr_tbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($sp, $phGr) = split /\s+/,$_;
    $tbl{$sp} = $phGr;
  }
  close IN;

  return %tbl;
}

# parse a clstr2 file
# output table: refId -> number of elements in the corresponding cluster
sub parseClstr2
{
  my $inFile = shift;

  my %tbl;
  open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
  foreach my $rec (<IN>)
  {
    chomp $rec;
    my @ids = split ",", $rec;
    my $refId = shift @ids;
    push @{ $tbl{$refId} }, @ids;
  }
  close IN;

  return %tbl;
}

##
## parse 3 column tbl
##

## 1642.V1_0	Lactobacillus_iners	0.93
## 0980.V2_1	Lactobacillus_iners	0.97
## 1670.V2_2	Lactobacillus_helveticus_acidophilus	0.98
## 0711.V3_3	Atopobium_vaginae	0.56
## 1149.V1_4	Lactobacillus_iners	0.94
## 1386.V1_5	Prevotella_buccalis	0.85
## 1119.V2_6	BVAB1	0.79
## 1449.V1_7	BVAB1	0.97
## 1600.V1_8	BVAB3	0.93

sub parse_pecan_tbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readQtxTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %spIDsTbl;
  ##my %phGrSppTbl;
  my %ppTbl;
  ##my %sp2phGrSppTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    next if /^$/;
    chomp;
    my ($id, $sp, $pp) = split /\s+/,$_;
    push @{$spIDsTbl{$sp}}, $id;
    $ppTbl{$id} = $pp;
  }
  close IN;

  return (\%spIDsTbl, \%ppTbl);
}


##
## parse 3 column species table
##

# file format

#   Akkermansia_muciniphila	2413	Verrucomicrobia_V3V4
#   Akkermansia_sp_5	14	Verrucomicrobia_V3V4
#   Akkermansia_sp_7	13	Verrucomicrobia_V3V4
#   Coraliomargarita_akajimensis	2	Verrucomicrobia_V3V4
#   Akkermansia_sp_4	1	Verrucomicrobia_V3V4
#   Fibrobacter_sp	51504	phyla_lessthen_1k_wOG_V3V4
#   Jonquetella_anthropi	1299	phyla_lessthen_1k_wOG_V3V4
#   Chlamydia_trachomatis	128	phyla_lessthen_1k_wOG_V3V4
#   Pyramidobacter_piscolens	62	phyla_lessthen_1k_wOG_V3V4
#   Cloacibacillus_evryensis	37	phyla_lessthen_1k_wOG_V3V4

sub parse_spp_tbl
{
  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in parse_spp_tbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %phGrSppTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    next if /^$/;
    chomp;
    my ($sp, $size, $gr) = split /\s+/,$_;
    push @{$phGrSppTbl{$gr}}, $sp;
  }
  close IN;

  return %phGrSppTbl;
}



# read 3 column clstrs table
sub readCltrsTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %vCltrTbl;
  my %txTbl;
  my %txTbl2;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $cl, $tx) = split /\s+/,$_;
    if (defined $id)
    {
      $vCltrTbl{$id} = "c" . $cl;
      $txTbl{$id} = $tx;
      if ($tx ne "NA")
      {
	$txTbl2{$id} = $tx;
      }
      else
      {
	$txTbl2{$id} = "c" . $cl;
      }
    }
  }
  close IN;

  return (\%vCltrTbl, \%txTbl, \%txTbl2);
}

# read 2 or 3 column table; create a table that assigns
# elements of the first column to the second column
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t, $r) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# read lineage table
sub readLineageTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readLineageTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
    ## test for '/' characters
    if ($t =~ /\//)
    {
      warn "\n\n\tERROR: Discovered '/' for id: $id\t$t";
      print "\n\n";
      exit 1;
    }
  }
  close IN;

  return %tbl;
}

sub get_seqIDs_from_fa
{
  my $file = shift;

  my $quiet = 1;
  my $startRun = time();
  my $endRun = time();

  open (IN, "<$file") or die "Cannot open $file for reading: $OS_ERROR";
  $/ = ">";
  my $junkFirstOne = <IN>;
  my $count = 1;
  my $timeStr = "";
  my @seqIDs;
  while (<IN>)
  {
    if ( !$quiet && ($count % 500 == 0) )
    {
      $endRun = time();
      my $runTime = $endRun - $startRun;
      if ( $runTime > 60 )
      {
	my $timeMin = int($runTime / 60);
	my $timeSec = sprintf("%02d", $runTime % 60);
	$timeStr = "$timeMin:$timeSec";
      }
      else
      {
	my $runTime = sprintf("%02d", $runTime);
	$timeStr = "0:$runTime";
      }
      print "\r$timeStr";
    }

    chomp;
    my ($id,@seqlines) = split /\n/, $_;
    push @seqIDs, $id;
    $count++;
  }
  close IN;
  $/ = "\n";

  return @seqIDs;
}

# common part of two arrays
sub comm{

  my ($a1, $a2) = @_;

  my @c; # common array
  my %count;

  foreach my $e (@{$a1}, @{$a2}){ $count{$e}++ }

  foreach my $e (keys %count)
  {
    push @c, $e if $count{$e} == 2;
  }

  return @c;
}

# read table with one column
sub read_array{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in read_array(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  if ( defined $hasHeader )
  {
    <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    push @rows, $_;
  }
  close IN;

  return @rows;
}

# read array of non-redundant seqIDs
sub read_NR_array{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in read_array(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  if ( defined $hasHeader )
  {
    <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    s/;size=\d+;//;
    push @rows, $_;
  }
  close IN;

  return @rows;
}

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}


# difference of two arrays
sub diff{

  my ($a1, $a2) = @_;

  my (%aa1, %aa2);

  foreach my $e (@{$a1}){ $aa1{$e} = 1; }
  foreach my $e (@{$a2}){ $aa2{$e} = 1; }

  my @d; # dfference array

  foreach my $e (keys %aa1, keys %aa2)
  {
    push @d, $e if exists $aa1{$e} && !exists $aa2{$e};
  }

  return @d;
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

# print elements of a hash table
sub printTbl{

  my $rTbl = shift;
  map {print "$_\t" . $rTbl->{$_} . "\n"} keys %$rTbl;
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTbl{

  my ($rTbl, $rSub) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  for (@args)
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print "$_$pad" . $rTbl->{$_} . "\n";
  }
  #print "\n";
}

# print elements of a hash table so that arguments are aligned
sub printFormatedTblToFile{

  my ($rTbl, $rSub, $fh) = @_; # the second argument is a subarray of the keys of the table

  my @args;
  if ($rSub)
  {
    @args = @{$rSub};
  }
  else
  {
    @args = keys %{$rTbl};
  }

  my $maxStrLen = 0;
  map { $maxStrLen = length($_) if( length($_) > $maxStrLen )} @args;

  for (@args)
  {
    my $n = $maxStrLen - length($_);
    my $pad = ": ";
    for (my $i=0; $i<$n; $i++)
    {
      $pad .= " ";
    }
    print $fh "$_$pad" . $rTbl->{$_} . "\n";
  }
  print $fh "\n";
}

# write array to a file (one column format)
sub writeArray
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify {
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}


## plot tree with clade colors
sub plot_tree
{
  my ($treeFile, $pdfFile, $title) = @_;

  my $showBoostrapVals = "T";

  if (!defined $title)
  {
    $title = "";
  }

  my $Rscript = qq~

source(\"$readNewickFile\")
require(phytools)
library(ade4)

tr <- read.newick(file=\"$treeFile\")
tr <- collapse.singles(tr)

(nLeaves <- length(tr\$tip.label))

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr, type=\"phylogram\", no.margin=FALSE, show.node.label=F, cex=0.8, main=\"$title\")
par(op)
dev.off()
~;

  run_R_script( $Rscript );
}

sub plot_tree2
{
  my ($treeFile, $clFile, $pdfFile, $title) = @_;

  my $showBoostrapVals = "F";

  if (!defined $title)
  {
    $title = "";
  }

  my $Rscript = qq~

clTbl <- read.table(\"$clFile\", header=F)
str(clTbl)

cltr <- clTbl[,2]
names(cltr) <- clTbl[,1]

source(\"$readNewickFile\")
require(phytools)

tr1 <- read.newick(file=\"$treeFile\")
tr1 <- collapse.singles(tr1)

tip.cltr <- cltr[tr1\$tip.label]

colIdx <- 1
tip.colors <- c()
tip.colors[1] <- colIdx
for ( i in 2:length(tip.cltr) )
{
    if ( tip.cltr[i] != tip.cltr[i-1] )
    {
        colIdx <- colIdx + 1
        if ( colIdx==9 )
        {
            colIdx <- 1
        }
    }
    tip.colors[i] <- colIdx
    if ( colIdx==7 )
    {
        tip.colors[i] <- "brown"
    }
}

(nLeaves <- length(tr1\$tip.label))

figH <- 8
figW <- 6
if ( nLeaves >= 50 )
{
    figH <- 6.0/50.0 * ( nLeaves - 50) + 10
    figW <- 6.0/50.0 * ( nLeaves - 50) + 6
}

pdf(\"$pdfFile\", width=figW, height=figH)
op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
plot(tr1,type=\"phylogram\", no.margin=FALSE, show.node.label=$showBoostrapVals, cex=0.8, tip.color=tip.colors, main=\"$title\")
par(op)
dev.off()
~;

  run_R_script( $Rscript );
}


# execute an R-script
sub run_R_script
{
  my $Rscript = shift;

  my ($fh, $inFile) = tempfile("rTmpXXXX", SUFFIX => '.R', OPEN => 1, DIR => $tmpDir);
  print $fh "$Rscript";
  close $fh;

  my $outFile = $inFile . "out";
  my $cmd = "$R CMD BATCH --no-save --no-restore-data $inFile $outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

  open IN, "$outFile" or die "Cannot open $outFile for reading: $OS_ERROR";
  my $exitStatus = 1;

  foreach my $line (<IN>)
  {
    if ( $line =~ /Error/ )
    {
      print "R script crashed at\n$line";
      print "check $outFile for details\n";
      $exitStatus = 0;
      exit 1;
    }
  }
  close IN;
}

# parse a CSV partition table
sub read_part_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  my $headerStr = <IN>;
  foreach my $line (<IN>)
  {
    chomp $line;

    ##  $ clustername        : int  426 426 432 432 432 432 432 432 432 449 ...
    ##  $ bootstrap          : num  0.904 0.904 0.908 0.908 0.908 0.908 0.908 0.908 0.908 0.976 ...
    ##  $ leafname           : chr  "Lactobacillus_hordei" "Lactobacillus_mali_tRT_2" "Lactobacillus_nagelii" "Lactobacillus_vini" ...
    ##  $ branchPath         : num  0.0462 0.0525 0.0547 0.0546 0.0526 ...
    ##  $ medianOfDistances  : num  0.00651 0.00651 0.01502 0.01502 0.01502 ...
    ##  $ sequencesperCluster: int  2 2 7 7 7 7 7 7 7 2 ...

    my @f = split ",", $line;
    my $cltrId = shift @f;
    my $boot   = shift @f;
    my $leafId = shift @f;
    $tbl{ $leafId } = $cltrId;
  }
  close IN;

  return %tbl;
}

# write hash table to a file
sub writeTbl
{
  my ($rTbl, $outFile) = @_;
  my %tbl = %{$rTbl};
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

## Testing if two arrays are identical in a set-theoretic sense. That is that
## they have exactly the same set of elements.
sub setEqual
{
  my ($rA, $rB) = @_;

  my @a = @{$rA};
  my @b = @{$rB};
  my @c = comm(\@a, \@b);

  my $ret = 1;

  if (@c != @a || @c != @b)
  {
    warn "\n\n\tERROR: Elements of the two arrays do not match";
    print "\n\tNumber of elements in the first array: " . @a . "\n";
    print "\tNumber of elements in the second array: " . @b . "\n";
    print "\tNumber of common elements: " . @c . "\n";

    # writeArray(\@a, "a.txt");
    # writeArray(\@b, "b.txt");
    #print "\n\tNew taxon keys and fasta IDs written to a.txt and b.txt, respectively\n\n";

    if (@a > @b)
    {
      my @d = diff(\@a, \@b);
      print "\nElements a but not b:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    if (@b > @a)
    {
      my @d = diff(\@b, \@a);
      print "\nElements in b that are not a:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    $ret = 0;
  }

  return $ret;
}

sub create_mothur_script
{
    my (@arr) = @{$_[0]};

    my ($fh, $inFile) = tempfile("mothur.XXXX", SUFFIX => '', OPEN => 1, DIR => $tmpDir);
    foreach my $c (@arr)
    {
        print $fh $c . "\n";
    }
    print $fh "quit()\n";
    close $fh;

    return $inFile;
}

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

# parse sp_to_phGr.txt whose format is
#
#   Lactobacillus_iners	1227355	Firmicutes_group_6_V3V4
#   Lactobacillus_helveticus_acidophilus	920798	Firmicutes_group_6_V3V4
#   Gardnerella_vaginalis	402803	Actinobacteria_group_0_V3V4
#
# and generate a table
# genus => number of _sp species within that genus

sub getSpIdx
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    my ($sp, $n, $phGr ) = split /\s+/, $_;
    my ( $ge, $suffix, $idx ) = split "_", $sp;

    if ( $suffix && $suffix eq "sp" )
    {
      if ( $idx )
      {
	if ( ( exists $tbl{$ge} && $tbl{$ge} < $idx ) || ( !exists $tbl{$ge} ) )
	{
	  $tbl{$ge} = $idx;
	}
      }
      else
      {
	$tbl{$ge} = 1;
      }
    }
  }
  close IN;

  return %tbl;
}

exit 0;
