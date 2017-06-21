#!/usr/bin/env perl

=head1 NAME

  phylo_split_by_genus.pl

=head1 DESCRIPTION

  The script identifies clades of subRDP+CX+HB sequences using a CX+HB+(a few
  genera of subRDP) tree. The groups are called in the sequal phylo-groups.

  Input:
  - Tree
  - genus taxonomy of the leaves of the tree
  - genus taxonomy of subRDP sequences
  - phylo-group size threshold

  Using the above genus taxonomy tables, genus => size table (for full and
  restricted species size selections) will be generated. Post-order traversal of
  the genus condensed tree will be run at each time testing if the size of the
  corresponding subtree is less than a threshold.

  Output:
  - phylo_groups.tx: a two column table with seqID and phylo-group ID
  - og.tx:           a two column OG table: seqID and phylo-group ID only for OG seq's
  - phylo-group dirs with og.seqIDs file in each

  After this script is run

  - multi_select_seqs_bigMem.pl will be run with -s phylo_groups.tx to generate
    fasta file of FL subRDP+CX+HB sequences per phylo-group.
  - do the same for lineage files

  Once fasta and lineage files are generated, ginsi alignments are build and
  then FastTree is run on them. The resulting trees are rerooted using OG seq's.

  Next, taxonomy_cleaup.pl is run.

  The resulting species taxonomy is fed into V3V4 phylo-group dirs.

  Taxonomy cleanup is run on V3V4 phylo-groups.

  We use the data to build models - using phylo-part lineage generated in the
  taxonomy cleanup process.

  M&M sequences are classified using the above models.
  NOTE that we have non-redundant files that were generated in the previous run
  of pecan validation and we can run the classifier on those files to speedup
  the whole process.

  PECAN validation is run and final taxonomy assignments generated.

=head1 SYNOPSIS

  phylo_split_by_genus.pl -m <max phGr size> -t <tree file> -s <cltr size file> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--tree-file, -t>
  Genus condensed tree file generated with vicut_classify.pl

=item B<--tree-genus-tx-file, -s>
  Genus taxonomy of the leaves of the input tree of vicut_classify.pl

=item B<--subRDP-genus-tx-file, -g>
  subRDP seq's genus taxonomy file.

=item B<--max-phGr-size, -m>
  Maximum size of a phylo-group.

=item B<--genus-tx-chg-file, -c>
  genus_taxonomy_changes2.txt output file from vicut_classify.pl

=item B<--output-dir, -o>
  Output dir with the following files.

=item B<--report>
  Prints report.

=item B<--run-all>
  Ignore if ( ! -e ... ) statements.

=item B<--verbose, -v>
  Prints content of some output files. Default value: 5000.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/Banfield_contax

  phylo_split_by_genus.pl -m 5000 -t Banfield_contax_FL_subRDP_supplemnt_algn_goodHBids_100_1776tr_nr_rr.tree -s Banfield_contax_subRDP_supplt_FL_min1250_100_1776tr_genus_dir/genus_classification_results2.txt -g .. -c Banfield_contax_subRDP_supplt_FL_min1250_100_1776tr_genus_dir/genus_taxonomy_changes2.txt -o cx_hb_ssubRDP_FL_5k_phGr_dir

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd qw(abs_path);
use File::Temp qw/ tempfile /;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "tree-file|t=s"            => \my $treeFile,
  "tree-genus-tx-file|s=s"   => \my $treeGenusTxFile,
  "subRDP-genus-tx-file|g=s" => \my $subRDPgenusTxFile,
  "genus-tx-chg-file|c=s"    => \my $genusTxChgFile,
  "max-phGr-size|m=s"        => \my $maxCltrSize,
  "output-dir|o=s"           => \my $outDir,
  "report"                   => \my $report,
  "run-all"                  => \my $runAll,
  "quiet"                    => \my $quiet,
  "igs"                      => \my $igs,
  "verbose|v"                => \my $verbose,
  "debug"                    => \my $debug,
  "dry-run"                  => \my $dryRun,
  "help|h!"                  => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$outDir )
{
  print "ERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$treeFile )
{
  print "ERROR: Missing tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$treeGenusTxFile )
{
  print "ERROR: Missing genus taxonomy of the leaves of the tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$subRDPgenusTxFile )
{
  print "ERROR: Missing subRDP genus taxonomy file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$genusTxChgFile )
{
  print "ERROR: Missing genus taxonomy change file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $treeFile || ! -s $treeFile )
{
  warn "ERROR: Tree file $treeFile does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

if ( ! -e $treeGenusTxFile || ! -s $treeGenusTxFile )
{
  warn "ERROR: Genus size file does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

if ( ! -e $subRDPgenusTxFile || ! -s $subRDPgenusTxFile )
{
  warn "ERROR: Genus taxonomy file does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

if ( ! -e $genusTxChgFile || ! -s $genusTxChgFile )
{
  warn "ERROR: Genus taxonomy change file does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

# creating output directory
my $cmd = "mkdir -p $outDir";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

my $tmpDir = $outDir . "/temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}


my $quietStr = "";
if ( $quiet )
{
  $quietStr = "--quiet";
}

my $debugStr = "";
if ( $debug )
{
  $debugStr = "--debug";
  $quietStr = "";
}

my $verboseStr = "";
if ( $verbose )
{
  $verboseStr = "--verbose";
}

my $nw_labels             = "nw_labels";
my $nw_order              = "nw_order";
my $nw_condense           = "nw_condense";
my $nw_rename             = "nw_rename";
my $nw_prune              = "nw_prune";
my $nw_reroot             = "nw_reroot";
my $nw_clade              = "nw_clade";
my $uc2clstr2             = "uc2clstr2.pl";
my $extract_seq_IDs       = "extract_seq_IDs.pl";
my $select_seqs           = "select_seqs.pl";
my $rmGaps                = "rmGaps";
my $FastTree              = "FastTree";
my $R                     = "R";
my $fix_fasta_headers     = "fix_fasta_headers.pl";
my $mothur                = "mothur";
my $usearch6              = "usearch6.0.203_i86osx32";
my $vicut                 = "vicut";
my $readNewickFile        = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";
my $ginsi                 = "/usr/local/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)

if ( defined $igs )
{
  $nw_labels             = "/usr/local/projects/pgajer/bin/nw_labels";
  $nw_order              = "/usr/local/projects/pgajer/bin/nw_order";
  $nw_condense           = "/usr/local/projects/pgajer/bin/nw_condense";
  $nw_rename             = "/usr/local/projects/pgajer/bin/nw_rename";
  $nw_prune              = "/usr/local/projects/pgajer/bin/nw_prune";
  $nw_reroot             = "/usr/local/projects/pgajer/bin/nw_reroot";
  $nw_clade              = "/usr/local/projects/pgajer/bin/nw_clade";
  $uc2clstr2             = "/home/pgajer/devel/MCclassifier/perl/uc2clstr2.pl";
  $extract_seq_IDs       = "/home/pgajer/devel/MCclassifier/perl/extract_seq_IDs.pl";
  $select_seqs           = "/home/pgajer/devel/MCclassifier/perl/select_seqs.pl";
  $rmGaps                = "/usr/local/projects/pgajer/bin/rmGaps";
  $FastTree              = "/home/pgajer/bin/FastTree_no_openMP";
  $R                     = "/home/pgajer/bin/R";
  $fix_fasta_headers     = "/home/pgajer/devel/MCclassifier/perl/fix_fasta_headers.pl";
  $mothur                = "/usr/local/projects/pgajer/bin/mothur";
  $usearch6              = "/local/projects/pgajer/bin/usearch6.0.203_i86linux32";
  $vicut                 = "/usr/local/projects/pgajer/bin/vicut";
  $readNewickFile        = "/local/projects/pgajer/devel/MCclassifier/perl/read.newick.R";
  $ginsi                 = "/home/pgajer/bin/mafft --maxiterate 1000 --globalpair"; # MAFFT v7.310 (2017/Mar/17)
}

####################################################################
##                               MAIN
####################################################################


print "--- Parsing genus taxonomy change file\n";
my %txChgTbl = parse_tx_chg_file( $genusTxChgFile );

print "--- Generating from => to table for lost genera\n";
my %fromTo;
for my $from ( keys %txChgTbl )
{
  my @tos = sort { $txChgTbl{$from}{$b} <=> $txChgTbl{$from}{$a} } keys %{$txChgTbl{$from}}
  $fromTo{$from} = shift @tos;
}

print "--- Parsing genus taxonomy of the leaves of the tree\n";
my %trGeTbl = read_tbl( $treeGenusTxFile ); # trGeTbl{$ge} = ref to array of seq's of that genus

print "--- Parsing subRDP genus taxonomy file\n";
my %rdpGeTbl = read_tbl( $subRDPgenusTxFile ); # rdpGeTbl{$ge} = ref to array of seq's of that genus

my %geTbl = %trGeTbl;
for my $ge ( keys %geTbl )
{
  if ( exists $rdpGeTbl{$ge} )
  {
    push @{$geTbl{$ge}}, @{$rdpGeTbl{$ge}};
    delete $rdpGeTbl{$ge};
  }
}

for my $ge ( keys %rdpGeTbl )
{
  if ( exists $fromTo{$ge} && exists $trGeTbl{ $fromTo{$ge} } )
  {
    push @{$geTbl{ $fromTo{$ge} }), @{$rdpGeTbl{$ge}};
  }
}

my %geSize;
for my $ge ( keys %geTbl )
{
  $geSize{$ge} = @{$geTbl{$ge}};
}

print "--- Generating for all leaves of the condensed tree\n";
## since this is a tree generated vicut_classify.pl, there will be leaves with
## c<N> label names corresponding to vicut clusters with only query sequences or
## clusters of some genera present in more than one place (only largest cluster
## of this genus will carry its name and the rest will be turned into c<N>
## unless they merge with other genera in the cleanup process).

my @leaves = get_leaves( $treeFile );

# Expanding $geSize to c<N> leaves, setting their size to 0
for my $ge ( @leaves )
{
  if ( !exists $geSize{$ge} )
  {
    $geSize{$ge} = 0;
  }
}

my $geSizeFile = $outDir . "/genus_size_tbl.txt";
write_tbl( \%geSize, $geSizeFile );

print "--- Running cut_tree\n";
$cmd = "cut_tree -m $maxCltrSize -s $geSizeFile -t $treeFile -o $outDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

print "--- Parsing cut_tree results\n";
my $cutFile = "$outDir/thld_" . $maxCltrSize . "_tree_cut_tbl.txt";
my ($rphGr, $rphGrSize) = parse_cut_file( $cutFile );

my %phGrTbl  = %{$rphGr};     # phGr => ref to array of tree leaves of that phylo-group
my %phGrSize = %{$rphGrSize}; # phGr => size


print "--- Selecting an OG for each phylo-group\n";
my %ogTbl = get_ogs( $treeFile, \%phGrTbl ); # phGr => seqID


print "--- Generating phGr table\n";
my $phGrFile = "$outDir/thld_" . $maxCltrSize . "_phGr_tbl.txt";
open OUT, ">$phGrFile" or die "Cannot open $phGrFile for writing: $OS_ERROR\n";
for my $phGr ( keys %phGrTbl )
{
  print "\nProcessing $phGr\n\n";

  # Creating phylo-group directory
  my $phGrDir = $outDir . "/$phGr" . "_dir";
  $cmd = "mkdir -p $phGrDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  # Getting an OG sequence of $phGr
  if ( ! exists $ogTbl{$phGr} )
  {
    warn "\n\n\tERROR: $phGr not found in ogTbl";
    print "\n\n";
    exit 1;
  }

  my $ogSeq = $ogTbl{$phGr};
  print OUT "$ogSeq\t$phGr\n";

  # Creating an outgroup file in the phGr dir
  my $ogFile = $phGrDir . "/$phGr" . "_outgroup.seqIDs";
  my @ogs = ($ogSeq);
  write_array( \@ogs, $ogFile );

  # Finishing seqID => phGr table for the given phGr
  my @ges = @{ $phGrTbl{$phGr} };
  my @uqGes = unique( \@ges );
  for my $ge ( @uqGes )
  {
    if ( exists $geTbl{$ge} )
    {
      for ( @{$geTbl{$ge}} )
      {
        print OUT "$_\t$phGr\n";
      }
    }
  }
}
close OUT;

cleanup_tmp_files();

####################################################################
##                               SUBS
####################################################################

## parse_tx_chg_file

## parsing 3 col tbl of the form

# HB1021	Mannheimia	Basfia
# ConTax140371	Thermocrinis	Caldimicrobium
# ConTax140332	Thermocrinis	Caldimicrobium
# ConTax1229	Halopenitus	Halorubrum

# 2nd col = from, 3rd = to
# generating table txChgTbl{from}{to} = seq count
sub parse_tx_chg_file
{
  my $file = shift;

  my %txChgTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    chomp;
    my ($id, $from, $to) = split /\s+/;
    $txChgTbl{$from}{$to}++;
  }
  close IN;

  return %txChgTbl;
}


sub generate_lineage_file
{
  my ($phGr, $ogSeq, $liFile) = @_;

  my @cxhgIDs = @{ $cxhbTbl{$phGr} };

  my %liTbl;
  check_key( $ogSeq, \%cxhgLiTbl, "cxhgLiTbl" );
  $liTbl{$ogSeq} = $cxhgLiTbl{$ogSeq};

  my %geTbl; # genera of the given phylo-group
  for my $id ( @cxhgIDs )
  {
    check_key( $id, \%cxhgLiTbl, "cxhgLiTbl" );
    # if ( !exists $cxhgLiTbl{$id} )
    # {
    #   warn "\n\n\tERROR: $id not found in cxhgLiTbl";
    #   print "\n\n";
    #   exit 1;
    # }
    my $li = $cxhgLiTbl{$id};
    $liTbl{$id} = $li;
    my $ge = get_genus($li);
    $geTbl{$ge}++;
  }

  for my $ge ( keys %geTbl )
  {
    my @ids = @{ $ourDBgeTbl{$ge} };
    for my $id ( @ids )
    {
      $liTbl{$id} = $ourDBliTbl{$id};
    }
  }

  open OUT, ">$liFile" or die "Cannot open $liFile for writing: $OS_ERROR";
  for my $id ( keys %liTbl )
  {
    print OUT "$id\t" . $liTbl{$id} . "\n";
  }
  close OUT;
}

sub get_genus
{
  my $li = shift;

  # $li is assumed to have the format
  # "Root;Archaea;Euryarchaeota;Methanococci;Methanococcales;Methanococcaceae;Methanothermococcus;Methanothermococcus_okinawensis"
  my @f = split ";", $li;
  pop @f;
  my $ge = pop @f;

  return $ge;
}

sub check_key
{
  my ($key, $rtbl, $tblName) = @_;

  if ( !exists $rtbl->{$key} )
  {
    warn "\n\n\tERROR: $key does not exist in $tblName";
    print "\n\n";
    exit 1;
  }
}

sub parse_geTx_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in parse_genusTx_tbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    push @{ $tbl{$t} }, $id;
  }
  close IN;

  return %tbl;
}

sub parse_spTx_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in parse_genusTx_tbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  foreach (<IN>)
  {
    chomp;
    my ($id, $sp) = split /\s+/,$_;
    my @f = split "_", $sp;
    my $ge = shift @f;
    push @{ $tbl{$ge} }, $id;
  }
  close IN;

  return %tbl;
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

sub ge_li_tbl
{
  my $file = shift;

  my %geFaTbl;
  my %geOrTbl;
  my %geClTbl;
  my %gePhTbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  my $header = <IN>;
  for ( <IN> )
  {
    chomp;
    my @f = split /\s+/;

    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;

    $geFaTbl{$ge} = $fa;
    $geOrTbl{$ge} = $or;
    $geClTbl{$ge} = $cl;
    $gePhTbl{$ge} = $ph;
  }
  close IN;

  return (\%geFaTbl, \%geOrTbl, \%geClTbl, \%gePhTbl);
}


# Parsing vicut clustering results on genus tree with family annotation coming
# from the medoids data
sub parse_cltr_file
{
  my $file = shift;

  # file format

  # Id            	clstr	annot
  # Acanthopleuribacter	0	Acanthopleuribacteraceae
  # Bacteria_FibrobacteresAcidobacteria_group_Acidobacteria_unclassified_Acidobacteria_Thermoanaerobaculum_aquaticum_MP_01	1	NA
  # Bacteria_unclassified_bacteria_CG_uncultured_03	2	NA
  # Chroococcidiopsis	3	Xenococcaceae
  # Symploca	4	Phormidiaceae

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  my $header = <IN>;
  for ( <IN> )
  {
    chomp;
    my ($id, $cltr, $ann) = split /\s+/;
    $cltr = "c$cltr";
    push @{ $tbl{$cltr} }, $id;
  }
  close IN;

  return %tbl;
}

# Parsing the output of vicut_and_plot.pl script
sub parse_cut_file
{
  my $cutFile = shift;

  my %phGr;
  my %phGrSize;

  open IN, "$cutFile" or die "Cannot open $cutFile for reading: $OS_ERROR";
  for ( <IN> )
  {
    chomp;
    my ($id, $cltr, $size) = split /\s+/;
    $cltr = "phGr$cltr";
    push @{ $phGr{$cltr} }, $id;
    $phGrSize{$cltr} = $size;
  }
  close IN;

  return (\%phGr, \%phGrSize);
}

# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}


# extract sequence count of a fasta file
sub seq_count
{
  my $file = shift;

  my $wcline = qx/ grep -c '>' $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
}

sub build_tree
{
  my ($algnFile, $treeFile) = @_;

  $cmd = "$FastTree -nt $algnFile > $treeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub reroot_tree
{
  my ($treeFile, $rog, $rrTreeFile ) = @_;

  my @ogs = @{$rog};

  $cmd = "$nw_reroot $treeFile @ogs | $nw_order -  > $rrTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub build_ann_seqID_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile2 | $nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}


sub build_ann_cltr_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  # $annFile file format

  # readId	clstr	annot
  # Chitinivibrio	0	Fibrobacteres
  # Ignavibacterium	1	Chlorobi
  # Caldisericum	2	Caldiserica

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$3\"__\"\$2}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile2 | $nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}



sub build_ann_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile | $nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub plot_color_tree
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

## plot tree with clade colors
sub plot_tree
{
  my ($treeFile, $title, $pdfFile) = @_;

  my $showBoostrapVals = "T";

  if (!defined $title)
  {
    $title = "";
  }

  my $Rscript = qq~

source(\"$readNewickFile\")
require(phytools)

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

# print elements of a hash table
sub print_tbl
{
  my ($rTbl, $r) = @_;

  map {print "$_\t" . $rTbl->{$_} . "\n"} @{$r};
}

# this is a version of print_formated_tbl() where sorting w/r values are
# performed within this routine
sub print_formated_freq_tbl
{
  my $rTbl = shift;

  my @args = sort { $rTbl->{$b} <=> $rTbl->{$a} } keys %{$rTbl};
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
    print "WARNING: tbl value not defined for $_\n" if !exists $rTbl->{$_};
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  print "\n";

  return @args;
}

# print elements of a hash table so that arguments are aligned
sub print_formated_tbl
{
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
    print "WARNING: tbl value not defined for $_\n" if !exists $rTbl->{$_};
    print "\t$_$pad" . $rTbl->{$_} . "\n";
  }
  #print "\n";
}

sub mothur_align_and_add
{
  my ($candidateFile, $templateFile, $nProc) = @_;

  my $seqCountBefore = seq_count( $templateFile );

  my @tmp;
  push (@tmp,"align.seqs(candidate=$candidateFile, template=$templateFile, flip=T, processors=$nProc)");

  print_array( \@tmp, "mothur commands" ) if ( $debug || $verbose );

  my $scriptFile = create_mothur_script( \@tmp );
  my $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;

  my $mothurAlgnFile = $candidateFile;
  $mothurAlgnFile =~ s/fa$/align/;
  if ( ! -e $mothurAlgnFile )
  {
    warn "\n\n\tERROR: Count not find $mothurAlgnFile";
    print "\n\n";
    exit 1;
  }

  $cmd = "cat $mothurAlgnFile >> $templateFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  my $seqCountAfter = seq_count( $templateFile );

  if ( !$quiet )
  {
    print "\nAlignment file line count BEFORE: $seqCountBefore\n";
    print   "Alignment file line count AFTER:  $seqCountAfter\n\n";
  }

  return ($seqCountBefore, $seqCountAfter);
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

# write array to a file (one column format)
sub write_array
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

# write hash table to a file
sub write_tbl
{
  my ($rTbl, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} sort keys %tbl;
  close OUT;
}

# write hash table to a file
sub write_sorted_tbl
{
  my ($rTbl, $r, $outFile) = @_;

  my %tbl = %{$rTbl};

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} @{$r};
  close OUT;
}

sub run_vicut_no_query
{
  my ($treeFile, $annFile, $vicutDir) = @_;

  my $cmd = "$vicut $quietStr -t $treeFile -a $annFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub cleanup_tmp_files
{
  my $cmd = "rm -f $tmpDir/*";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

## Create a condensed tree given a tree, translation table and a file name of the
## condensed tree
sub condense_tree
{
  my ($treeFile, $txFile, $condTreeFile) = @_;

  my %txTbl = read_tbl($txFile);
  my @txs = keys %txTbl;

  ## extracting leave IDs
  my $treeLeavesFile = "tmp_tree.leaves";
  $cmd = "rm -f $treeLeavesFile; $nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @treeLeaves = read_array ( $treeLeavesFile );
  unlink ( $treeLeavesFile );

  # if ( !setequal( \@txs, \@treeLeaves ) )
  # {
  #   warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
  #   print "\n\n";
  #   exit 1;
  # }

  my $sppTreeFile = "tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; $nw_rename $treeFile $txFile > $sppTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  $cmd = "rm -f $condTreeFile; $nw_condense $sppTreeFile > $condTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ( $debug )
  {
    print "\n\ntree:        $treeFile\n";
    print     "tx:          $txFile\n";
    print     "sppTree:     $sppTreeFile\n";
    print     "condTree: $condTreeFile\n\n";
  }
}

# read table with one column
sub read_array
{
  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -e $file )
  {
    warn "\n\n\nERROR: $file does not exist";
    print "\n\n";
    exit 1;
  }

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my ($file, $skipHeader) = @_;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  if ( $skipHeader )
  {
    my $header = <IN>;
  }
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub get_leaves
{
  my $treeFile = shift;

  my ($fh, $leavesFile) = tempfile("leaves.XXXX", SUFFIX => '', OPEN => 0, DIR => $tmpDir);
  my $cmd = "$nw_labels -I $treeFile > $leavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  my @a = read_array($leavesFile);

  return @a;
}

sub get_ogs
{
  my ( $treeFile, $rphGrTbl ) = @_;

  my %phGrTbl = %{$rphGrTbl};
  $" = ",";

  my $Rscript = qq~

tr <- read.tree(\"$treeFile\")
for (
ids <- c(
p <- first.common.ancestor( tr, ids )
p <- tree.parent.node( tr, p ) # get the grandparent of p
ids2.idx <- node.leaves( tr, p )
ids2 <- tr\$tip.label[ids2.idx]
cmd <- sprintf(\"nw_clade %s \", treeFile)
for ( id in ids2 )
{
    cmd <- paste(cmd, id)
}
cmd <- paste(cmd, \" > \", \"$cladeTreeFile\")
system(cmd)

cmd <- sprintf("nw_rename %s %s > %s", cladeTreeFile, \"$annFile\", pFile2)
system(cmd)

    rids  <- as.character(chb.tx[ids])
    rids2 <- as.character(chb.tx[ids2])

    clade.tr <- tree.read( cladeTreeFile2 )
    nLeaves <- length(clade.tr$tip.label)

    p <- first.common.ancestor( clade.tr, rids )
    cl <- strsplit(title, " ")[[1]][1]

    tip.colors <- rep( 1, length(clade.tr$tip.label) )
    names(tip.colors) <- clade.tr$tip.label
    tip.colors[rids] <- 4 # blue

    i <- grep("ConTax", rids)
    cid <- rids[i]
    tip.colors[cid] <- 2

    figH <- 8
    figW <- 12
    if ( nLeaves >= 50 )
    {
        figH <- 6.0/50.0 * ( nLeaves - 50) + 10
        figW <- 6.0/50.0 * ( nLeaves - 50) + 6
    }

    clade.tr$node.label <- rep("", length(clade.tr$node.label))

    pdf(pdfFile, width=figW, height=figH)
    op <- par(mar=c(0,0,1.5,0), mgp=c(2.85,0.6,0),tcl = -0.3)
    plot(clade.tr,type="phylogram", use.edge.length=FALSE, no.margin=FALSE, show.node.label=T, cex=0.8, tip.color=tip.colors, main=title)
    nodelabels(cl, p)
    par(op)
    dev.off()

        ~;

    unlink $idsFile;
}

exit 0;
