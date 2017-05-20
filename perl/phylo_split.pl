#!/usr/bin/env perl

=head1 NAME

  phylo_split.pl

=head1 DESCRIPTION

  The aim of the script is to generate a split of all bacterial (and Archeal)
  sequences into phylogenetically sound sub-groups.

  1. Build a lineage file of all sequences of a db.

  2. Build a table assigning to each genus a ref to an array of seqIDs of
  sequences of that family.

  3. For each family pick a random sequence from the most abundant species of
  that family.

  4. Get fasta file of these sequences.
  5. Align them.
  6. Build a tree.

  7. Add family size to the label of each leaf (family name) of the tree.

  8. Devise an algorithm (a version of vicut tree cutting) that cuts that tree at
  different nodes so that the total size of the cluster is less than the maxSize.

=head1 SYNOPSIS

  phylo_split.pl -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--output-dir, -o>
  Output dir with the following files.

  family.fa
  family_algn.fa
  family.tree
  family.size

=item B<--num-proc, -p>
  Number of processors to be used. Default value: 8.

=item B<--report>
  Prints report.

=item B<--max-cltr-size, m>
  Maximum size of a cluster of families.

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

  phylo_split.pl -o phylo_split_dir
  phylo_split.pl --igs -o phylo_split_dir

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

my $maxCltrSize = 3000;
my $nProc       = 8;

GetOptions(
  "output-dir|o=s"      => \my $outDir,
  "max-cltr-size|m=s"   => \$maxCltrSize,
  "num-proc|p=i"        => \$nProc,
  "report"              => \my $report,
  "run-all"             => \my $runAll,
  "quiet"               => \my $quiet,
  "igs"                 => \my $igs,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
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

my $nProcStr = "";
if ( $nProc )
{
  $nProcStr = "--thread $nProc";
}

my $spGeFile              = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/species_genus_tbl_may19_2017.txt";
my $geLiFile              = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/microcontax_genus_lineage_tbl.txt"; # NOTE: this one has a header !

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
my $readNewickFile        = "/Users/pgajer/.Rlocal/read.newick.R";
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
  $ginsi                 = "/home/pgajer/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)

  $spGeFile              = "/usr/local/projects/pgajer/devel/MCextras/data/species_genus_tbl_may19_2017.txt";
  $geLiFile              = "/usr/local/projects/pgajer/devel/MCextras/data/microcontax_genus_lineage_tbl.txt"; # NOTE: this one has a header !

  $nProcStr              = "";
  $nProc                 = 1;
}

####################################################################
##                               MAIN
####################################################################

print "--- Identifying lineage files\n";
my @liFiles;
if ( $igs )
{
  @liFiles = get_li_files_igs();
}
else
{
  @liFiles = get_li_files();
}

if ( $debug )
{
  print_array( \@liFiles, "\nliFiles" );
  print "\n\n";
  exit;
}

print "--- Parsing microcontax pkg genus-lineage tbl\n";
my ($rgeFaTbl, $rgeOrTbl, $rgeClTbl, $rgePhTbl) = ge_li_tbl( $geLiFile );

my %geFaTbl = %{$rgeFaTbl};
my %geOrTbl = %{$rgeOrTbl};
my %geClTbl = %{$rgeClTbl};
my %gePhTbl = %{$rgePhTbl};

$geFaTbl{"OG"} = "OG";
$geOrTbl{"OG"} = "OG";
$geClTbl{"OG"} = "OG";
$gePhTbl{"OG"} = "OG";


print "--- Parsing species-genus tbl\n";
my %spGeTbl = read_tbl( $spGeFile );


print "--- Extracting taxonomy of V3V4 ref seq's\n";

#my %txTbl;
my %spSeqIDs;
my %geFileFreq;
my %geSpFreq;
my %liFile2phGr;

foreach my $file ( @liFiles )
{
  print "\rProcessing $file        ";
  my @a = split "/", $file;
  my $phGr = pop @a;
  $phGr =~ s/\.lineage$//;
  $liFile2phGr{$file} = $phGr;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for my $lineage (<IN>)
  {
    chomp $lineage;
    my ($id, $li) = split /\s+/, $lineage;
    my @f = split ";", $li;

    my $sp = pop @f;
    #$txTbl{$id} = $sp;

    if ( ! exists $spGeTbl{$sp} )
    {
      my ($g, $suffix) = split "_", @f;
      print "\nWARNING: $sp not found in spGeTbl - using $g\n";
      $spGeTbl{$sp} = $g;
      #exit 1;
    }
    my $ge = $spGeTbl{$sp};
    $geFileFreq{$ge}{$file}++;
    $geSpFreq{$phGr}{$ge}{$sp}++;
    push @{ $spSeqIDs{$phGr}{$sp} }, $id;
  }
  close IN;
}

print "\r                                                                          \n";


print "--- Checking if all out db genera are found in the master genus-lineage table\n";
my @ges = values %spGeTbl;
@ges = unique( \@ges );
for my $ge ( @ges  )
{
  if ( ! exists $geFaTbl{$ge} )
  {
    warn "WARNING: $ge not found in geFaTbl\n";
    #print "$ge\n";
  }
}

if ( 0 )
{
  my @ourGenera = keys %geFileFreq;
  my @masterGenera = keys %geFaTbl;
  my @d = diff( \@masterGenera, \@ourGenera );
  my $geNotInOurDB = $outDir . "/genera_not_in_our_db.txt";
  write_array( \@d, $geNotInOurDB );
  print "\n\nGenera not in our db written to $geNotInOurDB\n\n";
}

# 3. For each genus pick a random sequence from the most abundant species of
# that genus.


my %geRefSeqID;
my %geTx;

my $geTxFile = $outDir . "/genus.tx";
my $algnFile = $outDir . "/genus_algn.fa";
my $ogSeqID  = "S000414080";

my @ogs   = ($ogSeqID);
my %ogInd = map { $_ => 1 } @ogs;

my $geFile = $outDir . "/genus.fa";
if ( ! -e $geFile || ! -s $geFile || $runAll )
{
  $runAll = 1;
  unlink( $geFile ) if $runAll;

  print "--- Picking a represetative sequence of each family\n";
  print "--- Building fasta file of genus represetative sequences\n";

  for my $ge ( keys %geFileFreq )
  {
    print "\rProcessing $ge                                      ";
    my %fileFreq = %{ $geFileFreq{$ge} };
    my @files = sort { $fileFreq{$b} <=> $fileFreq{$a} } keys %fileFreq;
    if ( @files > 1 )
    {
      print "\nWARNING: $ge is present in more than one phylogroups\n";
      #print "files: @files\n";
      for my $f ( @files )
      {
	my $phGr = $liFile2phGr{$f};
	print "$phGr\t" . $fileFreq{$f} . "\n";
      }
      print "\n";
    }
    my $file = $files[0];
    my $phGr = $liFile2phGr{$file};
    #print "\n\nli file: $file\n";
    $file =~ s/lineage$/ge/;
    #print "ge file: $file\n"; exit;

    my %spFreq = %{ $geSpFreq{$phGr}{$ge} };
    my @spp = sort{ $spFreq{$b} <=> $spFreq{$a} } keys %spFreq;
    my $sp = $spp[0]; # species with most representative seq's
    my @ids = @{ $spSeqIDs{$phGr}{$sp} };
    my $refSeqID = $ids[rand @ids];
    $geRefSeqID{$ge} = $refSeqID;
    $geTx{$refSeqID} = $ge;

    if ( ! -e $file )
    {
      my $ginsiFile = $file;
      $ginsiFile =~ s/\.ge$/_ginsi_algn.fa/;
      print "\n\nginsiFile not found: $ginsiFile\n" if ! -e $ginsiFile;
      $cmd = "$rmGaps -i $ginsiFile -o $file";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) geiled with exit code: $?" if !$dryRun;
    }

    # select ref sequence from $file and append it to the target genus ref seq's
    # gesta file
    ##my $refSeqID = add_rand_seq_to_gesta( $file, $geRefSeqID{$ge}, $geFile );
    add_seq_to_fasta( $file, $geRefSeqID{$ge}, $geFile );
  }

  print "--- Writing genus taxonomy to a file\n";
  $geTx{$ogSeqID} = "OG";
  write_tbl( \%geTx, $geTxFile );
}
else
{
  %geTx = read_tbl( $geTxFile );

  for my $refSeqID ( keys %geTx )
  {
    my $ge = $geTx{$refSeqID};
    $geRefSeqID{$ge} = $refSeqID;
  }
}

my $nGenera = keys %geFileFreq;
my $geFileSize = seq_count( $geFile );
if ( $geFileSize != $nGenera  )
{
  warn "\n\n\tERROR: Number of gemilies: $nGenera";
  print "but\n";
  print "$geFile has $geFileSize seq's\n\n";
  exit 1;
}

if ( ! -e $algnFile || ! -s $algnFile || $runAll )
{
  print "\r--- Generating ginsi alignment                              \n";

  $runAll = 1;
  unlink( $algnFile ) if $runAll;
  ginsi_algn( $geFile, $algnFile );

  ## adding outgroup seq to the alignment
  print "--- Aligning OG sequence to the alignment\n";
  my $ogGeFile = "/Users/pgajer/devel/MCextras/data/RDP/old/Archaea_Caldococcus_noboribetus_S000414080.fa";

  my ($seqCountBefore, $seqCountAfter) = mothur_align_and_add( $ogGeFile, $algnFile, $nProc );
}



print "--- Building a tree\n";
my $treeFile          = $outDir . "/genus_ref_seqs.tree";
my $geSeqIDsTreeFile  = $outDir . "/genus__seqID.tree";
my $geTreeFile        = $outDir . "/genus.tree";
if ( ! -e $treeFile || ! -s $treeFile || $runAll )
{
  my $nrTreeFile = $outDir . "/genus_not_rooted.tree";
  if ( $runAll )
  {
    unlink( $nrTreeFile );
    unlink( $treeFile );
    unlink( $geSeqIDsTreeFile );
  }

  build_tree( $algnFile, $nrTreeFile );

  print "--- Rerooting the tree using Archeal outgroup sequence\n";
  reroot_tree( $nrTreeFile, \@ogs, $treeFile );

  # print "--- Generating tree with <genus name>_<seqID> labels at leaves\n";
  # build_ann_seqID_tree( $treeFile, $geTxFile, $geSeqIDsTreeFile );

  print "--- Generating tree with <genus name> labels at leaves\n";
  build_ann_tree( $treeFile, $geTxFile, $geTreeFile );
}

my @leaves = get_leaves( $geTreeFile );

##
## Phylum analysis
##
print "--- Phylum analysis\n";

print "--- Genus => Phylum table in the order of the genus tree leaves\n";
my $gePhFile = $outDir . "/genus_vs_phylum.txt";
write_sorted_tbl( \%gePhTbl, \@leaves, $gePhFile );

print "--- Generating tree with <phylum name> labels at leaves\n";
my $phTreeFile = $outDir . "/phylum.tree";
build_ann_tree( $geTreeFile, $gePhFile, $phTreeFile );

print "--- Generating condensed phylum tree\n";
my $condPhTreeFile = $outDir . "/condensed_phylum.tree";
condense_tree_only( $phTreeFile, $condPhTreeFile );

print "--- Generating pdf of the condensed phylum tree\n";
$condPhTreeFile = abs_path( $condPhTreeFile );
my $pdfTreeFile = $outDir . "/condensed_phylum_tree.pdf";
plot_tree( $condPhTreeFile, "Condensed Phylum", $pdfTreeFile );

if ( $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

print "--- Running vicut using phylum annotation\n";
my $phVicutDir = $outDir . "/phylum_vicut_dir";
run_vicut_no_query( $geTreeFile, $gePhFile, $phVicutDir );


##
## Class analysis
##
print "--- Class analysis\n";

print "--- Genus => Class table in the order of the genus tree leaves\n";
print "\nGenus => Class on genus tree leaves\n";
my $geClFile = $outDir . "/genus_vs_class.txt";
write_sorted_tbl( \%geClTbl, \@leaves, $geClFile );

print "--- Generating tree with <class name> labels at leaves\n";
my $clTreeFile = $outDir . "/class.tree";
build_ann_tree( $geTreeFile, $geClFile, $clTreeFile );

print "--- Generating condensed class tree\n";
my $condClTreeFile = $outDir . "/condensed_class.tree";
condense_tree_only( $clTreeFile, $condClTreeFile );

print "--- Generating pdf of the condensed class tree\n";
$condClTreeFile = abs_path( $condClTreeFile );
$pdfTreeFile = $outDir . "/condensed_class_tree.pdf";
plot_tree( $condClTreeFile, "Condensed Class", $pdfTreeFile );

if ( $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

print "--- Running vicut using class annotation\n";
my $clVicutDir = $outDir . "/class_vicut_dir";
run_vicut_no_query( $geTreeFile, $geClFile, $clVicutDir );

##
## Order analysis
##
print "--- Order analysis\n";

print "--- Genus => Order table in the order of the genus tree leaves\n";
print "\nGenus => Order on genus tree leaves\n";
#print_formated_tbl( \%geOrTbl, \@leaves );
my $geOrFile = $outDir . "/genus_vs_order.txt";
write_sorted_tbl( \%geOrTbl, \@leaves, $geOrFile );

print "--- Generating tree with <order name> labels at leaves\n";
my $orTreeFile = $outDir . "/order.tree";
build_ann_tree( $geTreeFile, $geOrFile, $orTreeFile );

print "--- Generating condensed order tree\n";
my $condOrTreeFile = $outDir . "/condensed_order.tree";
condense_tree_only( $orTreeFile, $condOrTreeFile );

print "--- Generating pdf of the condensed order tree\n";
$condOrTreeFile = abs_path( $condOrTreeFile );
$pdfTreeFile = $outDir . "/condensed_order_tree.pdf";
plot_tree( $condOrTreeFile, "Condensed Order", $pdfTreeFile );

if ( $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

print "--- Running vicut using order annotation\n";
my $orVicutDir = $outDir . "/order_vicut_dir";
run_vicut_no_query( $geTreeFile, $geOrFile, $orVicutDir );


##
## Family analysis
##
print "--- Family analysis\n";

print "--- Genus => Family table in the order of the genus tree leaves\n";
my $geFaFile = $outDir . "/genus_vs_family.txt";
write_sorted_tbl( \%geFaTbl, \@leaves, $geFaFile );

print "--- Generating tree with <family name> labels at leaves\n";
my $faTreeFile = $outDir . "/family.tree";
build_ann_tree( $geTreeFile, $geFaFile, $faTreeFile );

print "--- Generating condensed family tree\n";
my $condFaTreeFile = $outDir . "/condensed_family.tree";
condense_tree_only( $faTreeFile, $condFaTreeFile );

print "--- Generating pdf of the condensed family tree\n";
$condFaTreeFile = abs_path( $condFaTreeFile );
$pdfTreeFile = $outDir . "/condensed_family_tree.pdf";
plot_tree( $condFaTreeFile, "Condensed Family", $pdfTreeFile );

if ( $OSNAME eq "darwin")
{
  $cmd = "open $pdfTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

print "--- Running vicut using family annotation\n";
my $faVicutDir = $outDir . "/family_vicut_dir";
run_vicut_no_query( $geTreeFile, $geFaFile, $faVicutDir );


print "\n\n\tSuccessfully finished all tasks\n";
print "\n\n";


if (0)
{
  print "--- Testing if OG seq's form a monophylectic clade at the top or bottom of the tree\n";
  if ( test_OG( $geSeqIDsTreeFile, \%ogInd ) != 0 )
  {
    warn "\n\n\tERROR: There is an issue with the tree outgroup seq's";
    print "\n\tPlease check out the following trees\n";
    print "\t$treeFile\n";
    print "\t$geSeqIDsTreeFile\n";
    print "\n\n";

    exit 1;
  }
}

cleanup_tmp_files();

####################################################################
##                               SUBS
####################################################################


# print array to stdout
sub print_array
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

sub get_li_files_igs
{
  my $baseDir = "/usr/local/projects/pgajer/devel/MCextras/data/RDP/V3V4/";

  my @liFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4.lineage",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4.lineage",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4.lineage",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4.lineage",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4.lineage",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4.lineage",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4.lineage",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4.lineage",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4.lineage",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4.lineage");

  my @liFiles = map{ $baseDir . $_ } @liFiles0;

  return @liFiles;
}

sub get_li_files
{
  my $baseDirExt = "/Users/pgajer/devel/MCextras/data/vag_exp_V3V4_phGrps_May16_dir/";

  my @liFilesExt0 = ("Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4.lineage",
		     "Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4.lineage",
		     "Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4.lineage",
		     "Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4.lineage",
		     "Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4.lineage",
		     "Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4.lineage",
		     "Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4.lineage",
		     "Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4.lineage",
		     "Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4.lineage",
		     "Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4.lineage",
		     "Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4.lineage",
		     "Fusobacteria_V3V4_dir/Fusobacteria_V3V4.lineage",
		     "phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4.lineage",
		     "Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4.lineage",
		     "Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4.lineage",
		     "Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4.lineage",
		     "Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4.lineage",
		     "Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4.lineage",
		     "Tenericutes_V3V4_dir/Tenericutes_V3V4.lineage");

  my @liFilesExt1 = map{ $_ = $baseDirExt . $_ } @liFilesExt0;


  my @extPhGrs = ("Actinobacteria_group_0_V3V4",
		  "Actinobacteria_group_1_V3V4",
		  "Actinobacteria_group_2_V3V4",
		  "Bacteroidetes_group_2_V3V4",
		  "Firmicutes_group_0_V3V4",
		  "Firmicutes_group_1_V3V4",
		  "Firmicutes_group_2_V3V4",
		  "Firmicutes_group_3_V3V4",
		  "Firmicutes_group_4_V3V4",
		  "Firmicutes_group_5_V3V4",
		  "Firmicutes_group_6_V3V4",
		  "Fusobacteria_V3V4",
		  "phyla_lessthen_1k_wOG_V3V4",
		  "Proteobacteria_group_10_V3V4",
		  "Proteobacteria_group_15_V3V4",
		  "Proteobacteria_group_17_V3V4",
		  "Proteobacteria_group_3_V3V4",
		  "Proteobacteria_group_9_V3V4",
		  "Tenericutes_V3V4");

  my %liExtTbl;
  for my $i ( 0..$#extPhGrs )
  {
    $liExtTbl{$extPhGrs[$i]} = $liFilesExt1[$i];
  }


  my $baseDir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";
  my @liFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4.lineage",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4.lineage",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4.lineage",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4.lineage",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4.lineage",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4.lineage",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4.lineage",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4.lineage",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4.lineage",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4.lineage");

  my @liFiles1 = map{ $_ = $baseDir . $_ } @liFiles0;
  ## print "lineageFile: @liFiles\n";

  my @allPhGrs = ("Actinobacteria_group_0_V3V4",
		  "Actinobacteria_group_1_V3V4",
		  "Actinobacteria_group_2_V3V4",
		  "Actinobacteria_group_3_V3V4",
		  "Actinobacteria_group_4_V3V4",
		  "Actinobacteria_group_5_V3V4",
		  "Bacteroidetes_group_0_V3V4",
		  "Bacteroidetes_group_1_V3V4",
		  "Bacteroidetes_group_2_V3V4",
		  "Bacteroidetes_group_3_V3V4",
		  "Chloroflexi_V3V4",
		  "Deinococcus_Thermus_V3V4",
		  "Fusobacteria_V3V4",
		  "Nitrospirae_V3V4",
		  "Planctomycetes_V3V4",
		  "Spirochaetes_V3V4",
		  "Tenericutes_V3V4",
		  "Verrucomicrobia_V3V4",
		  "phyla_lessthen_1k_wOG_V3V4",
		  "Firmicutes_group_0_V3V4",
		  "Firmicutes_group_1_V3V4",
		  "Firmicutes_group_2_V3V4",
		  "Firmicutes_group_3_V3V4",
		  "Firmicutes_group_4_V3V4",
		  "Firmicutes_group_5_V3V4",
		  "Firmicutes_group_6_V3V4",
		  "Proteobacteria_group_0_V3V4",
		  "Proteobacteria_group_10_V3V4",
		  "Proteobacteria_group_11_V3V4",
		  "Proteobacteria_group_12_V3V4",
		  "Proteobacteria_group_13_V3V4",
		  "Proteobacteria_group_14_V3V4",
		  "Proteobacteria_group_15_V3V4",
		  "Proteobacteria_group_17_V3V4",
		  "Proteobacteria_group_1_V3V4",
		  "Proteobacteria_group_2_V3V4",
		  "Proteobacteria_group_3_V3V4",
		  "Proteobacteria_group_4_V3V4",
		  "Proteobacteria_group_5_V3V4",
		  "Proteobacteria_group_6_V3V4",
		  "Proteobacteria_group_7_V3V4",
		  "Proteobacteria_group_8_V3V4",
		  "Proteobacteria_group_9_V3V4");

  my %liTbl;
  for my $i ( 0..$#allPhGrs )
  {
    $liTbl{$allPhGrs[$i]} = $liFiles1[$i];
  }

  my @liFiles;
  for ( @allPhGrs )
  {
    if ( exists $liExtTbl{$_} )
    {
      push @liFiles, $liExtTbl{$_};
    }
    else
    {
      push @liFiles, $liTbl{$_};
    }
  }

  return @liFiles;
}

# Select $seqID from $file and append it to $faFile
sub add_seq_to_fasta
{
  my ( $file, $seqID, $faFile ) = @_;

  # select sequence with $seqID from $file
  my ($seqIdFH, $seqIdFile) = tempfile( "tmp.XXXX", SUFFIX => 'seqID', OPEN => 1, DIR => $tmpDir );
  print $seqIdFH $seqID;
  close $seqIdFH;

  my ($seqFaFH, $seqFaFile) = tempfile( "tmp.XXXX", SUFFIX => 'fa', OPEN => 0, DIR => $tmpDir );
  my $cmd = "$select_seqs $quietStr -s $seqIdFile -i $file -o $seqFaFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  if ( seq_count( $seqFaFile ) != 1 )
  {
    warn "\n\n\tERROR: $seqID not found in $file";
    print "\n\n";
    exit 1;
  }

  # appending $seqFaFile to $faFile
  $cmd = "cat $seqFaFile >> $faFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  unlink( $seqIdFile );
  unlink( $seqFaFile );
}

# extract line count of a file
sub line_count
{
  my $file = shift;

  my $wcline = qx/ wc -l $file /;
  $wcline =~ s/^\s+//;
  my ($lcount, $str) = split /\s+/, $wcline;

  return $lcount;
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

sub ginsi_algn
{
  my ($faFile, $algnFile) = @_;

  my $cmd = "rm -f $algnFile; $ginsi --inputorder $quietStr $nProcStr $faFile > $algnFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
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

sub build_ann_tree
{
  my ($treeFile, $annFile, $annTreeFile ) = @_;

  $cmd = "rm -f $annTreeFile; $nw_rename $treeFile $annFile | $nw_order -  > $annTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
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
sub print_formated_tbl{

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

# test if OG seq's form one or more clusters in the tree
sub test_OG
{
  my ($treeFile, $rogInd) = @_;

  my %ogInd = %{$rogInd};

  my $debug_test_OG = 0;

  my $ret = 0;

  print "\t--- Extracting leaves from $treeFile\n" if $debug_test_OG;
  my $treeLeavesFile = "$outDir" . "/genus_sppSeqIDs.leaves";
  my $cmd = "rm -f $treeLeavesFile; $nw_labels -I $treeFile > $treeLeavesFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "\t--- Reading leaves\n" if $debug;
  my @leaves = read_array($treeLeavesFile);

  print "\t--- Checking the number of clusters formed by OG seqs\n" if $debug;
  my @ogIdx;
  for my $i (0..$#leaves)
  {
    if ( exists $ogInd{$leaves[$i]} )
    {
      push @ogIdx, $i;
    }
  }

  print_array(\@ogIdx, "\nPositions of OG seq's") if ($debug);

  ## identifying consecutive indices ranges
  my @start;
  my @end;

  push @start, $ogIdx[0];
  if (@ogIdx>1)
  {
    for my $i (1..$#ogIdx)
    {
      ##if ( !$foundEnd && $ogIdx[$i] != $start[$rIdx] )
      if ( $ogIdx[$i-1]+1 != $ogIdx[$i] )
      {
	#$foundEnd = 1;
	push @end, $ogIdx[$i-1];
	push @start, $ogIdx[$i];
      }
      if ($i==$#ogIdx)
      {
	push @end, $ogIdx[$i];
      }

      if (0 && $debug)
      {
	print "\ni: $i\n";
	print_array(\@start, "start");
	print_array(\@end, "end");
      }
    }
  }
  else
  {
    push @end, $ogIdx[0];
  }

  my @ogPos1;  # OG positions
  for my $i (0..$#start)
  {
    push @ogPos1, ($start[$i] .. $end[$i]);
  }

  my @og = @leaves[@ogPos1];
  my @ogBig = @leaves[($start[0] .. $end[$#end])];

  print_array( \@og, "\nOutgroup elements" ) if $debug;

  if ( scalar(@start) != scalar(@end) )
  {
    warn "\n\n\tERROR: start and end arrays have different lengths!";
    print "length(start): " . @start . "\n";
    print "length(end): " . @end . "\n\n";
    $ret = 1;
  }

  my @rangeSize;
  for my $i (0..$#start)
  {
    push @rangeSize, ($end[$i] - $start[$i]+1);
  }

  if ($debug)
  {
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";
  }

  if (@rangeSize>1)
  {
    warn "\n\n\tERROR: Detected multiple OG clusters";
    print "\n\n";

    my $imax = argmax( \@rangeSize );
    print "imax: $imax\n";
    print "Maximal range size: " . $rangeSize[$imax] . "\n";

    my $minCladeSize = @leaves;
    my $minCladeSizeIdx = $imax;
    print "Clade size of each cluster of maximal range size\n";
    print "\nidx\tstart\tend\trgSize\tcladeSize\n";
    for my $i (0..$#rangeSize)
    {
      #if ($rangeSize[$i] == $rangeSize[$imax])
      if (1)
      {
	my @pos = ($start[$i] .. $end[$i]);
	my @og = @leaves[@pos];

	#print "\t--- Extracting the clade of OG sequences\n";
	my $ogCladeTreeFile = "$outDir/genus" . "_clade.tree";
	$cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$outDir/genus" . "_clade.leaves";
	$cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Reading the leaves\n" if $debug;
	my @ogCladeLeaves = read_array($ogCladeTreeLeavesFile);

	print "i: $i\t$start[$i]\t$end[$i]\t$rangeSize[$i]\t" . @ogCladeLeaves . "\n";
	if ( @ogCladeLeaves < $minCladeSize )
	{
	  $minCladeSize = @ogCladeLeaves;
	  $minCladeSizeIdx = $i;
	}
      }
    }

    # $imax = $minCladeSizeIdx;
    # print "\nUpdated imax: $imax\n";
    $ret = 1;
  }
  elsif ( !( $start[0] == 0 || $end[0] == $#leaves) )
  {
    warn "\n\n\tERROR: In the pruned tree outgroups sequences are not at the top or bottom of the tree!";

    print "\n\nNumber of leaves: " . @leaves . "\n";
    print "\nstart\tend\tsize\n";
    for my $i (0..$#start)
    {
      print "$start[$i]\t$end[$i]\t$rangeSize[$i]\n";
    }
    print "\n";

    print_array(\@og, "og");
    print "\n";

    my $maxOGbigSize = 100;
    if ( @ogBig < $maxOGbigSize )
    {
      print_array(\@ogBig, "Leaves from first to last OG seq");
    }

    print "\t--- Extracting the clade of OG sequences\n";
    my $ogCladeTreeFile = "$outDir/genus" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$outDir/genus" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug;
    my @ogCladeLeaves = read_array($ogCladeTreeLeavesFile);

    my $maxCladeSize = 100;
    if ( @ogCladeLeaves < $maxCladeSize )
    {
      print_array( \@ogCladeLeaves, "OG Clade Leaves" );
    }
    else
    {
      print "\n\tLeaves of the OG clade written to $ogCladeTreeLeavesFile\n"
    }

    print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
    print   "\tNumber of OG sequences: " . @og . "\n\n";

    $ret = 1;
  }
  else
  {
    print "\t--- Extracting the clade of OG sequences\n" if $debug;
    my $ogCladeTreeFile = "$outDir/genus" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$outDir/genus" . "_OG_clade.leaves";
    $cmd = "rm -f $ogCladeTreeLeavesFile; $nw_labels -I $ogCladeTreeFile > $ogCladeTreeLeavesFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Reading the leaves\n" if $debug;
    my @ogCladeLeaves = read_array( $ogCladeTreeLeavesFile );

    if ( @ogCladeLeaves != @og )
    {
      warn "\n\n\tERROR: The outgroup sequences do not form a monophyletic clade!";

      my $maxCladeSize = 100;
      if ( @ogCladeLeaves < $maxCladeSize )
      {
	print_array(\@ogCladeLeaves, "OG Clade Leaves");
      }
      else
      {
	print "\n\tLeaves of the OG clade written to $ogCladeTreeLeavesFile\n"
      }

      print "\n\tNumber of leaves of the OG clade: " . @ogCladeLeaves . "\n";
      print   "\tNumber of OG sequences: " . @og . "\n\n";
      $ret = 1;
    }
  }

  print "\n\tOG seq's form a monophylectic clade at the top or bottom of the tree\n\n" if $debug;

  return $ret;
}

sub argmax
{
  my $r = shift;

  my $index = undef;
  my $max   = undef;

  my $count = 0;
  foreach my $val(@{$r})
  {
    if ( not defined $max or $val > $max )
    {
      $max   = $val;
      $index = $count;
    }
    $count++;
  }
  return $index;
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

# common part of two arrays
sub comm
{
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

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl
{
  my $file = shift;

  if ( ! -e $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

# print elements of a hash table whose values are reference to a hash table so
sub print_tbl_valued_tbl
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

  for (@args)
  {
    print "$_\n";
    for my $e ( keys %{$rTbl->{$_}} )
    {
      print "\t$e\n";
    }
  }
  print "\n";
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

# write array to a file (one column format)
sub write_array
{
  my ($a, $outFile) = @_;
  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT "$_\n"} @{$a};
  close OUT
}

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

## Create a condensed tree given a tree with certain labels at the leaves
sub condense_tree_only
{
  my ($treeFile, $condTreeFile) = @_;

  $cmd = "rm -f $condTreeFile; $nw_condense $treeFile > $condTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
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

  if ( !setequal( \@txs, \@treeLeaves ) )
  {
    warn "\n\n\tERROR: Discrepancy between sequence IDs of the leaves of $treeFile and the keys of the taxon table";
    print "\n\n";
    exit 1;
  }

  my $sppTreeFile = "tmp_spp.tree";
  $cmd = "rm -f $sppTreeFile; $nw_rename $treeFile $txFile | $nw_order -c n  - > $sppTreeFile";
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

## are two arrays equal set-theoretically
sub setequal
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

    print "\na: ";
    map {print "$_ "} @a;
    print "\n\n";

    print "b: ";
    map {print "$_ "} @b;
    print "\n\n";

    # writeArray(\@a, "a.txt");
    # writeArray(\@b, "b.txt");
    #print "\n\tNew taxon keys and fasta IDs written to a.txt and b.txt, respectively\n\n";

    if (@a > @b)
    {
      my @d = diff(\@a, \@b);
      print "\nElements in a, but not b:\n";
      for (@d)
      {
	print "\t$_\n";
      }
      print "\n\n";
    }

    if (@b > @a)
    {
      my @d = diff(\@b, \@a);
      print "\nElements in b that are not in a:\n";
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

exit 0;
