#!/usr/bin/env perl

=head1 NAME

  vicut_and_plot.pl

=head1 DESCRIPTION

  Run vicut on a tree given an annotation file and generate the corresponding
  annotation tree, its condensed version and also a version with cluster IDs and
  the cluster sizes if an optional lineage file is supplied (taxons sizes).

=head1 SYNOPSIS

  vicut_and_plot.pl -t <tree file> -a <ann file> -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--li-file, -l>
  Master lineage file (seqID => lineage)

=item B<--ann-file, -a>
  Annotation file (seqID => labels)

=item B<--tree-file, -t>
  Consensus tree file

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

  cd ~/projects/PECAN/data/microcontax
  $ microcontax_to_banfield.pl -o contax_trim_banfield_itr_mothur_algn.fa
  $ FastTree -nt contax_trim_banfield_itr_mothur_algn.fa > contax_trim_banfield_itr_mothur.tree

  vicut_and_plot.pl -t contax_trim_banfield_itr_mothur.tree -a contax_trim_phylum.tx -l contax_trim.lineage -o contax_trim_banfield_phylum_dir

  vicut_and_plot.pl -t microcontax_mediods_V3V4_algn.tree -l medoids.lineage -o microcontax_mediods_V3V4_vicut_dir

  vicut_and_plot.pl -t medoids_ginsi.tree -l medoids.lineage -o mediods_ginsi_vicut_dir

  vicut_and_plot.pl -t Banfield_medoids_FL_algn.tree -l medoids.lineage -o banfield_medoids_FL_vicut_dir

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
my $spSizeThld  = 100;
GetOptions(
  #"tx-file|a=s"         => \my $txFile,
  "li-file|l=s"         => \my $liFile,
  "tree-file|t=s"       => \my $treeFile,
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
# elsif ( !$txFile )
# {
#   print "ERROR: Missing master lineage file\n\n";
#   pod2usage(verbose => 2,exitstatus => 0);
#   exit 1;
# }
elsif ( !$liFile )
{
  print "ERROR: Missing master fasta file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$treeFile )
{
  print "ERROR: Missing tree file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}


# if ( ! -e $txFile || ! -s $txFile )
# {
#   warn "ERROR: Master taxon file $txFile does not exist (or has size 0)";
#   print "\n\n";
#   exit 1;
# }

if ( ! -e $liFile || ! -s $liFile )
{
  warn "ERROR: Master fasta file does not exist (or has size 0)";
  print "\n\n";
  exit 1;
}

if ( ! -e $treeFile || ! -s $treeFile )
{
  warn "ERROR: Tree file does not exist (or has size 0)";
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

my $nProcStr = "";
if ( $nProc )
{
  $nProcStr = "--thread $nProc";
}

my $baseDir               = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
my $spGeFile              = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/species_genus_tbl_may19_2017.txt";
my $geLiFile              = "/Users/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this file has a header !
# my $ogFaFile              = "/Users/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4.fa";
# my $ogTxFile              = "/Users/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4_v2.tx";

my $txFile                = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/master_V3V4_no_outliers.tx";

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
  $ginsi                 = "/home/pgajer/bin/mafft --maxiterate 1000 --globalpair"; # MAFFT v7.310 (2017/Mar/17)

  $spGeFile              = "/usr/local/projects/pgajer/devel/MCextras/data/species_genus_tbl_may19_2017.txt";
  $geLiFile              = "/usr/local/projects/pgajer/devel/MCextras/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this one has a header !
  #$ogFaFile              = "/usr/local/projects/pgajer/devel/MCextras/data/phylo_split/Archaea_Caldococcus_noboribetus_S000414080.fa";

  $baseDir               = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
  $spGeFile              = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_groups/v0.2/species_genus_tbl_may19_2017.txt";
  $geLiFile              = "/usr/local/projects/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this file has a header !
  # $ogFaFile              = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4.fa";
  # $ogTxFile              = "/usr/local/projects/pgajer/projects/PECAN/data/phylo_split/bacterial_OGs/Archaea_OG_for_Bacterial_V3V4_v2.tx";

  $nProcStr              = "";
  $nProc                 = 1;
}

####################################################################
##                               MAIN
####################################################################

print "--- Generating our db taxon size tables\n";
my ($rsizeTbl, $rsizeTblR) = get_size_tbls();

my %sizeTbl  = %{$rsizeTbl};  # taxon => size
my %sizeTblR = %{$rsizeTblR}; # taxon => size with the uppler limit on the number of representatives of at the species level



## tree associated lineage file will be used to construct
## a tree with genus names at the leaves
## and then run vicut with annotation files
## mapping genera to higher taxonomic ranks

## for this I have build helper routines that build appropriate tables from the
## lineage file

print "--- Parsing contax lineage file\n";
my ($rgeTbl, $rgePhTbl, $rgeClTbl, $rgeOrTbl, $rgeFaTbl) = get_tbls_from_lineage( $liFile );

my %geTbl   = %{$rgeTbl};
my %gePhTbl = %{$rgePhTbl};
my %geClTbl = %{$rgeClTbl};
my %geOrTbl = %{$rgeOrTbl};
my %geFaTbl = %{$rgeFaTbl};

my $geTblFile = $outDir . "/genus.tx";
write_tbl( \%geTbl, $geTblFile );

print "--- Building genus names tree\n";
my $geTreeFile = $outDir . "/genus.tree";
build_ann_tree( $treeFile, $geTblFile, $geTreeFile );

$treeFile = $geTreeFile;
my @leaves = get_leaves( $treeFile );

##
## Phylum analysis
##
print "--- Phylum analysis\n";

print "--- Genus => Phylum table in the order of the genus tree leaves\n";
my $gePhFile = $outDir . "/genus_phylum.txt";
write_sorted_tbl( \%gePhTbl, \@leaves, $gePhFile );

print "--- Generating tree with <phylum name> labels at leaves\n";
my $phTreeFile = $outDir . "/phylum.tree";
build_ann_tree( $treeFile, $gePhFile, $phTreeFile );

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

print "--- Running vicut on genus tree using phylum annotation\n";
my $phVicutDir = $outDir . "/phylum_vicut_dir";
#run_vicut_no_query( $treeFile, $gePhFile, $phVicutDir );
run_vicut( $treeFile, $gePhFile, $phVicutDir );

print "--- Generating table of phylum-cluster-size labels\n";
my $phCltrFile = $phVicutDir . "/minNodeCut.cltrs";
my ($rphCltrSizeTbl, $rphCltrSizeTblR) = get_cltrSizes_labels( \%sizeTbl, \%sizeTblR, $phCltrFile );

# print "\n\nAfter get_cltrSizes_labels()\n";
# my $counter = 0;
# for ( keys %{$rphCltrSizeTbl} )
# {
#   print "$_\t" . $rphCltrSizeTbl->{$_} . "\n";
#   last if $counter == 10;
#   $counter++;
# }

# exit;

my $phCltrSizeTblFile = $outDir . "/genus_phVicut.cltrSize";
write_tbl( $rphCltrSizeTbl, $phCltrSizeTblFile );

my $phCltrSizeTblFileR = $outDir . "/genus_phVicut.cltrSizeR";
write_tbl( $rphCltrSizeTblR, $phCltrSizeTblFileR );


print "--- Generating tree with <cltr>_<size> labels\n";

my $condPhCltrSizeTreeFile = $outDir . "/cond_phCltrSize.tree";
condense_tree( $treeFile, $phCltrSizeTblFile, $condPhCltrSizeTreeFile );

my $condPhCltrSizeTreeFileR = $outDir . "/cond_phCltrSizeR.tree";
condense_tree( $treeFile, $phCltrSizeTblFileR, $condPhCltrSizeTreeFileR );


##
## Class analysis
##
print "--- Class analysis\n";

print "--- Genus => Class table in the order of the genus tree leaves\n";
print "\nGenus => Class on genus tree leaves\n";
my $geClFile = $outDir . "/genus_class.txt";
write_sorted_tbl( \%geClTbl, \@leaves, $geClFile );

print "--- Generating tree with <class name> labels at leaves\n";
my $clTreeFile = $outDir . "/class.tree";
build_ann_tree( $treeFile, $geClFile, $clTreeFile );

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
#run_vicut_no_query( $treeFile, $geClFile, $clVicutDir );
run_vicut( $treeFile, $geClFile, $clVicutDir );

print "--- Generating table of class-cluster-size labels\n";
my $clCltrFile = $clVicutDir . "/minNodeCut.cltrs";
my ($rclCltrSizeTbl, $rclCltrSizeTblR) = get_cltrSizes_labels( \%sizeTbl, \%sizeTblR, $clCltrFile );

my $clCltrSizeTblFile = $outDir . "/genus_clVicut.cltrSize";
my $clCltrSizeTblFileR = $outDir . "/genus_clVicut.cltrSizeR";

write_tbl( $rclCltrSizeTbl, $clCltrSizeTblFile );
write_tbl( $rclCltrSizeTblR, $clCltrSizeTblFileR );


print "--- Generating tree with <cltr>_<size> labels\n";

my $condClCltrSizeTreeFile = $outDir . "/cond_clCltrSize.tree";
condense_tree( $treeFile, $clCltrSizeTblFile, $condClCltrSizeTreeFile );

my $condClCltrSizeTreeFileR = $outDir . "/cond_clCltrSizeR.tree";
condense_tree( $treeFile, $clCltrSizeTblFileR, $condClCltrSizeTreeFileR );


##
## Order analysis
##
print "--- Order analysis\n";

print "--- Genus => Order table in the order of the genus tree leaves\n";
print "\nGenus => Order on genus tree leaves\n";
#print_formated_tbl( \%geOrTbl, \@leaves );
my $geOrFile = $outDir . "/genus_order.txt";
write_sorted_tbl( \%geOrTbl, \@leaves, $geOrFile );

print "--- Generating tree with <order name> labels at leaves\n";
my $orTreeFile = $outDir . "/order.tree";
build_ann_tree( $treeFile, $geOrFile, $orTreeFile );

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
#run_vicut_no_query( $treeFile, $geOrFile, $orVicutDir );w
run_vicut( $treeFile, $geOrFile, $orVicutDir );

print "--- Generating table of order-cluster-size labels\n";
my $orCltrFile = $orVicutDir . "/minNodeCut.cltrs";
my ($rorCltrSizeTbl, $rorCltrSizeTblR) = get_cltrSizes_labels( \%sizeTbl, \%sizeTblR, $orCltrFile );

my $orCltrSizeTblFile = $outDir . "/genus_orVicut.cltrSize";
my $orCltrSizeTblFileR = $outDir . "/genus_orVicut.cltrSizeR";

write_tbl( $rorCltrSizeTbl, $orCltrSizeTblFile );
write_tbl( $rorCltrSizeTblR, $orCltrSizeTblFileR );


print "--- Generating tree with <cltr>_<size> labels\n";

my $condOrCltrSizeTreeFile = $outDir . "/cond_orCltrSize.tree";
condense_tree( $treeFile, $orCltrSizeTblFile, $condOrCltrSizeTreeFile );

my $condOrCltrSizeTreeFileR = $outDir . "/cond_orCltrSizeR.tree";
condense_tree( $treeFile, $orCltrSizeTblFileR, $condOrCltrSizeTreeFileR );


##
## Family analysis
##
print "--- Family analysis\n";

print "--- Genus => Family table in the order of the genus tree leaves\n";
my $geFaFile = $outDir . "/genus_family.txt";
write_sorted_tbl( \%geFaTbl, \@leaves, $geFaFile );

print "--- Generating tree with <family name> labels at leaves\n";
my $faTreeFile = $outDir . "/family.tree";
build_ann_tree( $treeFile, $geFaFile, $faTreeFile );

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
#run_vicut_no_query( $treeFile, $geFaFile, $faVicutDir );
run_vicut( $treeFile, $geFaFile, $faVicutDir );


print "--- Generating table of family-cluster-size labels\n";
my $faCltrFile = $faVicutDir . "/minNodeCut.cltrs";
my ($rfaCltrSizeTbl, $rfaCltrSizeTblR) = get_cltrSizes_labels( \%sizeTbl, \%sizeTblR, $faCltrFile );

my $faCltrSizeTblFile = $outDir . "/genus_faVicut.cltrSize";
my $faCltrSizeTblFileR = $outDir . "/genus_faVicut.cltrSizeR";

write_tbl( $rfaCltrSizeTbl, $faCltrSizeTblFile );
write_tbl( $rfaCltrSizeTblR, $faCltrSizeTblFileR );


print "--- Generating tree with <cltr>_<size> labels\n";
my $condFaCltrSizeTreeFile = $outDir . "/cond_faCltrSize.tree";
condense_tree( $treeFile, $faCltrSizeTblFile, $condFaCltrSizeTreeFile );

my $condFaCltrSizeTreeFileR = $outDir . "/cond_faCltrSizeR.tree";
condense_tree( $treeFile, $faCltrSizeTblFileR, $condFaCltrSizeTreeFileR );


print "\n\n\tSuccessfully finished all tasks\n";
print "\n\n";

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

sub get_li_files
{
  my @liFiles0 = ("Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4.lineage",
		  "Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4.lineage",
		  "Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4.lineage",
		  "Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4.lineage",
		  "Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4.lineage",
		  "Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4.lineage",
		  "Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4.lineage",
		  "Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4.lineage",
		  "Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4.lineage",
		  "Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4.lineage",
		  "Chloroflexi_V3V4_dir/Chloroflexi_V3V4.lineage",
		  "Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4.lineage",
		  "Fusobacteria_V3V4_dir/Fusobacteria_V3V4.lineage",
		  "Nitrospirae_V3V4_dir/Nitrospirae_V3V4.lineage",
		  "Planctomycetes_V3V4_dir/Planctomycetes_V3V4.lineage",
		  "Spirochaetes_V3V4_dir/Spirochaetes_V3V4.lineage",
		  "Tenericutes_V3V4_dir/Tenericutes_V3V4.lineage",
		  "Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4.lineage",
		  "phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4.lineage",
		  "Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4.lineage",
		  "Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4.lineage",
		  "Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4.lineage",
		  "Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4.lineage",
		  "Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4.lineage",
		  "Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4.lineage",
		  "Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4.lineage",
		  "Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4.lineage",
		  "Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4.lineage",
		  "Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4.lineage",
		  "Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4.lineage",
		  "Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4.lineage",
		  "Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4.lineage",
		  "Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4.lineage",
		  "Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4.lineage",
		  "Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4.lineage",
		  "Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4.lineage",
		  "Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4.lineage",
		  "Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4.lineage",
		  "Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4.lineage",
		  "Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4.lineage",
		  "Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4.lineage",
		  "Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4.lineage",
		  "Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4.lineage");

  my @liFiles = map{ $baseDir . $_ } @liFiles0;

  ## testing for existence
  for my $file ( @liFiles )
  {
    if ( ! -e $file )
    {
      warn "\n\n\tERROR: $file does not exist";
      print "\n\n";
      exit 1;
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

  ## get leaves and test if they have the same number of elements as the table or
  ## the table has those plus more


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
  #map {print OUT $_ . "\t" . $tbl{$_} . "\n"} @{$r};
  for ( @{$r} )
  {
    if ( exists $tbl{$_} )
    {
      print OUT $_ . "\t" . $tbl{$_} . "\n"
    }
  }
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
sub diff
{
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

sub run_vicut
{
  my ($treeFile, $annFile, $vicutDir) = @_;

  my %annTbl = read_tbl( $annFile );
  my @ann = keys %annTbl;
  my @leaves = get_leaves( $treeFile );

  my @d = diff( \@leaves, \@ann );
  my ($fh, $qFile) = tempfile("query.XXXX", SUFFIX => 'txt', OPEN => 1, DIR => $tmpDir);
  for ( @d )
  {
    print $fh "$_\n";
  }
  close $fh;

  my $cmd = "$vicut $quietStr -t $treeFile -a $annFile -q $qFile -o $vicutDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
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

sub get_cltrSizes_labels
{
  my ( $rsizeTbl, $rsizeTblR, $txCltrFile ) = @_;

  my %cltrTbl  = read_tbl( $txCltrFile );

  my %sizeTbl  = %{ $rsizeTbl };
  my %sizeTblR = %{ $rsizeTblR };

  my %clSize;
  my %clSizeR;

  for my $tx ( keys %cltrTbl )
  {
    my $cl = "c" . $cltrTbl{$tx};
    if ( ! exists $sizeTbl{ $tx } )
    {
      $sizeTbl{ $tx }++;
    }
    if ( ! exists $sizeTblR{ $tx } )
    {
      $sizeTblR{ $tx }++;
    }
    $clSize{ $cl }  += $sizeTbl{ $tx };
    $clSizeR{ $cl } += $sizeTblR{ $tx };
  }

  my %cltrSizeTbl;
  my %cltrSizeTblR;

  for my $tx ( keys %cltrTbl )
  {
    my $cl = "c" . $cltrTbl{$tx};
    $cltrSizeTbl{$tx} = $cl . "__" . $clSize{ $cl };
    $cltrSizeTblR{$tx} = $cl . "__" . $clSizeR{ $cl };
  }

  return (\%cltrSizeTbl, \%cltrSizeTblR);
}


## constructing tables reporting the number of sequences we have in our db at
## each taxonomic rank
sub get_size_tbls
{
  print "--- Parsing seqID => species table\n";
  my %spTbl = read_tbl( $txFile );

  print "--- Parsing species-genus tbl\n";
  my %spGeTbl = read_tbl( $spGeFile );


  print "--- Parsing microcontax pkg genus-lineage tbl\n";
  my ($rgeFaTbl, $rgeOrTbl, $rgeClTbl, $rgePhTbl) = ge_li_tbl( $geLiFile );

  my %geFaTbl = %{$rgeFaTbl};
  my %geOrTbl = %{$rgeOrTbl};
  my %geClTbl = %{$rgeClTbl};
  my %gePhTbl = %{$rgePhTbl};

  # for ( keys %ogTx )
  # {
  #   my $tx = $ogTx{$_};
  #   $geFaTbl{$tx} = "OG";
  #   $geOrTbl{$tx} = "OG";
  #   $geClTbl{$tx} = "OG";
  #   $gePhTbl{$tx} = "OG";
  # }

  print "--- Building species-genus tables\n";
  my %spSeqIDs;
  my %geFileFreq;
  my %geSpFreq;
  my %liFile2phGr;
  my $spGeTblNeedsUpdate = 0;
  foreach my $id ( keys %spTbl )
  {
    my $sp = $spTbl{$id};

    if ( ! exists $spGeTbl{$sp} )
    {
      my ($g, $suffix) = split "_", $sp;
      print "WARNING: $sp not found in spGeTbl - using $g\n";
      $spGeTbl{$sp} = $g;
      $spGeTblNeedsUpdate = 1;
      #exit 1;
    }
    my $ge = $spGeTbl{$sp};
    $geSpFreq{$ge}{$sp}++;
    push @{ $spSeqIDs{$sp} }, $id;
  }


  if ( $spGeTblNeedsUpdate )
  {
    print "--- Updating $spGeFile\n";
    write_tbl( \%spGeTbl, $spGeFile );
  }

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


  print "--- Computing taxon sizes (with and without capping of species sizes)\n";

  my %sizeTbl;  # taxon => size
  my %sizeTblR; # taxon => size

  for my $ge ( keys %geSpFreq )
  {
    print "\rProcessing $ge                                      ";
    my %spFreq = %{ $geSpFreq{$ge} };
    my @spp = keys %spFreq;
    my $geSize = 0;
    map { $geSize += $spFreq{$_} } @spp;
    $sizeTbl{$ge} = $geSize;

    for ( @spp )
    {
      if ( $spFreq{$_} > $spSizeThld )
      {
	$spFreq{$_} = $spSizeThld;
      }
    }

    $geSize = 0;
    map { $geSize += $spFreq{$_} } @spp;
    $sizeTblR{$ge} = $geSize;
  }

  return (\%sizeTbl, \%sizeTblR);
}


sub get_tbls_from_lineage
{
  my $file = shift;

  my %geTbl;   # seqID => genus
  my %geFaTbl; # genus => family
  my %geOrTbl; # genus => order
  my %geClTbl; # genus => class
  my %gePhTbl; # genus => phylum

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    chomp;
    my @f = split /\s+/;

    my $id = shift @f;
    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;
    #my $do = shift @f;

    $geTbl{$id}   = $ge;
    $geFaTbl{$ge} = $fa;
    $geOrTbl{$ge} = $or;
    $geClTbl{$ge} = $cl;
    $gePhTbl{$ge} = $ph;
  }
  close IN;

  return (\%geTbl, \%gePhTbl, \%geClTbl, \%geOrTbl, \%geFaTbl);
}

exit 0;
