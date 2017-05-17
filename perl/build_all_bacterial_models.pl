#!/usr/bin/env perl

=head1 NAME

  build_all_bacterial_models.pl

=head1 DESCRIPTION

  Generating MC models for all bacterial species using 16S rRNA sequences
  restricted to the V3V4 amplicon.

  Data from 43 phyl-groups is being used here to build one classifier.

=head1 SYNOPSIS

  build_all_bacterial_models.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

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

  build_all_bacterial_models.pl

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $offsetCoef   = 0.9;
my $txSizeThld   = 10;

GetOptions(
  "offset-coef|o=f"     => \$offsetCoef,
  "tx-size-thld|t=i"    => \$txSizeThld,
  "pp-embedding"        => \my $ppEmbedding,
  "skip-err-thld"       => \my $skipErrThld,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

####################################################################
##                               MAIN
####################################################################

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

my $skipErrThldStr = "";
if ($skipErrThld)
{
  $skipErrThldStr = "--skip-err-thld";
}

my $ppEmbeddingStr = "";
if ($ppEmbedding)
{
  $ppEmbeddingStr = "--pp-embedding"
}

## Gather final lineage, fasta (ungapped) and taxon files for all phylo-groups


my $baseDirExt = "/Users/pgajer/devel/MCextras/data/vag_exp_V3V4_phGrps_May16_dir/";

my @faFilesExt0 = ("Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.fa",
		   "Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.fa",
		   "Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.fa",
		   "Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.fa",
		   "Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.fa",
		   "Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.fa",
		   "Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.fa",
		   "Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.fa",
		   "Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.fa",
		   "Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.fa",
		   "Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.fa",
		   "Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.fa",
		   "phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.fa",
		   "Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.fa",
		   "Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.fa",
		   "Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.fa",
		   "Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.fa",
		   "Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.fa",
		   "Tenericutes_V3V4_dir/Tenericutes_V3V4_final.fa");

my @faFilesExt1 = map{ $_ = $baseDirExt . $_ } @faFilesExt0;



my @txFilesExt0 = ("Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.tx",
		   "Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.tx",
		   "Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.tx",
		   "Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.tx",
		   "Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.tx",
		   "Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.tx",
		   "Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.tx",
		   "Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.tx",
		   "Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.tx",
		   "Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.tx",
		   "Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.tx",
		   "Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.tx",
		   "phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.tx",
		   "Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.tx",
		   "Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.tx",
		   "Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.tx",
		   "Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.tx",
		   "Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.tx",
		   "Tenericutes_V3V4_dir/Tenericutes_V3V4_final.tx");

my @txFilesExt1 = map{ $_ = $baseDirExt . $_ } @txFilesExt0;

my @spLiFilesExt0 = ("Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.spLineage",
			  "Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.spLineage",
			  "Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.spLineage",
			  "Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.spLineage",
			  "Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.spLineage",
			  "Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.spLineage",
			  "Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.spLineage",
			  "Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.spLineage",
			  "Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.spLineage",
			  "Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.spLineage",
			  "Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.spLineage",
			  "Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.spLineage",
			  "phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.spLineage",
			  "Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.spLineage",
			  "Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.spLineage",
			  "Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.spLineage",
			  "Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.spLineage",
			  "Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.spLineage",
			  "Tenericutes_V3V4_dir/Tenericutes_V3V4_final.spLineage");

my @spLiFilesExt1 = map{ $_ = $baseDirExt . $_ } @spLiFilesExt0;


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

my %faExtTbl;
my %txExtTbl;
my %spLiExtTbl;
for my $i ( 0..$#extPhGrs )
{
  $faExtTbl{$extPhGrs[$i]}   = $faFilesExt1[$i];
  $txExtTbl{$extPhGrs[$i]}   = $txFilesExt1[$i];
  $spLiExtTbl{$extPhGrs[$i]} = $spLiFilesExt1[$i];
}


my $baseDir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";

my @faFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.fa",
		"Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.fa",
		"Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.fa",
		"Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4_final.fa",
		"Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4_final.fa",
		"Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4_final.fa",
		"Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4_final.fa",
		"Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4_final.fa",
		"Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.fa",
		"Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4_final.fa",
		"final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4_final.fa",
		"final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4_final.fa",
		"final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.fa",
		"final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4_final.fa",
		"final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4_final.fa",
		"final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4_final.fa",
		"final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4_final.fa",
		"final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4_final.fa",
		"final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.fa",
		"Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4_final.fa",
		"Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.fa");

my @faFiles1 = map{ $_ = $baseDir . $_ } @faFiles0;
## print "faFiles: @faFiles\n";

my @spLiFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.spLineage",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.spLineage",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.spLineage",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4_final.spLineage",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4_final.spLineage",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4_final.spLineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4_final.spLineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4_final.spLineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.spLineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4_final.spLineage",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.spLineage",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4_final.spLineage",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.spLineage");

my @spLiFiles1 = map{ $_ = $baseDir . $_ } @spLiFiles0;
## print "spLineageFile: @spLiFiles\n";


my @txFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.tx",
		"Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.tx",
		"Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.tx",
		"Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4_final.tx",
		"Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4_final.tx",
		"Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4_final.tx",
		"Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4_final.tx",
		"Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4_final.tx",
		"Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.tx",
		"Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4_final.tx",
		"final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4_final.tx",
		"final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4_final.tx",
		"final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.tx",
		"final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4_final.tx",
		"final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4_final.tx",
		"final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4_final.tx",
		"final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4_final.tx",
		"final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4_final.tx",
		"final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.tx",
		"Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4_final.tx",
		"Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.tx");

my @txFiles1 = map{ $_ = $baseDir . $_ } @txFiles0;
## print "txFile: @txFiles\n";

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

my %fastaTbl;
my %txTbl;
my %spLiTbl;
for my $i ( 0..$#allPhGrs )
{
  $fastaTbl{$allPhGrs[$i]}   = $faFiles1[$i];
  $txTbl{$allPhGrs[$i]}   = $txFiles1[$i];
  $spLiTbl{$allPhGrs[$i]} = $spLiFiles1[$i];
}

my @faFiles;
my @txFiles;
my @spLiFiles;
for ( @allPhGrs )
{
  if ( exists $faExtTbl{$_} )
  {
    push @faFiles,   $faExtTbl{$_};
    push @txFiles,   $txExtTbl{$_};
    push @spLiFiles, $spLiExtTbl{$_};
  }
  else
  {
    push @faFiles,   $fastaTbl{$_};
    push @txFiles,   $txTbl{$_};
    push @spLiFiles, $spLiTbl{$_};
  }
}


if ( $debug )
{
  print_array( \@faFiles, "\nfaFiles" );
  print "\n\n";

  print_array( \@txFiles, "\ntxFiles" );
  print "\n\n";

  print_array( \@spLiFiles, "\nspLiFiles" );
  print "\n\n";

  exit;
}


## checking that each taxon and fasta files have the same sequence IDs

## checking that all unique taxons from the taxon files are represented in the species spLineage file

##
## concatenating all fasta files
##

print "--- Concatenating all fasta files\n";
my $faFile = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/all_bacteria_V3V4.fa";

my $cmd = "rm -f $faFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# copying the first file in the list to outFile
my $file = shift @faFiles;
$cmd = "cp $file $faFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
# concatenating the rest of files to outFile
foreach $file (@faFiles)
{
  print "\rProcessing $file                  ";
  $cmd = "cat $file >> $faFile";
  print "cmd=$cmd\n" if $dryRun;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

##
## concatenating all taxon files
##

print "\r--- Concatenating all taxon files                           \n";
my $txFile = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/all_bacteria_V3V4.tx";

$cmd = "rm -f $txFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

# copying the first file in the list to outFile
$file = shift @txFiles;
$cmd = "cp $file $txFile";
print "cmd=$cmd\n" if $dryRun;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
# concatenating the rest of files to outFile
foreach $file (@txFiles)
{
  print "\rProcessing $file            ";
  $cmd = "cat $file >> $txFile";
  print "cmd=$cmd\n" if $dryRun;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

##
## concatenating all spLineage files
##

print "\r--- Concatenating all spLineage files                              \n";

my %spTbl;
my %sgeTbl;
my %geTbl;
my %faTbl;
my %orTbl;
my %clTbl;
my %phTbl;

##my %spLineage; # species => lineage of the species (recorded as a string) with the corresponding phyo-group name at the end
foreach $file (@spLiFiles)
{
  print "\rProcessing $file        ";
  my @tokens = split "/", $file;
  my $phGr = $tokens[7];
  $phGr =~ s/_dir//;
  #print "tokens: @tokens\nphGr: $phGr\n"; exit;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for my $lineage (<IN>)
  {
    chomp $lineage;
    my @f = split "\t", $lineage;
    my $sp = shift @f;
    my $sge = shift @f;
    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;

    ##$spLineage{$sp} = $lineage . "\t$phGr";

    # print "lineage: $lineage\n";
    # print "sp : $sp\n";
    # print "sge : $sge\n";
    # print "ge : $ge\n";
    # print "fa : $fa\n";
    # print "or : $or\n";
    # print "cl : $cl\n";
    # print "ph : $ph\n";
    # exit;
    ##print "\nExtended lineage: $spLineage{$sp}\n";  exit;

    $spTbl{$sp}{$phGr}   = 1;
    $sgeTbl{$sge}{$phGr} = 1;
    $geTbl{$ge}{$phGr}   = 1;
    $faTbl{$fa}{$phGr}   = 1;
    $orTbl{$or}{$phGr}   = 1;
    $clTbl{$cl}{$phGr}   = 1;
    $phTbl{$ph}{$phGr}   = 1;
  }
  close IN;
}

if (0)
{
  print "\rTesting if there are species with the same name in different phylo-groups. There should not be any";
  for my $sp (keys %spTbl)
  {
    my @gps = keys %{$spTbl{$sp}};
    if ( @gps > 1 )
    {
      print "$sp: @gps\n";
    }
  }

  print "\rTesting if there are sub-genera with the same name in different phylo-groups. There should not be any";
  for my $sge (keys %sgeTbl)
  {
    my @gps = keys %{$sgeTbl{$sge}};
    if ( @gps > 1 )
    {
      print "$sge: @gps\n";
    }
  }

  print "\rTesting if there are genera with the same name in different phylo-groups. There should not be any";
  for my $ge (keys %geTbl)
  {
    my @gps = keys %{$geTbl{$ge}};
    if ( @gps > 1 )
    {
      print "$ge: @gps\n";
    }
  }

  print "\rTesting if there are families with the same name in different phylo-groups.                            ";
  for my $fa (keys %faTbl)
  {
    my @gps = keys %{$faTbl{$fa}};
    if ( @gps > 1 )
    {
      print "$fa: @gps\n";
    }
  }

  print "\rTesting if there are orders with the same name in different phylo-groups.                               ";
  for my $or (keys %orTbl)
  {
    my @gps = keys %{$orTbl{$or}};
    if ( @gps > 1 )
    {
      print "$or: @gps\n";
    }
  }

  print "\rTesting if there are classes with the same name in different phylo-groups.                               ";
  for my $cl (keys %clTbl)
  {
    my @gps = keys %{$clTbl{$cl}};
    if ( @gps > 1 )
    {
      print "$cl: @gps\n";
    }
  }

  print "\rTesting if there are phyla with the same name in different phylo-groups.                                  ";
  for my $ph (keys %phTbl)
  {
    my @gps = keys %{$phTbl{$ph}};
    if ( @gps > 1 )
    {
      print "$ph: @gps\n";
    }
  }
}

##print "\r                                                                                                            \n\n";

# Testing if there are orders with the same name in different phylo-groups.
# o_Actinomycetales_2: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales: Actinobacteria_group_4_V3V4 Actinobacteria_group_2_V3V4 Actinobacteria_group_3_V3V4
# o_Actinomycetales_1: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_23: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_12: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_22: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_5: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_17: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_7: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_8: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_6: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Pseudomonadales: Proteobacteria_group_17_V3V4 Proteobacteria_group_15_V3V4
# o_Actinomycetales_24: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_4: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_18: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_21: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_13: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_16: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_9: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_3: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_11: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_20: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4
# o_Actinomycetales_15: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Lactobacillales: Firmicutes_group_5_V3V4 Firmicutes_group_6_V3V4 Firmicutes_group_4_V3V4
# o_Actinomycetales_14: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Burkholderiales: Proteobacteria_group_4_V3V4 Proteobacteria_group_3_V3V4
# o_Actinomycetales_10: Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4
# o_Actinomycetales_19: Actinobacteria_group_1_V3V4 Actinobacteria_group_5_V3V4

# Testing if there are classes with the same name in different phylo-groups.
# c_Alphaproteobacteria_5: Proteobacteria_group_5_V3V4 Proteobacteria_group_8_V3V4
# c_Actinobacteria: Actinobacteria_group_1_V3V4 Actinobacteria_group_3_V3V4 Actinobacteria_group_2_V3V4 Actinobacteria_group_0_V3V4 Actinobacteria_group_4_V3V4
# c_Betaproteobacteria: Proteobacteria_group_4_V3V4 Proteobacteria_group_2_V3V4 Proteobacteria_group_3_V3V4
# c_Bacilli: Firmicutes_group_5_V3V4 Firmicutes_group_3_V3V4 Firmicutes_group_4_V3V4 Firmicutes_group_6_V3V4
# c_Alphaproteobacteria_4: Proteobacteria_group_5_V3V4 Proteobacteria_group_8_V3V4
# c_Alphaproteobacteria_3: Proteobacteria_group_5_V3V4 Proteobacteria_group_8_V3V4
# c_Clostridia: Firmicutes_group_2_V3V4 Firmicutes_group_1_V3V4
# c_Alphaproteobacteria: Proteobacteria_group_6_V3V4 Proteobacteria_group_7_V3V4 Proteobacteria_group_9_V3V4
# c_Alphaproteobacteria_1: Proteobacteria_group_8_V3V4 Proteobacteria_group_5_V3V4
# c_Alphaproteobacteria_2: Proteobacteria_group_8_V3V4 Proteobacteria_group_5_V3V4
# c_Gammaproteobacteria: Proteobacteria_group_14_V3V4 Proteobacteria_group_12_V3V4 Proteobacteria_group_15_V3V4 Proteobacteria_group_11_V3V4 Proteobacteria_group_17_V3V4 Proteobacteria_group_13_V3V4

# Testing if there are phyla with the same name in different phylo-groups.
# p_Bacteroidetes: Bacteroidetes_group_3_V3V4 Bacteroidetes_group_2_V3V4 Bacteroidetes_group_0_V3V4 Bacteroidetes_group_1_V3V4
# p_Actinobacteria: Actinobacteria_group_0_V3V4 Actinobacteria_group_2_V3V4 Actinobacteria_group_4_V3V4 Actinobacteria_group_5_V3V4 Actinobacteria_group_1_V3V4 Actinobacteria_group_3_V3V4
# p_Proteobacteria: Proteobacteria_group_13_V3V4 Proteobacteria_group_1_V3V4 Proteobacteria_group_5_V3V4 Proteobacteria_group_4_V3V4 Proteobacteria_group_9_V3V4 Proteobacteria_group_2_V3V4 Proteobacteria_group_0_V3V4 Proteobacteria_group_12_V3V4 Proteobacteria_group_7_V3V4 Proteobacteria_group_15_V3V4 Proteobacteria_group_11_V3V4 Proteobacteria_group_8_V3V4 Proteobacteria_group_14_V3V4 Proteobacteria_group_10_V3V4 Proteobacteria_group_6_V3V4 Proteobacteria_group_3_V3V4 Proteobacteria_group_17_V3V4
# p_Firmicutes: Firmicutes_group_4_V3V4 Firmicutes_group_6_V3V4 Firmicutes_group_2_V3V4 Firmicutes_group_1_V3V4 Firmicutes_group_0_V3V4 Firmicutes_group_3_V3V4 Firmicutes_group_5_V3V4


## use the following tables to figure out if order, class or phylum appears in
## multiple phylo-groups. Only orders, classes and phyla with numerical suffix
## '_\d+' are to be made unique by attaching to them a suffix consisting of the
## first letter of the phylogroup and the group index. For example, if the taxon
## is present in Proteobacteria_group_9_V3V4 (among other groups) it will have
## the suffix _P9 attached to it.

##  %spTbl;
##  %sgeTbl;
##  %geTbl;
##  %faTbl;
##  %orTbl;
##  %clTbl;
##  %phTbl;

my %spLineage; # species => lineage of the species (recorded as a string) with the corresponding phyo-group name at the end
foreach $file (@spLiFiles)
{
  my @tokens = split "/", $file;
  my $phGr = $tokens[7];
  $phGr =~ s/_dir//;
  print "\rProcessing $phGr     ";
  my $phGrID;
  if ( $phGr =~ /_group/ )
  {
    my $phGrFirstLetter = substr($phGr, 0, 1);
    my @phGrTokens = split "_", $phGr;
    if (!defined $phGrFirstLetter)
    {
      warn "\n\n\tWARNING: Undefined phGrFirstLetter for phGr: $phGr";
      print "\n\n";
      exit 1;
    }
    if (!defined $phGrTokens[2])
    {
      warn "\n\n\tWARNING: Undefined phGrTokens[2] for phGr: $phGr";
      print "\n\n";
      exit 1;
    }
    $phGrID = "_$phGrFirstLetter" . $phGrTokens[2];
  }
  else
  {
    $phGrID = "_" . substr($phGr, 0, 2);
  }

  ##print "tokens: @tokens\t"; print "phGr: $phGr\n"; exit;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for my $lineage (<IN>)
  {
    chomp $lineage;
    my @f = split "\t", $lineage;
    my $sp = shift @f;
    my $sge = shift @f;
    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;

    my @gps = keys %{$orTbl{$or}};
    if ( @gps > 1 )
    {
      $or .= $phGrID;
      ##print "modified order: $or\n"; exit;
    }

    @gps = keys %{$clTbl{$cl}};
    if ( @gps > 1 )
    {
      $cl .= $phGrID;
    }

    @gps = keys %{$phTbl{$ph}};
    if ( @gps > 1 )
    {
      $ph .= $phGrID;
    }

    $spLineage{$sp} = "$sp\t$sge\t$ge\t$fa\t$or\t$cl\t$ph\td_Bacteria";
  }
  close IN;
}

my $spLineageFile = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/all_bacteria_V3V4.spLineage";
open OUT, ">$spLineageFile" or die "Cannot open $spLineageFile for writing: $OS_ERROR";
for my $sp (keys %spLineage)
{
  print OUT $spLineage{$sp} . "\n";
}
close OUT;

## Identifying taxonomic ranks with the same name that appear in different phylo-groups
## and resolve the duplicated attaching phylo-group identifier

print "\r--- Building model tree and creating taxon's reference fasta files       ";
my $mcDir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/all_bacteria_V3V4_MC_models_dir";
$cmd = "rm -rf $mcDir; buildModelTree $quietStr -l $spLineageFile -i $faFile -t $txFile -o $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "\r--- Building MC models                                                    ";
$cmd = "buildMC -t $mcDir/spp_paths.txt -k 8 -d $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "\r--- Estimating error thresholds                                      ";
$cmd = "est_error_thlds --offset-coef $offsetCoef --tx-size-thld $txSizeThld -d $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "\r--- Running classify on the ref fasta file                               ";
$cmd = "classify $skipErrThldStr -d $mcDir -i $faFile -o $mcDir"; # $ppEmbeddingStr
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "\r--- Comparing ref seq's taxonomy with the classification results                ";
$cmd = "cmp_tx.pl --verbose $quietStr -i $txFile -j $mcDir/MC_order7_results.txt -o $mcDir";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

print "\r                                                                                 ";
print "\n\t fasta file written to $faFile\n";
print "\t taxon file written to $txFile\n";
print "\t spLineage written to $spLineageFile\n";
print "\t MC master models written to $mcDir\n\n";


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

exit 0;
