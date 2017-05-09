#!/usr/bin/env perl

=head1 NAME

  mm_validate_pecan.pl

=head1 DESCRIPTION

  A version of mm_validate_pecan.pl that instead of processing all species in a
  one big loop, processes only selected species passed to the script in a file.

  This script attempts to validate PECAN taxonomic assignment of the M&M
  sequences by generating vicut clustering results report and creating
  phylogenetic trees of representative sequences to PECAN species together with
  ref seq's of of the corresponding phylo group.

=head1 SYNOPSIS

  mm_validate_pecan.pl

=head1 OPTIONS

=over

=item B<--spp-file, -i>
  Maximal number of non-redundant seq's.

=item B<--out-dir, -o>
  Output directory. Optional parameter. If not specified the output is written to
  /Users/pgajer/projects/M_and_M/new_16S_classification_data/mm_validate_reports_dir_<timeStamp>.

=item B<--max-no-nr-seqs, -n>
  Maximal number of non-redundant seq's.

=item B<--perc-coverage, -p>
  Percentage coverage: the percentage, say 80%, of the total number of non-redundant sequences
  Used to reduce the number of non-redundant seq's

=item B<--phylo-part-perc-thld, -t>
  Percentile threshold for phylo-partitioning specified as a decimal between 0 and 1. Default: 0.1

=item B<--num-proc, -m>
  Number of processors to be used. Default value: 8.

=item B<--run-all>
  Ignore if ( ! -e ... ) statements.

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

  in ~/projects/M_and_M/new_16S_classification_data

  mm_validate_pecan.pl --spp-file mm_valid_spp_part_1.txt --max-no-nr-seqs 500 --perc-coverage 80

  mm_validate_pecan.pl --num-proc 8 --spp-file spp_missing_record.txt --max-no-nr-seqs 500 --perc-coverage 80

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
use Parallel::ForkManager;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $nProc             = 1;
my $maxNumNRseqs      = 100;
my $percCoverage      = 80;
my $phyloPartPercThld = 0.1;

GetOptions(
  "spp-file|i=s"             => \my $sppFile,
  "out-dir|o=s"              => \my $outDir,
  "max-no-nr-seqs|n=i"       => \$maxNumNRseqs,
  "perc-coverage|p=i"        => \$percCoverage,
  "num-proc|m=i"             => \$nProc,
  "igs"                      => \my $igs,
  "run-all"                  => \my $runAll,
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

if ( !$sppFile )
{
  warn "\n\n\tERROR: Missing species list file";
  print "\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $baseDir        = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";
my $mmDir          = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/";
my $mmMasterFaFile = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mm_all_good_seqIDs.fa";
my $mmPECANfile    = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mmDir_May5/mm_pecan_may5_no_controls.tx";

my $mothur         = "/Users/pgajer/bin/mothur";
my $usearch6       = "/Users/pgajer/bin/usearch6.0.203_i86osx32";
my $readNewickFile = "/Users/pgajer/.Rlocal/read.newick.R";

if ( defined $igs )
{
  $baseDir        = "/local/scratch/MM/V3V4_unzipped/";
  $mmDir          = "/local/scratch/MM/";
  $mmMasterFaFile = "/local/scratch/MM/MM_all_updated_seqID_no_ctrls.fa";
  $mmPECANfile    = "/local/scratch/MM/MM_MC_order7_results.txt";
  $mothur         = "/usr/local/packages/mothur-1.36.1/mothur";
  $usearch6       = "/local/projects/pgajer/bin/usearch6.0.203_i86linux32";
  $readNewickFile = "/home/pgajer/.Rlocal/read.newick.R";
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

my @faFiles = map{ $_ = $baseDir . $_ } @faFiles0;
## print "faFiles: @faFiles\n";

my @algnFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_algn_trimmed_final.fa",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_algn_trimmed_final.fa",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_algn_trimmed_final.fa",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4_algn_trimmed_final.fa",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4_algn_trimmed_final.fa",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4_algn_trimmed_final.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4_algn_trimmed_final.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4_algn_trimmed_final.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_algn_trimmed_final.fa",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4_algn_trimmed_final.fa",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_algn_trimmed_final.fa",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4_algn_trimmed_final.fa",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_algn_trimmed_final.fa");

my @algnFiles = map{ $_ = $baseDir . $_ } @algnFiles0;
## print "algnFiles: @algnFiles\n";


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

my @txFiles = map{ $baseDir . $_ } @txFiles0;
# print "txFile: @txFiles\n"; exit;


##
## Creating output reports directory
##

if ( !defined $outDir )
{
  # my @now = localtime();
  # my $timeStamp = sprintf("%04d%02d%02d_%02d%02d%02d",
  # 			  $now[5]+1900, $now[4]+1, $now[3],
  # 			  $now[2],      $now[1],   $now[0]);
  # $outDir = $mmDir . "mm_validate_reports_dir_$timeStamp";
  $outDir = $mmDir . "mm_validate_reports_dir";
}

if ( ! -e $outDir )
{
  my $cmd = "mkdir -p $outDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?";# if !$dryRun;
}

my $treesDir = $outDir . "/trees_dir";

if ( ! -e $treesDir )
{
  my $cmd = "mkdir -p $treesDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}

my $tmpDir = $outDir . "/temp_dir";
if ( ! -e $tmpDir )
{
  my $cmd = "mkdir -p $tmpDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}


##
## MAIN section
##

print "--- Parsing table of species to be processed\n";
my %phGrSppTbl = parse_spp_tbl( $sppFile ); # phGr => ref to array of species from that phylo-group

print "--- Parsing PECAN classification results on M&M's sequences\n";
my ($rspIDsTbl, $rppTbl) = parse_pecan_tbl( $mmPECANfile );

my %spIDsTbl = %{ $rspIDsTbl };   # sp   => ref of array with seqIDs of seq's classified to sp
my %ppTbl    = %{ $rppTbl };      # seqID => posterior probability of the best model

my $forkMgr = new Parallel::ForkManager( $nProc );

for my $phGr ( keys %phGrSppTbl )
{
  print "\r--- Processing $phGr species                                    \n";

  # print $ROUT "================================================\n";
  # print $ROUT "\t$phGr\n";
  # print $ROUT "================================================\n\n";

  ## create mm phylo-group dir
  my $phGrDir = $outDir . "/mm_" . $phGr . "_dir/";
  ## print "\n\nphGrDir: $phGrDir\n"; exit;

  ##my $cmd = "rm -rf $phGrDir; mkdir $phGrDir";
  my $cmd = "mkdir -p $phGrDir";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  ## Identifying fa finale file of the given phylo-group
  my @f = grep { $_ =~ /$phGr/ } @faFiles;
  my $phGrFaFile = $f[0];
  ## print "\nphGr: $phGr; phGrFaFile: $phGrFaFile\n";

  ## Identifying algn file of the given phylo-group
  @f = grep { $_ =~ /$phGr/ } @algnFiles;
  my $phGrAlgnFile = $f[0];

  if ( -l $phGrAlgnFile )
  {
    $phGrAlgnFile = readlink( $phGrAlgnFile );
  }
  ## print "\nphGr: $phGr; phGrAlgnFile: $phGrAlgnFile\n";
  ## Final alignment has OG seq's !!!!

  ## Identifying tx finale file of the given phylo-group
  @f = grep { $_ =~ /$phGr/ } @txFiles;
  my $phGrTxFile = $f[0];
  ## print "\nphGr: $phGr; phGrTxFile: $phGrTxFile\n";

  #my %phGrTxTbl = readTbl($phGrTxFile);

  ## file with the given phylo-group's outgroup seq's
  my $phGrOGseqIDsFile = $phGrFaFile;
  $phGrOGseqIDsFile =~ s/_final\.fa/_outgroup\.seqIDs/;
  ## print "phGrOGseqIDsFile: $phGrOGseqIDsFile\n";

  my @ogSeqIDs = readArray($phGrOGseqIDsFile);
  print "\nNo. OG seq's: " . scalar(@ogSeqIDs) . "\n";

  ## file with the given phylo-group's fa file of all seq's before curation including outgroup seq's
  my $phGrBigFaFile = $phGrFaFile;
  $phGrBigFaFile =~ s/_final//;
  ## print "phGrBigFaFile: $phGrBigFaFile\n";

  my $phGrOGfaFile = $phGrDir . "og.fa";
  if ( ! -e $phGrOGfaFile || $runAll )
  {
    print "--- Generating fa file of outgroup seq's of $phGr       ";
    $cmd = "select_seqs.pl $quietStr -s $phGrOGseqIDsFile -i $phGrBigFaFile -o $phGrOGfaFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }

  my @uqSpp = unique($phGrSppTbl{$phGr});
  my @spp = sort { @{$spIDsTbl{$b}} <=> @{$spIDsTbl{$a}} } @uqSpp;

  for my $spIdx ( 0..$#spp )
  {
    $forkMgr->start and next; # do the fork

    my $sp = $spp[$spIdx];
    my @ids = @{$spIDsTbl{$sp}}; # seq IDs of $sp
    my $nSp = scalar(@ids); # number of sequences of the given species

    print "\n--- Processing $sp (n=" . commify($nSp) . ")                                    \n";

    my $spDir = $phGrDir . $sp . "_dir";
    my $cmd = "mkdir -p $spDir";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    my $spReport = $spDir . "/report.txt";
    open my $ROUT, ">$spReport" or die "Cannot open $spReport for writing: $OS_ERROR";

    my $spSeqsFile = "$spDir/$sp" . ".seqIDs";
    writeArray(\@ids, $spSeqsFile);

    my $spFaFile = "$spDir/$sp" . ".fa";
    if ( ! -e $spFaFile || $runAll )
    {
      print "\r\t\tCreating sp fa file                  ";
      $cmd = "select_seqs.pl $quietStr -s $spSeqsFile -i $mmMasterFaFile -o $spFaFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    }

    my $spNRfaFile  = "$spDir/$sp" . "_nr.fa";
    my $spUCfile    = "$spDir/$sp" . ".uc";
    my $spUCfilelog = "$spDir/$sp" . "_uc.log";
    if ( ! -e $spNRfaFile || $runAll )
    {
      print "\r\t\tDereplicating species fasta file                   ";
      $cmd = "$usearch6 -cluster_fast $spFaFile -id 1.0 -uc $spUCfile -centroids $spNRfaFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;
    }

    my $nrSeqIDsFile = "$spDir/$sp" . "_nr.seqIDs";
    if ( ! -e $nrSeqIDsFile || $runAll )
    {
      print "\r\t\tExtracting non-redundant seq IDs               ";
      ## extracting seq IDs from the alignment file and selecting those IDs from the taxon file
      $cmd = "extract_seq_IDs.pl -i $spNRfaFile -o $nrSeqIDsFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    my @nrSeqIDs = readArray($nrSeqIDsFile);
    my $nnrSp = @nrSeqIDs;
    print "\nNo. nr seq IDs: " . commify($nnrSp) . "\n";

    print $ROUT "------------------------------------------------\n\n";
    print $ROUT "$sp\n\n";
    print $ROUT "n:     " . commify($nSp) . "\n";
    print $ROUT "n(nr): " . commify($nnrSp) . "\n";

    my $spClstr2File = "$spDir/$sp" . "_nr.clstr2";
    if ( ! -e $spClstr2File || $runAll )
    {
      print "\r\t\tCreating clstr2 file                               ";
      $cmd = "uc2clstr2.pl -i $spUCfile -o $spClstr2File";
      print "cmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    print "\r\t\tParsing clstr2 file                   ";
    my %cTbl = parseClstr2($spClstr2File);

    # sort cluster reference sequence IDs w/r cluster size
    @nrSeqIDs = sort { $cTbl{$b} <=> $cTbl{$a} } keys %cTbl;

    # print "\n\nFirst 10 ref seq's and the corresponding cluster sizes\n";
    # map { print "\t$_\t$cTbl{$_}\n" } @nrSeqIDs[0..10];

    my @clSizes = map { $cTbl{$_} } @nrSeqIDs;
    my $total   = sum (@clSizes); # total number of sequences
    my @clPercs = map { 100.0*$_/$total } @clSizes;

    my %clSizeTbl = map { $_ => $cTbl{$_} } @nrSeqIDs;
    my %clPercsTbl = map{ $_ => 100.0 * $clSizeTbl{$_} / $total } @nrSeqIDs;

    ## for validation of the above calculations and getting a sense what the thld for cumulative percentage to choose
    # my $outDir = "/Users/pgajer/projects/M_and_M/new_16S_classification_data/mm_speciation_validation_data/";
    # writeArray(\@clSizes, $outDir . "Li_cltr_sizes.txt");
    # writeArray(\@clPercs, $outDir . "Li_cltr_percs.txt");

    my $cumPerc = 0;
    my $percCovIdx = 0; # index of the sequence in @nrSeqIDs so that sum(clPercs[0..percCovIdx]) gives percCoverage
    for my $i (0..$#clPercs)
    {
      $cumPerc += $clPercs[$i];
      if ( $cumPerc > $percCoverage )
      {
	$percCovIdx = $i-1;
	last;
      }
    }

    $percCovIdx = 0 if $percCovIdx < 0;
    ##print "\npercCovIdx: $percCovIdx\n";

    my $covSuffix = "";
    if ( @nrSeqIDs > $maxNumNRseqs )
    {
      $covSuffix = "_nr_cov" . $percCoverage;

      my @nrSeqIDsMax = @nrSeqIDs[0..($maxNumNRseqs-1)];
      my @nrSeqIDsPer = @nrSeqIDs[0..$percCovIdx];

      ## If the number of seq's constituting $percCoverage is less than
      ## $maxNumNRseqs, we go with the later.
      if ( @nrSeqIDsPer > @nrSeqIDsMax )
      {
	@nrSeqIDs = @nrSeqIDsPer;
	print $ROUT "n(nr seq's covering $percCoverage" . "% of seq's classified to $sp): " . commify(scalar(@nrSeqIDs)) . "\n";
	print "\nNo. nr seq IDs covering $percCoverage" . "% of seq's classified to $sp: " . commify(scalar(@nrSeqIDs)) . "\n";
      }
      else
      {
	my $percCov = 0;
	for my $i (0..($maxNumNRseqs-1))
	{
	  $percCov += $clPercs[$i];
	}
	#my $perc = 100 *

	@nrSeqIDs = @nrSeqIDsMax;

	print $ROUT "$maxNumNRseqs no. of nr seq's covering $percCov" . "% of seq's classified to $sp): " . commify(scalar(@nrSeqIDs)) . "\n";
	print "\n$maxNumNRseqs no. of nr seq's covering $percCov" . "% of seq's classified to $sp: " . commify(scalar(@nrSeqIDs)) . "\n";
      }

      $nrSeqIDsFile = "$spDir/$sp" . $covSuffix . ".seqIDs";
      writeArray(\@nrSeqIDs, $nrSeqIDsFile);

      ## Restricting nr fa file to only nr ref seq's covering $percCoverage of all seq's
      my $spNRfaFile2 = "$spDir/$sp" . $covSuffix . ".fa";
      print "\r\t\tCreating restricted $sp fa file                  ";
      $cmd = "select_seqs.pl $quietStr -s $nrSeqIDsFile -i $spNRfaFile -o $spNRfaFile2";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

      $spNRfaFile = $spNRfaFile2;

    } # end of      if ( @nrSeqIDs > $maxNumNRseqs )

    ##
    ## 2. Generate alignment
    ##
    my $bigAlgnFile = "$spDir/$sp" .  $covSuffix . "_algn.fa";
    if ( ! -e $bigAlgnFile || $runAll )
    {
      print "\r\t\tAligning phGr ref seq's (includeing OG seq's) and the selected seq's of $sp           ";

      my @tmp;
      push (@tmp,"align.seqs(candidate=$spNRfaFile, template=$phGrAlgnFile, processors=8, flip=T)");

      printArray(\@tmp, "mothur commands") if ($debug || $verbose);

      my $scriptFile = create_mothur_script(\@tmp);
      $cmd = "$mothur < $scriptFile; rm -f $scriptFile mothur.*.logfile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?" if !$dryRun;


      my $mothurAlgnFile = $spNRfaFile; # "$trDir/" . $candBasename . ".align";
      $mothurAlgnFile =~ s/fa$/align/;
      print "mothurAlgnFile: $mothurAlgnFile\n" if $debug;

      $cmd = "rm -f $bigAlgnFile; cat $mothurAlgnFile $phGrAlgnFile > $bigAlgnFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    ##
    ## 3. Generate phylo tree
    ##
    my $bigNotRootedTreeFile = "$spDir/$sp" . $covSuffix . "_not_rooted_with_OGs.tree";
    if ( ! -e $bigNotRootedTreeFile || $runAll )
    {
      print "\r\t\tGenerating phylo tree of the above alignment                                    ";
      $cmd = "rm -f $bigNotRootedTreeFile; FastTree -nt $bigAlgnFile > $bigNotRootedTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    ## Rerooting the tree
    my $bigTreeWithOGsFile = "$spDir/$sp" . $covSuffix . "_with_OGs.tree";
    if ( ! -e $bigTreeWithOGsFile || $runAll )
    {
      print "\r\t\tRerooting the tree using outgroup sequences                          ";
      $cmd = "rm -f $bigTreeWithOGsFile; nw_reroot $bigNotRootedTreeFile @ogSeqIDs | nw_order -  > $bigTreeWithOGsFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    ## Pruning tree froom OG seq's
    my $bigTreeFile = "$spDir/$sp" . $covSuffix . ".tree";
    if ( ! -e $bigTreeFile || $runAll )
    {
      print "\r\t\tPruning the tree from OG seq's                                          ";
      $cmd = "rm -f $bigTreeFile; nw_prune $bigTreeWithOGsFile @ogSeqIDs | nw_order -  > $bigTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    ##
    ## 4. Running vicut on the tree using seq's of $sp as query nodes
    ##
    my $vicutDir    = "$spDir/$sp" . $covSuffix . "_vicut_dir";
    my $annFile     = $phGrTxFile;
    my $queryFile   = $nrSeqIDsFile;
    ##my $vicutTxFile = $vicutDir . "/minNodeCut_NAge1_TXge1_querySeqs.taxonomy"; # minNodeCut.cltrs
    my $vicutCltrsFile = $vicutDir . "/minNodeCut.cltrs";
    if ( ! -e $vicutCltrsFile || $runAll )
    {
      print "\r\t\tRunning vicut                                                              ";
      $cmd = "vicut $quietStr -t $bigTreeFile -a $annFile -q $queryFile -o $vicutDir";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    }

    ##
    ## 5. Reporting results of vicut
    ##

    my ($rvCltrTbl, $rvTxTbl, $rvExtTxTbl) = readCltrsTbl($vicutCltrsFile);

    my %vCltrTbl   = %{$rvCltrTbl};  # seqID => vicut cluster ID
    my %vTxTbl     = %{$rvTxTbl};    # seqID => taxonomy (NA for query seq's)
    my %vExtTxTbl  = %{$rvExtTxTbl}; # seqID => taxonomy of seqID if seqID is a phGr ref seq and c<vicut cluster ID of seqID> if seqID is a query seq

    my $vExtTxTblFile = "$spDir/$sp" . $covSuffix . "_ext.tx";
    writeTbl(\%vExtTxTbl, $vExtTxTblFile);

    ## vicut-cltr/tx frequency table
    my %vCltrvTxFreq;
    my %vCltrvTxIds;  # $vCltrvTxIds{cltr}{tx} = ref to seqID of the cluster's, cltr, taxon, tx.
    my %vCltrIds;
    for my $id ( keys %vCltrTbl )
    {
      $vCltrvTxFreq{$vCltrTbl{$id}}{$vTxTbl{$id}}++;
      push @{$vCltrvTxIds{$vCltrTbl{$id}}{$vTxTbl{$id}}}, $id;
      push @{$vCltrIds{$vCltrTbl{$id}}}, $id;
    }

    ## Identifing clusters that contain query sequences
    my @nrSeqCltrs; # = @vCltrTbl{@nrSeqIDs};
    ##print "\n\nvCltrTbl\n";
    for (@nrSeqIDs)
    {
      if ( exists $vCltrTbl{$_} )
      {
	push @nrSeqCltrs, $vCltrTbl{$_};
      }
      else
      {
	print "\n\nWARNING: $_ undefined in vCltrTbl\n";
      }
    }
    my @nrCltrs = unique(\@nrSeqCltrs);

    print "\nnrCltrs: @nrCltrs\n";

    ## size of each cluster
    my %vicutCltrSize;
    for my $cl ( @nrCltrs )
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
	warn "\nWARNING $cl not found in vCltrvTxFreq";
	print "\n";
      }
    }

    #print "\nFrequency table of vicut taxonomic assignments on selected nr seq's of $sp\n";
    my @nrSortedCltrs = sort { $vicutCltrSize{$b} <=> $vicutCltrSize{$a} } @nrCltrs;
    for my $cl (@nrSortedCltrs)
    {
      print "\nCluster $cl (" . $vicutCltrSize{$cl} . ")\n";
      print $ROUT "\n\nCluster $cl (" . $vicutCltrSize{$cl} . ")\n\n";

      ## Generating a list of species present in $cl sorted by size and with NA
      ## at the end (igoring the size of NA when sorting
      my @txs = keys %{$vCltrvTxFreq{$cl}};
      @txs = sort { $vCltrvTxFreq{$cl}{$b} <=> $vCltrvTxFreq{$cl}{$a} } @txs;
      ## putting NA at the end
      my @na = ("NA");
      @txs = diff(\@txs, \@na);
      push @txs, "NA";
      my %txSizes;
      for my $tx ( @txs )
      {
	$txSizes{$tx} = $vCltrvTxFreq{$cl}{$tx};
	#print "\t$tx\t" . $vCltrvTxFreq{$cl}{$tx} . "\n";
      }
      #print "\n";

      printFormatedTbl(\%txSizes, \@txs);
      printFormatedTblToFile(\%txSizes, \@txs, $ROUT);

      ## Reporting some characteristics of query seq's

      ## Coverage: percentage of seq's of the 100% identity clusters of NAs
      ## within all NAs' clusters
      my @clNRids = @{$vCltrvTxIds{$cl}{"NA"}};
      my $nCov = sum( @cTbl{ @clNRids } );
      my $pCov = sprintf( "%.1f%%", 100.0 * $nCov/ $nSp );
      print "Coverage: $pCov (" . commify($nCov) . " out of " . commify($nSp) . " seq's)\n";
      print $ROUT "Coverage: $pCov (" . commify($nCov) . " out of " . commify($nSp) . " seq's)\n";

      ## Size ranks
      my %clNRidsTbl = map { $_ => 1 } @clNRids;
      my @sizeRanks = grep { exists $clNRidsTbl{$nrSeqIDs[$_ - 1]}  } 1..($#nrSeqIDs+1);
      my @sizeRanks0 = grep { exists $clNRidsTbl{$nrSeqIDs[$_]}  } 0..$#nrSeqIDs;
      print "Size ranks: @sizeRanks\n";
      print $ROUT "Size ranks: @sizeRanks\n";

      ## Size percentage
      my @clSizePercs = @clPercsTbl{ @nrSeqIDs[ @sizeRanks0 ] };
      @clSizePercs = map { sprintf("%.2f", $_) } @clSizePercs;
      print "Size %'s: @clSizePercs\n";
      print $ROUT "Size %'s: @clSizePercs\n";

      ## pp's
      my @pps = @ppTbl{ @nrSeqIDs };
      @pps = @pps[ @sizeRanks0 ];
      print "pp's: @pps\n";
      print $ROUT "pp's: @pps\n";

      ## range(pp)
      my $minpp = min @pps;
      my $maxpp = max @pps;
      print "range(pp): [$minpp, $maxpp]\n";
      print $ROUT "range(pp): [$minpp, $maxpp]\n";
    }
    print "\n";
    print $ROUT "\n";

    ##
    ## 6. Generating and maybe viewing the tree
    ##
    ## Collapsing the tree using ref tx and vicut tx on query seq's

    print "\r\t\tGenerating a condensed tree of ref seq's species and vicut tx clades collapsed to a single node  ";
    my $condTreeFile2 = "$spDir/$sp" . $covSuffix . "_spp_cond2.tree";
    if ( ! -e $condTreeFile2 || $runAll )
    {
      print "\r\t\tGenerating a tree with species names at leaves ";
      my $sppTreeFile = "$spDir/$sp" . $covSuffix . "_spp.tree";
      $cmd = "rm -f $sppTreeFile; nw_rename $bigTreeFile $vExtTxTblFile | nw_order -c n  - > $sppTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

      my $condTreeFile = "$spDir/$sp" . $covSuffix . "_spp_cond1.tree";
      $cmd = "rm -f $condTreeFile; nw_condense $sppTreeFile > $condTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

      ## Relabeling tree so that only $sp and vicut tx taxonomy is left

      my $condSppLeavesFile = "$spDir/$sp" . $covSuffix . "_cond_spp.leaves";
      $cmd = "rm -f $condSppLeavesFile; nw_labels -I $condTreeFile > $condSppLeavesFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

      my @condSppTreeLeaves = readArray($condSppLeavesFile);

      ##print "condSppTreeLeaves: @condSppTreeLeaves\n";

      my %spThatMingleWithQuery;
      for my $cl (@nrCltrs)
      {
	my @txs = keys %{$vCltrvTxFreq{$cl}};
	for my $tx (@txs)
	{
	  if ( $tx ne "NA" )
	  {
	    $spThatMingleWithQuery{$tx} = 1;
	  }
	}
      }

      ##my @matches = grep { /pattern/ } @condSppTreeLeaves;
      my %newLeafNames;
      for my $l (@condSppTreeLeaves)
      {
	if ( exists $spThatMingleWithQuery{$l} || $l =~ /^c\d+/ )
	{
	  $newLeafNames{$l} = $l;
	}
	else
	{
	  $newLeafNames{$l} = "*";
	}
      }

      my $condSppLeavesFile2 = "$spDir/$sp" . $covSuffix . "_spp_cond.leaves2";
      writeTbl(\%newLeafNames, $condSppLeavesFile2);

      $cmd = "rm -f $sppTreeFile; nw_rename $condTreeFile $condSppLeavesFile2 | nw_order -c n  - > $condTreeFile2";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    }

    ## Producing pdf file of the tree sending it to either dir of species where
    ## vicut taxonomy agrees with PECAN's or to dir with spp for which there is
    ## a disagreement.

    my $pdfTreeFile = $treesDir . "/$sp" . "_tree.pdf";
    if ( ! -e $pdfTreeFile || $runAll )
    {
      my $treeAbsPath = abs_path( $condTreeFile2 );
      plot_tree($treeAbsPath, $pdfTreeFile, $sp);
    }

    if ( $showTree && $OSNAME eq "darwin")
    {
      $cmd = "open $pdfTreeFile";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
    }

    close $ROUT;

    $forkMgr->finish;

  } ## end of    for my $spIdx (0..

  # $forkMgr->wait_all_children; # this is so that we do not advance to the next
  # 			 # phylo-group before all species of that group were
  # 			 # processed

} ## end of   for my $phGr ( keys %phGrSppTbl )


## report timing
$endRun = time();
$runTime = $endRun - $startRun;
if ( $runTime > 60 )
{
  $timeMin = int($runTime / 60);
  $timeSec = sprintf("%02d", $runTime % 60);
  print "\rCompleted in $timeMin:$timeSec\n"
}
else
{
  print "\rCompleted in $runTime seconds\n"
}

print "\n\n\tOutput written to $outDir\n\n";

####################################################################
##                               SUBS
####################################################################

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
    ##$tbl{$refId} = \@ids;
    $tbl{$refId} = @ids; # we are only interested in the size of the cluseter
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

# read two column table; create a table that assigns
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
sub readArray{

  my ($file, $hasHeader) = @_;
  my @rows;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readArray(): $file does not exist";
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
    print "\t$_$pad" . $rTbl->{$_} . "\n";
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

  runRscript( $Rscript );
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

  runRscript( $Rscript );
}


# execute an R-script
sub runRscript
{
  my $Rscript = shift;

  my ($fh, $inFile) = tempfile("rTmpXXXX", SUFFIX => '.R', OPEN => 1, DIR => $tmpDir);
  print $fh "$Rscript";
  close $fh;

  my $outFile = $inFile . "out";
  my $cmd = "R CMD BATCH $inFile $outFile";
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

  if ( ! -f $file )
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

sub create_mothur_script{

    my (@arr) = @{$_[0]};
    my $file = "mothur_script.txt"; ##tmpnam();
    open OUT, ">$file" or die "Cannot open file $file to write: $!\n";
    foreach my $c (@arr){
        print OUT $c . "\n";
    }
    print OUT "quit()\n";

    return $file;
}

# print array to stdout
sub printArray
{
  my ($a, $header) = @_;
  print "\n$header\n" if $header;
  map {print "$_\n"} @{$a};
}

exit 0;
