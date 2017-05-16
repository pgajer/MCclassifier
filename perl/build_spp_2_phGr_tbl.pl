#!/usr/bin/env perl

=head1 NAME

  build_spp_2_phGr_tbl.pl

=head1 DESCRIPTION

  Generating <species> <size> <phylo-group> table

  Data from 43 phyl-groups is being used here to build master models.

=head1 SYNOPSIS

  build_spp_2_phGr_tbl.pl -o <output file> [Options]

=head1 OPTIONS

=over

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

  build_spp_2_phGr_tbl.pl -o spp_2_phGr_ext_may12_models.txt

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

GetOptions(
  "out-file|o=s"        => \my $outFile,
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

if ( !$outFile )
{
  warn "\n\n\tERROR: Missing output file";
  print "\n\n";
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

## Gather final lineage, fasta (ungapped) and taxon files for all phylo-groups

my $baseDirExt = "/Users/pgajer/devel/MCextras/data/vag_exp_V3V4_phGrps_May16_dir/";

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

my @spLiFilesExt1 = map{ $baseDirExt . $_ } @spLiFilesExt0;


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


my %spLiExtTbl;
for my $i ( 0..$#extPhGrs )
{
  $spLiExtTbl{$extPhGrs[$i]} = $spLiFilesExt1[$i];
}


my $baseDir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";

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

my @spLiFiles1 = map{ $baseDir . $_ } @spLiFiles0;
## print "spLineageFile: @spLiFiles\n";


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


# get phylo dir of each group
my @phDir;
my %phDirTbl;
for my $i ( 0..$#allPhGrs )
{
  my @f = split /\//, $spLiFiles0[$i];
  #print "f: @f\n"; exit;
  my $dir = shift @f;
  push @phDir, $dir;
  $phDirTbl{ $allPhGrs[$i] } = $dir;
}

# print_array(\@phDir, "phDir");
# print "\n";
# exit;

# list species in allPhGrs but not in extPhGrs

my @d = diff( \@allPhGrs, \@extPhGrs );
print "\nNumber of non-extended phylo-groups: " . @d . "\n";
print_array(\@d);
print "\n\n";

if ( 0 )
{
  my $dir = "/Users/pgajer/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir/";
  my $count = 1;
  for my $phGr ( @d )
  {
    print "--- [$count] Processing $phGr\n";
    $count++;
    my $phDir = $phDirTbl{$phGr};
    my $cmd = "scp $dir$phDir/$phGr"  ."_dir.tgz cadbane.igs.umaryland.edu:/usr/local/projects/pgajer/devel/MCextras/data/RDP/V3V4";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?";# if !$dryRun;
  }
  exit;
}


my %spLiTbl;
for my $i ( 0..$#allPhGrs )
{
  $spLiTbl{$allPhGrs[$i]} = $spLiFiles1[$i];
}

my @spLiFiles;
my %spLiFile2phGr;
for ( @allPhGrs )
{
  if ( exists $spLiExtTbl{$_} )
  {
    push @spLiFiles, $spLiExtTbl{$_};
    $spLiFile2phGr{ $spLiExtTbl{$_} } = $_;
  }
  else
  {
    push @spLiFiles, $spLiTbl{$_};
    $spLiFile2phGr{ $spLiTbl{$_} } = $_;
  }
}


if ( $debug )
{
  print_array( \@spLiFiles, "\nspLiFiles" );
  print "\n\n";

  exit;
}


open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
for my $file ( @spLiFiles )
{
  my $phGr = $spLiFile2phGr{$file};
  my @spp = get_spp_from_spLi( $file );
  for ( @spp )
  {
    print OUT "$_\t$phGr\n";
  }
}
close OUT;


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

# extracting species (first column) from the spLineage file

# $ head Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.spLineage
# Ilumatobacter_sp_3	sg_Aciditerrimonas_Euzebya_Ilumatobacter_etal	g_Gaiella_Rubrobacter_Ferrimicrobium_etal	f_Acidimicrobiaceae_Gaiellaceae_Rubrobacteraceae	o_Acidimicrobiales_Gaiellales_Rubrobacterales	c_Actinobacteria	p_Actinobacteria	d_Bacteria

sub get_spp_from_spLi
{
  my $file = shift;

  my @spp;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    my @f = split /\s+/;
    push @spp, shift @f;
  }
  close IN;

  return @spp;
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

exit 0;
