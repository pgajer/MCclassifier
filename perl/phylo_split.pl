#!/usr/bin/env perl

=head1 NAME

  phylo_split.pl

=head1 DESCRIPTION

  The aim of the script is to generate a split of all bacterial (and Archeal)
  sequences into phylogenetically sound sub-groups.

=head1 SYNOPSIS

  phylo_split.pl

=head1 OPTIONS

=over

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

  phylo_split.pl

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
  "offset-coef|o=f"     => \$offsetCoef,
  "tx-size-thld|t=i"    => \$txSizeThld,
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

my @liFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4_final.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4_final.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4_final.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4_final.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4_final.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4_final.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4_final.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4_final.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4_final.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4_final.lineage",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4_final.lineage",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4_final.lineage",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4_final.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4_final.lineage");

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


if ( $debug )
{
  print_array( \@liFiles, "\nliFiles" );
  print "\n\n";

  exit;
}


##
## concatenating all lineage files
##

print "\r--- Concatenating all lineage files                              \n";

my %spTbl;
my %sgeTbl;
my %geTbl;
my %faTbl;
my %orTbl;
my %clTbl;
my %phTbl;

##my %lineage; # species => lineage of the species (recorded as a string) with the corresponding phyo-group name at the end
my %badSppName;
foreach $file (@liFiles)
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
    #my $sge = shift @f;
    my $ge = shift @f;
    my $fa = shift @f;
    my $or = shift @f;
    my $cl = shift @f;
    my $ph = shift @f;

    #
    if ( $sp =~ /MG/ || $sp =~ /FN/ || $sp =~ /\./)
    {
      warn "\n\nERROR: $sp has strange name";
      print "$file\n\n";
      $badSppName{$sp} = $file;
    }

    # are there any backslashes in the species name
    @f = split "/", $sp;
    if ( @f > 1 )
    {
      warn "\n\nERROR: $sp seems to have a backslash in its name";
      print "$file\n\n";
      $badSppName{$sp} = $file;
    }

     print "lineage: $lineage\n";
     print "sp : $sp\n";
    # print "sge : $sge\n";
     print "ge : $ge\n";
     print "fa : $fa\n";
     print "or : $or\n";
     print "cl : $cl\n";
     print "ph : $ph\n";
     exit;

    $spTbl{$sp}{$phGr}   = 1;
    $geTbl{$ge}{$phGr}   = 1;
    $faTbl{$fa}{$phGr}   = 1;
    $orTbl{$or}{$phGr}   = 1;
    $clTbl{$cl}{$phGr}   = 1;
    $phTbl{$ph}{$phGr}   = 1;
  }
  close IN;
}


if ( keys %badSppName > 0 )
{
  print "\n\nDiscovered the following species with suspecious names\n";
  my @a = sort { $badSppName{$a} cmp $badSppName{$b} } keys %badSppName;
  print_tbl( \%badSppName, \@a );

  print "\n\nPlease fix these names and rerun taxonomy_cleanup.pl on the corresonding phylo-groups\n\n";

  exit;
}


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

# print elements of a hash table
sub print_tbl
{
  my ($rTbl, $r) = @_;

  map {print "$_\t" . $rTbl->{$_} . "\n"} @{$r};
}

exit 0;
