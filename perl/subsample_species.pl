#!/usr/bin/env perl

=head1 NAME

  subsample_species.pl


=head1 DESCRIPTION

  Given a master lineage file (fields separated by tabs, spaces, or semicolons),
  produce a subset of this lineage file based on the number of sequences in each
  species. If the number of sequences is greater than the number provided, choose 
  randomly that number of sequences of that species and place their lineage into 
  the subset lineage file. Species with less than 100 sequences will have all 
  sequences included in the subset lineage file.


=head1 SYNOPSIS

    subsample_species.pl -l <lineage file> -n <maximum number of sequences per species>

=head1 OPTIONS

=over

=item B<--lineage-file, -l>
  Source lineage file.

=item B<--max-seq, -n>
  Maximum number of sequences to represent each species.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/local/projects/pgajer/projects/PECAN/data/phylo_groups/v0.2

  subsample_species.pl -l master_V3V4_no_outliers.lineage -n 100

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Statistics::Descriptive;
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);


$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "lineage-file|l=s" 	  => \my $lineageFile,
  "max-seq|n=s"  		  => \my $n,
  "quiet"           	  => \my $quiet,
  "verbose"           	  => \my $verbose,
  "debug"            	  => \my $debug,
  "dry-run"          	  => \my $dryRun,
  "help|h!"          	  => \my $help,
  "igs"                   => \my $igs,
  "johanna"               => \my $johanna,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!lineage)
{print "\n\nERROR: Missing source lineage file!\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!maxSeq)
{print "\n\nMissing desired maximum number of sequences.\n\n\n";
  print "\n\nDefaulting to 100 sequences per species.\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  $n=100;
}

####################################################################
##                               MAIN
####################################################################

#open master lineage file
#Read master lineage file (tabs, spaces, and ;)
print "--- Parsing lineage table\n";
my %lineageTbl = read2colTbl($lineageFile);

my $lineage;
my @levels;
my @species;
#For each entry of the array @uniq_species
foreach my $seqID (keys %lineageTbl)
{
	#Make an array of the 8th field of the file (species)
	$lineage = $lineageTbl{$_};
	@levels = split(/[; ]/, $lineage);
	push @species, $levels[8];
}

#make an array of this 8th field that is sorted and uniq'd, @uniq_species
my @uniqSpecies = uniq(@species);

my @allmatches;
my @randMatches;
#For each unique species, grep the master lineage file
foreach my $spp(@uniqSpecies)
{
	@allmatches = grep { $lineageTbl{$_} eq $spp} keys %lineageTbl;
	# count the number of hits, length()
	# if the number of hits > 100
	if scalar @allmatches > $n 
	{
		# randomly choose 100 of the seqIDs (keys)
		for (1 .. $n)
		{
			# write the seqIDs to array
			push @randMatches, splice @allmatches, rand @allmatches, 1;
		}
	}
	else
	{
		# write all of the seqIDs to array
		splice @randMatches, 1, 0, @allmatches;
	} 
}

#use the array of seqIDs to get the subset of lineages

my $subMasterLineage = "subsampled_species.lineage";
open OUT, ">$subMasterLineage" or die "Cannot open $subMasterLineage for writing: $OS_ERROR\n";

foreach $seq(@randMatches)
{
	$lineage = $lineageTbl{$_};
	print OUT "$seq\t$lineage\n";

}
close OUT;


####################################################################
##                               SUBS
####################################################################
sub read2colTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\nERROR in read2colTbl(): $file does not exist\n\n\n";
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

exit 0;