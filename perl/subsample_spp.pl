#!/usr/bin/env perl

=head1 NAME

  subsample_spp.pl


=head1 DESCRIPTION

  Given subsample lineage file so that none of species has more than maxSp sequences.


=head1 SYNOPSIS

  subsample_spp.pl --max-sp-size <maxSp> -l <lineage file> -o <output lineage dir> [Options]

=head1 OPTIONS

=over

=item B<--lineage-file, -l>
  Lineage file.

=item B<--max-sp-size, -s>
  Maximal allowed number of sequences within a species.

=item B<--output-file, -o>
  Output lineage file.

=item B<--verbatim, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--quiet>
  Do not print progress messages.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/phylo_groups/v0.2

  subsample_spp.pl -m 100 -l master_V3V4_no_outliers.lineage -o master_V3V4_no_outliers_sp100.lineage

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
  "lineage-file|l=s"   => \my $liFile,
  "max-sp-size|m=i"    => \my $maxSpSize,
  "output-file|o=s"    => \my $outFile,
  "verbose|v"          => \my $verbose,
  "quiet"              => \my $quiet,
  "debug"              => \my $debug,
  "dry-run"            => \my $dryRun,
  "help|h!"            => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$liFile)
{
  print "\n\nERROR: Missing lineage file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$maxSpSize)
{
  print "\n\nERROR: Missing sample size\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif (!$outFile)
{
  print "\n\nERROR: Missing output file\n\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( ! -e $liFile || ! -s $liFile )
{
  print "\n\nERROR: $liFile does not exist\n\n\n";
  exit 1;
}


####################################################################
##                               MAIN
####################################################################

print "--- Parsing lineage table\n";
my ($rliTbl, $rspTbl) = parse_li_tbl( $liFile );

my %spTbl = %{$rspTbl}; # spTbl{$sp} ref to seq's of that species
my %liTbl = %{$rliTbl}; # liTbl{seqID} = lineage string; For example, S001097805 => Root;Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Streptomycetaceae;Streptomyces;Streptomyces_sp

print "--- Computing species frequencies\n";
my %spSize;
for my $sp ( keys %spTbl )
{
  $spSize{$sp} = @{$spTbl{$sp}};
}

print "--- Generating reduced size species table\n";
my %spTblR;
my %spSizeR;
for my $sp ( keys %spTbl )
{
  if ( $spSize{$sp} <= $maxSpSize )
  {
    $spTblR{$sp} = $spTbl{$sp};
    $spSizeR{$sp} = $spSize{$sp};
  }
  else
  {
    $spTblR{$sp} = sample( $spTbl{$sp}, $maxSpSize );
    $spSizeR{$sp} = @{$spTblR{$sp}};
  }
}

my @ids;
for my $sp ( keys %spTblR )
{
  push @ids, @{$spTblR{$sp}};
}

write_tbl2( \%liTbl, \@ids, $outFile );

my $n = commify( scalar(keys %liTbl) );
print "\nNumber of elements in the input lineage file:       $n\n";
print   "Number of elements in the restricted  lineage file: " . commify( scalar(@ids) ) . "\n";

my $spSizeFile = $liFile;
$spSizeFile =~ s/lineage/spSize/;
write_tbl( \%spSize, $spSizeFile );
print "\nSpecies sizes written to $spSizeFile\n";

my $spSizeRFile = $liFile;
$spSizeRFile =~ s/lineage/spSizeR/;
write_tbl( \%spSizeR, $spSizeRFile );
print "Restricted species sizes written to $spSizeRFile\n";
print "New lineage file written to $outFile\n\n";

####################################################################
##                               SUBS
####################################################################

sub parse_li_tbl
{
  my $file = shift;

  my %liTbl;
  my %spTbl;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  for ( <IN> )
  {
    chomp;
    my ($id, $li) = split /\s+/;
    $liTbl{$id} = $li;

    my @f = split ";", $li;
    my $sp = pop @f;
    push @{$spTbl{$sp}}, $id;
  }
  close IN;

  return (\%liTbl, \%spTbl);
}


# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle
{
    my $array = shift;

    my $i;
    for ($i = @$array; --$i; )
    {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

sub sample
{
  my ( $r, $max ) = @_;

  fisher_yates_shuffle( $r );
  my @ids = @{$r};
  $#ids = $max - 1;

  return \@ids;
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

# write hash table over some slice of its keys
sub write_tbl2
{
  my ($rTbl, $r, $outFile) = @_;

  open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR";
  map {print OUT $_ . "\t" . $rTbl->{$_} . "\n"} sort @{$r};
  close OUT;
}

## put commas in numbers for better readability
## lifted from
## http://www.perlmonks.org/?node_id=2145
sub commify
{
   local $_  = shift;
   s{(?<!\d|\.)(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     $n;
    }eg;
   return $_;
}

exit 0;
