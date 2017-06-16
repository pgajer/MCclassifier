#! /usr/bin/perl

=head1 NAME

  chg_fasta_seqIDs.pl

=head1 DESCRIPTION

  Change sequence IDs in a given fasta file to names supplied by an translation
  table. Only seqIDs present in the table are changed.

=head1 SYNOPSIS

  chg_fasta_seqIDs.pl -i <fasta file> -a <translation table file> -o <output file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input fasta file.

=item B<--tr-file, -a>
  Translation table file

=item B<--output-file, -o>
  Output file.

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  ~/projects/PECAN/data/Banfield_contax/FL

  mv Banfield_medoids_FL_algn.fa Banfield_medoids_FL_algn0.fa
  chg_fasta_seqIDs.pl -i Banfield_medoids_FL_algn0.fa -a ../Hug_Banfield_old_to_new.seqIDs -o Banfield_medoids_FL_algn.fa

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
  "input-file|i=s"  => \my $inFile,
  "tr-file|a=s"     => \my $trFile,
  "output-file|o=s" => \my $outFile,
  "help|h!"         => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

if (!$inFile)
{
  print "ERROR: Missing input file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$trFile)
{
  print "ERROR: Missing taxon file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}
elsif (!$outFile)
{
  print "ERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################

## parsing annotation file
my %trTbl = read_tbl($trFile);  # $trTbl{seqId} = taxonomic assignment of seqId

open (FASTA, "<$inFile") or die "Cannot open $inFile for reading: $OS_ERROR\n";
$/ = ">";
my $junkFirstOne = <FASTA>;

open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
my $count = 1;
while (<FASTA>)
{
  if ($count % 100 == 0)
  {
    print "\r$count";
  }

  chomp;
  my ($def,@seqlines) = split /\n/, $_;
  my ($id) = split /\s+/, $def;
  my $seq = join '', @seqlines;

  if ( exists $trTbl{$id} )
  {
    print OUT ">" . $trTbl{$id} . "\n$seq\n";
  }
  else
  {
    print OUT ">$id\n$seq\n";
  }
  $count++;
}
close OUT;

print "\rOutput written to $outFile\n";

####################################################################
##                               SUBS
####################################################################

# read two column table; create a table that assigns
# elements of the first column to the second column
sub read_tbl{

  my $file = shift;

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

exit;
