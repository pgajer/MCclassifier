#! /usr/bin/perl

=head1 NAME

  clstr_to_clstr2.pl

=head1 DESCRIPTION

  transfer cd-hit clstr file into clstr2 file
  clstr2 file format
  seedId,s1,s2,...

=head1 SYNOPSIS

  clstr_to_clstr2.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS


=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /Users/pgajer/Shared_VM/33w_preg_poolFiles

  clstr_to_clstr2.pl -i AV730_w50/results.clstr -o AV730_w50/results.clstr2

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
  "output-file|o=s" => \my $outFile,
  "echo"            => \my $echo,
  "dry-run"         => \my $dryRun,
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
elsif (!$outFile)
{
  print "ERROR: Missing output file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
}

####################################################################
##                               MAIN
####################################################################

my %cltr = readCdHitCltrFile2("$inFile");

open OUT, ">$outFile" or die "Cannot open $outFile for reading: $OS_ERROR\n";
while (my ($k,$v) = each %cltr)
{
  print OUT "$k";
  foreach (@{$v})
  {
    print OUT ",$_";
  }
  print OUT "\n";
}
close OUT;
print "Output written to $outFile\n" if $echo;


####################################################################
##                               SUBS
####################################################################

# parse cd-hit-est sequence clustering file
sub readCdHitCltrFile2{

  my $cltrFile = shift;

  my %cltrTbl; # $cltrTbl{$id} = ref to array of reads associated with $id


  $INPUT_RECORD_SEPARATOR=">Cluster";
  open IN, "$cltrFile" or die "Cannot open $cltrFile for reading: $OS_ERROR\n";
  my @in = <IN>;
  shift @in;			# got rid of first '>'

  ##my $count = 0;

  foreach my $rec (@in)
  {
    my @lines = split /\n/, $rec;
    shift @lines;		# got rid of >Cluster line
    my $refId;
    my @reads;

    foreach my $line (@lines)
    {
      my ($readId) = ($line =~ />(\S+)/);

      next if $readId =~ /Cluster/;

      $readId =~ s/\.\.\.$//;

      if ( $line =~ /\*$/ )
      {
	$refId = $readId;
      }

      push @reads, $readId;
    }

    ##print "ref=$refId\treads=@reads\n";
    $cltrTbl{$refId} = \@reads;

    # last if $count > 3;
    # $count++;
  }

  close IN;
  $INPUT_RECORD_SEPARATOR="\n";

  return %cltrTbl;
}

exit;
