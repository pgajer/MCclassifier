#!/usr/bin/env perl

=head1 NAME

  clean_master_lineage_v2.pl

=head1 DESCRIPTION

  This script modifies species names so that they are of the form Genus_speciesname
  where Genus is capitilized and speciesname is a lowercase

  '/', '-' and '_' characters (except the '_' between genus and speciesname) are changed to '|'

  _gp1/2/3 changes to |gp1/2/3

  Escherichia-Shigella_sp => Escherichia|Shigella_sp

  sp.suffix => sp

  Genus name and prefix of species name have to be the same

=head1 SYNOPSIS

  clean_master_lineage_v2.pl -i <lineage file> -o <output lineage file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input lineage file.

=item B<--output-file, -o>
  Output lineage file.

=item B<--verbose, -v>
  Prints content of some output files. Default value: 5000.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd ~/projects/PECAN/data/RDP/r2.2

  clean_master_lineage_v2.pl -i rdp_Bacteria_fp_seqlen_amb_filtered_wBVAB_no_incertae_sedis_nr.lineage -o rdp_Bacteria_curated.lineage

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

GetOptions(
  "input-file|i=s"      => \my $inFile,
  "output-file|o=s"      => \my $outFile,
  "quiet"               => \my $quiet,
  "igs"                 => \my $igs,
  "verbose|v"           => \my $verbose,
  "dry-run"             => \my $dryRun,
  "debug"               => \my $debug,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);


if ( $help )
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if ( !$inFile )
{
  print "ERROR: Missing input lineage file\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}
elsif ( !$outFile )
{
  print "ERROR: Missing output file name\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

my $debugStr = "";
my $quietStr = "--quiet";
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

####################################################################
##                               MAIN
####################################################################

print "--- Parsing lineage table\n";
my $suffix;
my %liTbl;
open IN, "$inFile" or die "Cannot open $inFile for reading: $OS_ERROR\n";
for ( <IN> )
{
  chomp;
  my ($id, $li) = split /\s+/;

  my @f = split ";", $li;
  my $sp = pop @f;
  my $ge = pop @f;
  my $fa = pop @f;

  my $origSp = $sp;
  #$sp = "AAA-BBB_gp1_sp.sldur-23X";
  #$sp = "Curtobacterium_MG-2011-84-GV";
  #$sp = "Mycobacterium_n";
  #$sp = "Granulicatella_para|adiacens";
  #$sp = "Selenomonas_flueggei|like";
  #print "sp BEFORE: $sp\n";
  $sp =~ s/\//\|/g;
  $sp =~ s/\-/\|/g;
  $sp =~ s/_gp\d//;
  $sp =~ s/\..+//;
  $sp =~ s/_[[:upper:]].+/_sp/;
  $sp =~ s/_\w$/_sp/;
  $sp =~ s/Ruminococcus2/Ruminococcus/;
  $sp =~ s/Armatimonas\|Armatimonadetes/Armatimonas/;
  $sp =~ s/Chthonomonas\|Armatimonadetes/Chthonomonas/;
  $sp =~ s/Escherichia\|Shigella/EscherichiaShigella/;
  $sp =~ s/_para\|/_/;
  $sp =~ s/\|like$//;
  $sp =~ s/lxb\|14/lxb14/;
  $sp =~ s/lxb\|3/lxb3/;
  #print "sp AFTER: $sp\n";  exit 1;

  if ( $sp !~ /[[:upper:]][\||\w]+_[[:lower:]]\w+/ && $sp !~ /BVAB/ )
  {
    warn "\n\n\tERROR: incorrect format sp: $sp";
    print "origSp: $origSp\n\n";
    exit 1;
  }

  if ( $sp !~ /BVAB/ )
  {
    ($ge, $suffix) = split "_", $sp;
  }

  if ( $ge eq "Thiomonas" || $ge eq "Ideonella" || $ge eq "Paucibacter" || $ge eq "Methylibium"
       || $ge eq "Mitsuaria" || $ge eq "Thiobacter" || $ge eq "Xylophilus" || $ge eq "Piscinibacter"
       || $ge eq "Aquincola" || $ge eq "Sphaerotilus" || $ge eq "Inhella" || $ge eq "Leptothrix"
       || $ge eq "Tepidimonas" || $ge eq "Rubrivivax" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Burkholderiales.family";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Cellvibrio" || $ge eq "Marinimicrobium"
          || $ge eq "Teredinibacter" || $ge eq "Umboniibacter"
          || $ge eq "Saccharophagus" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Cellvibrionales;Cellvibrionaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Caldalkalibacillus" || $ge eq "Calditerricola" || $ge eq "Viridibacillus"
          || $ge eq "Microaerobacter" || $ge eq "Lysinibacillus" )
  {
    $li = "Root;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Haliea" || $ge eq "Agarivorans"
          || $ge eq "Alishewanella" || $ge eq "Marinobacterium" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Alteromonadaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Butyricicoccus" || $ge eq "Alkaliphilus"
          || $ge eq "Clostridiisalibacter" || $ge eq "Saccharofermentans" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Alkaliflexus" || $ge eq "Anaerophaga" || $ge eq "Marinilabilia" )
  {
    $li = "Root;Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Marinilabiliaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Caldanaerobius" || $ge eq "Syntrophaceticus" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Thermoanaerobacterales;Thermoanaerobacteraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Ruminococcus" || $ge eq "Pseudoflavonifractor" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Nesiotobacter" || $ge eq "Labrenzia" )
  {
    $li = "Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhodobacterales;Rhodobacteraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Wohlfahrtiimonas" || $ge eq "Ignatzschineria" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Xanthomonadaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Methylocaldum" || $ge eq "Methylococcus" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Methylococcales;Methylococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Methylocaldum" || $ge eq "Methylococcus" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Methylococcales;Methylococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Nevskia" || $ge eq "Hydrocarboniphaga" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Nevskiales;Nevskiaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Vasilyevaea" || $ge eq "Amorphus" )
  {
    $li = "Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiales.family";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Salicola" || $ge eq "Spongiispira" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Oceanospirillales;Oceanospirillales.family";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Atopobacter" )
  {
    $li = "Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Carnobacteriaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Flavonifractor" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales.family";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Meganema" )
  {
    $li = "Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Methylobacteriaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Nitrosococcus" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Chromatiales;Chromatiaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Maritalea" )
  {
    $li = "Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Hyphomicrobiaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Pleomorphomonas" )
  {
    $li = "Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Methylocystaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Dethiobacter" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Syntrophomonadaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Bavariicoccus" )
  {
    $li = "Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Thiobacillus" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Hydrogenophilales;Hydrogenophilaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Aquabacterium" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Thioalkalispira" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Chromatiales;Thioalkalispiraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Howardella" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Akkermansia" )
  {
    $li = "Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Dasania" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Cellvibrionales;Spongiibacteraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Sharpea" )
  {
    $li = "Root;Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Caedibacter" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichales.family";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Oceanotoga" )
  {
    $li = "Root;Bacteria;Thermotogae;Thermotogae;Thermotogales;Thermotogaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Rummeliibacillus" )
  {
    $li = "Root;Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Catellicoccus" )
  {
    $li = "Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Blautia" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Steroidobacter" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Caldithrix" )
  {
    $li = "Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacterales.family";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Ferritrophicum" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Ferritrophicales;Ferritrophicaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Rugamonas" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Oxalobacteraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Roseiflexus" )
  {
    $li = "Root;Bacteria;Chloroflexi;Chloroflexia;Chloroflexales;Roseiflexaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Proteocatella" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Peptostreptococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Alkanibacter" )
  {
    $li = "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Lutispora" )
  {
    $li = "Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Sinosporangium" )
  {
    $li = "Root;Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Streptosporangineae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Desulfurobacteriaceae" )
  {
    $li = "Root;Bacteria;Aquificae;Aquificae;Desulfurobacteriales";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Solibacillus" )
  {
    $li = "Root;Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }
  elsif ( $ge eq "Shinella" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Rhodocyclales;Rhodocyclaceae";
    @f = split ";", $li;
    $fa = pop @f;
  }

  ## Family
  if ( $fa eq "Gallionellaceae" )
  {
    $li = "Root;Bacteria;Proteobacteria;Betaproteobacteria;Nitrosomonadales";
    @f = split ";", $li;
  }
  elsif ( $fa eq "Rhodothermaceae" )
  {
    $li = "Root;Bacteria;Bacteroidetes;Cytophagia;Cytophagales";
    @f = split ";", $li;
  }
  elsif ( $fa eq "Nitrospinaceae" )
  {
    $li = "Root;Bacteria;Nitrospinae;Nitrospinia;Nitrospinales";
    @f = split ";", $li;
  }
  elsif ( $fa eq "Cyclobacteriaceae" )
  {
    $li = "Root;Bacteria;Bacteroidetes;Cytophagia;Cytophagales";
    @f = split ";", $li;
  }
  elsif ( $fa eq "Desulfurobacteriaceae" )
  {
    $li = "Root;Bacteria;Aquificae;Aquificae;Desulfurobacteriales";
    @f = split ";", $li;
  }

  my $or = pop @f;
  my $cl = pop @f;
  my $ph = pop @f;


  ## Order
  if ( $or eq "Planctomycetales" )
  {
    $li = "Root;Bacteria;Planctomycetes;Planctomycea";
    @f = split ";", $li;
    $cl = pop @f;
    $ph = pop @f;
  }
  elsif ( $or eq "Thermoleophilales" )
  {
    $li = "Root;Bacteria;Actinobacteria;Thermoleophilia";
    @f = split ";", $li;
    $cl = pop @f;
    $ph = pop @f;
  }
  elsif ( $or eq "Chlorobiales" )
  {
    $li = "Root;Bacteria;Chlorobi;Chlorobia";
    @f = split ";", $li;
    $cl = pop @f;
    $ph = pop @f;
  }
  elsif ( $or eq "Chlamydiales" )
  {
    $li = "Root;Bacteria;Chlamydiae;Chlamydiia";
    @f = split ";", $li;
    $cl = pop @f;
    $ph = pop @f;
  }
  elsif ( $or eq "Spirochaetales" )
  {
    $li = "Root;Bacteria;Spirochaetes;Spirochaetia";
    @f = split ";", $li;
    $cl = pop @f;
    $ph = pop @f;
  }

  ## Class
  if ( $cl eq "Ignavibacteria" )
  {
    $li = "Root;Bacteria;Chlorobi";
    @f = split ";", $li;
    $ph = pop @f;
  }


  $li = join ";", (@f, $ph, $cl, $or, $fa, $ge, $sp);
  $liTbl{$id} = $li;
  #print "$id\t$li\n"; exit 1;
}
close IN;

print "--- Checking uniquness of parents\n";
if ( each_tx_has_unique_parent(\%liTbl) )
{
  exit 1;
}

write_tbl( \%liTbl, $outFile );


####################################################################
##                               SUBS
####################################################################

## check if each node of the lineage structure has only one parent
sub each_tx_has_unique_parent
{
  my $r = shift;
  my %liTbl = %{$r};

  my  %prt;
  for my $id ( keys %liTbl )
  {
    my $lineage = $liTbl{$id};
    my @f = split ";", $lineage;
    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    $ge = "g_$ge";
    $fa = "f_$fa";
    $or = "o_$or";
    $cl = "c_$cl";
    $ph = "p_$ph";

    $prt{$sp}{$ge}++;
    $prt{$ge}{$fa}++;
    $prt{$fa}{$or}++;
    $prt{$or}{$cl}++;
    $prt{$cl}{$ph}++;
  }

  my $ret = 0;
  for my $tx (keys %prt)
  {
    my $nPrts = keys %{$prt{$tx}};
    if ( $nPrts > 1 )
    {
      $ret = 1;
      warn "\n\n\tERROR: $tx has more than one parent";
      print "\n\t$tx parents\n";
      for ( keys %{$prt{$tx}} )
      {
        print "\t\t$_\t" . $prt{$tx}{$_} . "\n";
      }
      print "\n\n";

      #print "Attempting to fix the problem. Going with more abundant lineage\n";
      #my $tmpLiFile = "tmp.lineage";
      #write_tbl( \%liTbl, $tmpLiFile );
      #print "\tLineage file written to $tmpLiFile\n\n";
    }
  }

  return $ret;
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

exit 0;
