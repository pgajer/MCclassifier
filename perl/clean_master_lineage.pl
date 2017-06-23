#!/usr/bin/env perl

=head1 NAME

  clean_master_lineage.pl

=head1 DESCRIPTION

  This script reads lineage files from the full length phylo-groups and checks for

  1. consistency of species and genus names
  2. equality of family and order,
  3. order and class
  4. class and phylum ranks
  5. endings of families
  6. endings of orders

  It corrects inconsistencies and saves the curated species lineage file and V3V4 taxonomy.

=head1 SYNOPSIS

  clean_master_lineage.pl -o <output dir> [Options]

=head1 OPTIONS

=over

=item B<--output-dir, -o>
  Output dir with the following files.

  sp.lineage
  v3v4.tx

=item B<--verbose, -v>
  Prints content of some output files. Default value: 5000.

=item B<--debug>
  Prints system commands

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  clean_master_lineage.pl -o master_lineage_dir

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
  "output-dir|o=s"      => \my $outDir,
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

if ( !$outDir )
{
  print "ERROR: Missing output directory\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
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
my $readNewickFile        = "/Users/pgajer/organizer/programming/R/libs/read.newick.R";
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
  $ginsi                 = "/home/pgajer/bin/ginsi"; # MAFFT v7.310 (2017/Mar/17)
}

####################################################################
##                               MAIN
####################################################################

print "--- Get V3V4 lineage files\n";
my @liFiles = get_li_files();

if ( $debug )
{
  print_array( \@liFiles, "\nliFiles" );
  print "\n\n";
  exit;
}

my $baseDir  = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.1/";
#my $spGeFile = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/species_genus_tbl_may19_2017.txt";
my $geLiFile = "/Users/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this file has a header !

print "--- Building taxonomic parent tables\n";

#my $baseDir = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.1/";
my @mLiFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_dir/Actinobacteria_group_0.lineage",
		 "Actinobacteria_dir/Actinobacteria_group_1_dir/Actinobacteria_group_1.lineage",
		 "Actinobacteria_dir/Actinobacteria_group_2_dir/Actinobacteria_group_2.lineage",
		 "Actinobacteria_dir/Actinobacteria_group_3_dir/Actinobacteria_group_3.lineage",
		 "Actinobacteria_dir/Actinobacteria_group_4_dir/Actinobacteria_group_4.lineage",
		 "Actinobacteria_dir/Actinobacteria_group_5_dir/Actinobacteria_group_5.lineage",
		 "Bacteroidetes_dir/Bacteroidetes_group_0_dir/Bacteroidetes_group_0.lineage",
		 "Bacteroidetes_dir/Bacteroidetes_group_1_dir/Bacteroidetes_group_1.lineage",
		 "Bacteroidetes_dir/Bacteroidetes_group_2_dir/Bacteroidetes_group_2.lineage",
		 "Bacteroidetes_dir/Bacteroidetes_group_3_dir/Bacteroidetes_group_3.lineage",
		 "small_phyla/Chloroflexi_dir/Chloroflexi.lineage",
		 "small_phyla/Deinococcus-Thermus_dir/Deinococcus-Thermus.lineage",
		 "small_phyla/Fusobacteria_dir/Fusobacteria.lineage",
		 "small_phyla/Nitrospirae_dir/Nitrospirae.lineage",
		 "small_phyla/Planctomycetes_dir/Planctomycetes.lineage",
		 "small_phyla/Spirochaetes_dir/Spirochaetes.lineage",
		 "small_phyla/Tenericutes_dir/Tenericutes.lineage",
		 "small_phyla/Verrucomicrobia_dir/Verrucomicrobia.lineage",
		 "small_phyla/phyla_lessthen_1k_wOG_dir/phyla_lessthen_1k_wOG.lineage",
		 "Firmicutes_dir/Firmicutes_group_0_dir/Firmicutes_group_0.lineage",
		 "Firmicutes_dir/Firmicutes_group_1_dir/Firmicutes_group_1.lineage",
		 "Firmicutes_dir/Firmicutes_group_2_dir/Firmicutes_group_2.lineage",
		 "Firmicutes_dir/Firmicutes_group_3_dir/Firmicutes_group_3.lineage",
		 "Firmicutes_dir/Firmicutes_group_4_dir/Firmicutes_group_4.lineage",
		 "Firmicutes_dir/Firmicutes_group_5_dir/Firmicutes_group_5.lineage",
		 "Firmicutes_dir/Firmicutes_group_6_dir/Firmicutes_group_6.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_0_dir/Proteobacteria_group_0.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_10_dir/Proteobacteria_group_10.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_11_dir/Proteobacteria_group_11.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_12_dir/Proteobacteria_group_12.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_13_dir/Proteobacteria_group_13.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_14_dir/Proteobacteria_group_14.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_15_dir/Proteobacteria_group_15.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_17_dir/Proteobacteria_group_17.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_1_dir/Proteobacteria_group_1.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_2_dir/Proteobacteria_group_2.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_3_dir/Proteobacteria_group_3.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_4_dir/Proteobacteria_group_4.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_5_dir/Proteobacteria_group_5.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_6_dir/Proteobacteria_group_6.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_7_dir/Proteobacteria_group_7.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_8_dir/Proteobacteria_group_8.lineage",
		 "Proteobacteria_dir/Proteobacteria_group_9_dir/Proteobacteria_group_9.lineage");
my @mLiFiles1 = map{ $_ = $baseDir . $_ } @mLiFiles0;


print "--- Extracting species => genus table\n";

my %badGenera = (
  unknown => 1,
  Kosakonia => 0,
  Moheibacter => 0,
  Thermoflexus => 0,
  Chitinivibrio => 0,
  Pluralibacter => 0,
  "Escherichia-Shigella" => 1,
  Ferrovum => 0,
  Planococcaceae => 1,
  "Chthonomonas-Armatimonadetes" => 1,
  Erysipelotrichaceae => 1,
  Peptostreptococcaceae => 1,
  Lelliottia => 0,
  Ruminococcus2 => 1,
  Lachnospiracea => 1,
  Psychrosinus => 0,
  "Armatimonas/Armatimonadetes_gp1" => 1,
  Sporolactobacillaceae => 1,
  "Armatimonas-Armatimonadetes" => 1,
  Desulfosporomusa => 0,
  "Chthonomonas/Armatimonadetes_gp3" => 1,
  );

my %ourGeLiTbl;
my %spParent;
foreach my $file ( @mLiFiles1 )
{
  my @a = split "/", $file;
  my $phGr = pop @a;
  $phGr =~ s/\.lineage$//;
  print "\rProcessing $phGr        ";

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for my $lineage (<IN>)
  {
    chomp $lineage;
    my ($id, $li) = split /\s+/, $lineage;
    my @f = split ";", $li;

    my $sp = pop @f;
    my $ge = pop @f;

    if ( exists $badGenera{$ge} && $badGenera{$ge} == 1 )
    {
      print "\n$file\n";
      print "$li\n";
      exit 1;
    }
    elsif ( exists $badGenera{$ge} && $badGenera{$ge} == 0 )
    {
      $ourGeLiTbl{$ge} = $li;
    }
    $spParent{$sp} = $ge;
  }
  close IN;
}

my %liFile2phGr;
foreach my $file ( @liFiles )
{
  print "\rProcessing $file        ";
  my @a = split "/", $file;
  my $phGr = pop @a;
  $phGr =~ s/\.lineage$//;
  $liFile2phGr{$file} = $phGr;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for my $lineage (<IN>)
  {
    chomp $lineage;
    my ($id, $li) = split /\s+/, $lineage;
    my @f = split ";", $li;

    my $sp = pop @f;

    if ( ! exists $spParent{$sp} )
    {
      my $ge = pop @f;
      my ($g, $suffix) = split "_", @f;
      if ( (defined $suffix && $ge eq $g) )
      {
        print "\nWARNING: $sp not found in spParent tbl - using $g\n";
        $spParent{$sp} = $g;
      }
      elsif ( !defined $suffix )
      {
        print "\nWARNING: $sp not found in spParent tbl - using $ge\n";
        $spParent{$sp} = $ge;
      }
      else
      {
       	warn "\nWARNING: $sp not found in spParent tbl";
        print "Two possible candidates: \n$g\n$ge\n\n";
        exit 1;
      }
    }
  }
  close IN;
}

print "\r                                                      ";

#my $geLiFile = "/Users/pgajer/projects/PECAN/data/microcontax/microcontax_genus_lineage_tbl.txt"; # NOTE: this one has a header !

if (0)
{
  print "--- Appending our genus-lineage to the microcontax pkg genus-lineage tbl\n";
  open OUT, ">>$geLiFile" or die "Cannot open $geLiFile for writing: $OS_ERROR\n";
  for my $ge ( keys %ourGeLiTbl )
  {
    my $li = $ourGeLiTbl{$ge};
    my @f = split ";", $li;
    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;
    print OUT "$ge $fa $or $cl $ph Bacteria\n";

  }
  close OUT;
}

print "--- Parsing microcontax pkg genus-lineage tbl\n";
my ($rgeFaTbl, $rgeOrTbl, $rgeClTbl, $rgePhTbl) = ge_li_tbl( $geLiFile );
my %geFaTbl = %{$rgeFaTbl};



print "--- Checking if all out db genera are found in the master genus-lineage table\n";
my @ges = values %spParent;
@ges = unique( \@ges );
for my $ge ( @ges  )
{
  if ( ! exists $geFaTbl{$ge} )
  {
    #warn "WARNING: $ge not found in geFaTbl\n";
    print "$ge\n";
  }
}

my $spGeFile = $outDir . "/species_genus_tbl_june19_2017.txt";
my @spp = sort keys %spParent;
write_sorted_tbl( \%spParent, \@spp, $spGeFile );

print "\r\n\nSpecies to genus tbl written to $spGeFile\n\n";

exit;

my %geParent;
my %faParent;
my %orParent;
my %clParent;

my %suspeciousSpGePair; # species is not single term name and its first term is not the same as the genus name
my %suspeciousFaName;   # family not ending with 'eae
my %suspeciousOrName;   # order not ending with 'les

my %faEQor;
my %orEQcl;
my %clEQph;


# array and table listing classes with identical taxonomy as their phyla
my @goodCls = ("Thermodesulfobacteria", "Aquificae", "Gemmatimonadetes", "Actinobacteria", "Chrysiogenetes", "Elusimicrobia", "Deferribacteres", "Thermotogae");
my %goodClsTbl = map { $_ => 1 } @goodCls;

# [~/devel/MCextras/data/RDP/rdp_Bacteria_phylum_dir]$ awk 'FS=";" {print $2 " => " "\""$1"\""}' tmp.txt | sort | uniq

my %geParent1 = (
  Spongiispira => "Oceanospirillaceae",
  Salicola     => "Halomonadaceae",
  Vasilyevaea  => "Hyphomicrobiaceae",
  Amorphus     => "Rhodobiaceae",
  Teredinibacter => "Cellvibrionaceae",
  Dasania      => "Spongiibacteraceae",
  Oceanotoga   => "Thermotogaceae",
  Caldithrix   => "DeferribacteralesIS",
  Phocaeicola  => "BacteroidalesIS",
  Fangia       => "ThiotrichalesIS",
  Caedibacter  => "ThiotrichalesIS",
  Thiomonas    => "BurkholderialesIS",
  Aquabacterium => "BurkholderialesIS",
  Aquincola => "BurkholderialesIS",
  Ideonella => "BurkholderialesIS",
  Inhella => "BurkholderialesIS",
  Leptothrix => "BurkholderialesIS",
  Methylibium => "BurkholderialesIS",
  Mitsuaria => "BurkholderialesIS",
  Paucibacter => "BurkholderialesIS",
  Piscinibacter => "BurkholderialesIS",
  Rubrivivax => "BurkholderialesIS",
  Sphaerotilus => "BurkholderialesIS",
  Tepidimonas => "Burkholderiales",
  Tepidimonas => "BurkholderialesIS",
  Thiobacter => "BurkholderialesIS",
  Thiomonas => "BurkholderialesIS",
  Xylophilus => "BurkholderialesIS",
  Microaerobacter => "Bacillaceae",
  Calditerricola => "Bacillaceae",
  Caldalkalibacillus => "Bacillaceae",
);


foreach my $file ( @mLiFiles1 )
{
  my @a = split "/", $file;
  my $phGr = pop @a;
  $phGr =~ s/\.lineage$//;
  print "\rProcessing $phGr        ";

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR";
  for my $lineage (<IN>)
  {
    chomp $lineage;
    my ($id, $li) = split /\s+/, $lineage;
    my @f = split ";", $li;

    my $sp = pop @f;
    my $ge = pop @f;
    my $fa = pop @f;
    my $or = pop @f;
    my $cl = pop @f;
    my $ph = pop @f;

    if ( $fa !~ /eae$/ )
    {
      $suspeciousFaName{$fa}{$li}++;
    }

    if ( $or !~ /les$/ )
    {
      $suspeciousOrName{$or}{$li}++;
    }

    if ( $fa eq $or )
    {
      $faEQor{$fa}{$li}++;
    }

    if ( $or eq $cl )
    {
      $orEQcl{$or}{$li}++;
    }

    if ( $cl eq $ph )
    {
      $clEQph{$cl}{$li}++;
    }

    my @sf = split "_", $sp;
    my $g = $sf[0];
    if ( @sf > 1 && $g ne $ge )
    {
      $suspeciousSpGePair{$sp} = $ge;
    }

    $spParent{$sp} = $ge;
    $geParent{$ge} = $fa;
    $faParent{$fa} = $or;
    $orParent{$or} = $cl;
    $clParent{$cl} = $ph;
  }
  close IN;
}

if ( keys %suspeciousSpGePair > 0 )
{
  print "\n\nSuspecious species-genus pairs\n";
  print_tbl( \%suspeciousSpGePair );
  print "\n\n";
}

if ( keys %suspeciousFaName > 0 )
{
  print "\n\nSuspecious family name\n";
  print_tbl_valued_tbl( \%suspeciousFaName );
  print "\n\n";
}

if ( keys %suspeciousOrName > 0 )
{
  print "\n\nSuspecious order name\n";
  print_tbl_valued_tbl( \%suspeciousOrName );
  print "\n\n";
}

if ( keys %faEQor > 0 )
{
  print "\n\nFamily and order names the same\n";
  print_tbl_valued_tbl( \%faEQor );
  print "\n\n";
}

if ( keys %orEQcl > 0 )
{
  print "\n\nOrder and class names the same\n";
  print_tbl_valued_tbl( \%orEQcl );
  print "\n\n";
}


if ( keys %clEQph > 0 )
{
  print "\n\nClass and phylum  names the same\n";
  print_tbl_valued_tbl( \%clEQph );
  print "\n\n";
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

  my $baseDir  = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/";
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

sub old_get_li_files
{
  #my $baseDirExt = "/Users/pgajer/devel/MCextras/data/vag_exp_V3V4_phGrps_May16_dir/";
  my $baseDirExt = "/Users/pgajer/projects/PECAN/data/phylo_groups/v0.2/";

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
  my @liFiles0 = ("Actinobacteria_dir/Actinobacteria_group_0_V3V4_dir/Actinobacteria_group_0_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_1_V3V4_dir/Actinobacteria_group_1_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_2_V3V4_dir/Actinobacteria_group_2_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_3_V3V4_dir/Actinobacteria_group_3_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_4_V3V4_dir/Actinobacteria_group_4_V3V4.lineage",
		  "Actinobacteria_dir/Actinobacteria_group_5_V3V4_dir/Actinobacteria_group_5_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_0_V3V4_dir/Bacteroidetes_group_0_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_1_V3V4_dir/Bacteroidetes_group_1_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_2_V3V4_dir/Bacteroidetes_group_2_V3V4.lineage",
		  "Bacteroidetes_dir/Bacteroidetes_group_3_V3V4_dir/Bacteroidetes_group_3_V3V4.lineage",
		  "final_small_phyla_V3V4/Chloroflexi_V3V4_dir/Chloroflexi_V3V4.lineage",
		  "final_small_phyla_V3V4/Deinococcus_Thermus_V3V4_dir/Deinococcus_Thermus_V3V4.lineage",
		  "final_small_phyla_V3V4/Fusobacteria_V3V4_dir/Fusobacteria_V3V4.lineage",
		  "final_small_phyla_V3V4/Nitrospirae_V3V4_dir/Nitrospirae_V3V4.lineage",
		  "final_small_phyla_V3V4/Planctomycetes_V3V4_dir/Planctomycetes_V3V4.lineage",
		  "final_small_phyla_V3V4/Spirochaetes_V3V4_dir/Spirochaetes_V3V4.lineage",
		  "final_small_phyla_V3V4/Tenericutes_V3V4_dir/Tenericutes_V3V4.lineage",
		  "final_small_phyla_V3V4/Verrucomicrobia_V3V4_dir/Verrucomicrobia_V3V4.lineage",
		  "final_small_phyla_V3V4/phyla_lessthen_1k_wOG_V3V4_dir/phyla_lessthen_1k_wOG_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_0_V3V4_dir/Firmicutes_group_0_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_1_V3V4_dir/Firmicutes_group_1_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_2_V3V4_dir/Firmicutes_group_2_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_3_V3V4_dir/Firmicutes_group_3_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_4_V3V4_dir/Firmicutes_group_4_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_5_V3V4_dir/Firmicutes_group_5_V3V4.lineage",
		  "Firmicutes_dir/Firmicutes_group_6_V3V4_dir/Firmicutes_group_6_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_0_V3V4_dir/Proteobacteria_group_0_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_10_V3V4_dir/Proteobacteria_group_10_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_11_V3V4_dir/Proteobacteria_group_11_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_12_V3V4_dir/Proteobacteria_group_12_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_13_V3V4_dir/Proteobacteria_group_13_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_14_V3V4_dir/Proteobacteria_group_14_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_15_V3V4_dir/Proteobacteria_group_15_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_17_V3V4_dir/Proteobacteria_group_17_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_1_V3V4_dir/Proteobacteria_group_1_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_2_V3V4_dir/Proteobacteria_group_2_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_3_V3V4_dir/Proteobacteria_group_3_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_4_V3V4_dir/Proteobacteria_group_4_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_5_V3V4_dir/Proteobacteria_group_5_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_6_V3V4_dir/Proteobacteria_group_6_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_7_V3V4_dir/Proteobacteria_group_7_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_8_V3V4_dir/Proteobacteria_group_8_V3V4.lineage",
		  "Proteobacteria_dir/Proteobacteria_group_9_V3V4_dir/Proteobacteria_group_9_V3V4.lineage");

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

sub build_spp_seqID_tree
{
  my ($treeFile, $annFile, $ssTreeFile ) = @_;

  my ($fh, $annFile2) = tempfile("tmp.XXXX", SUFFIX => '.tx', OPEN => 0, DIR => $tmpDir);
  $cmd = "awk '{print \$1\"\\t\"\$2\"__\"\$1}' $annFile > $annFile2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

  $cmd = "rm -f $ssTreeFile; $nw_rename $treeFile $annFile2 | $nw_order -  > $ssTreeFile";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;
}

sub build_spp_tree
{
  my ($treeFile, $annFile, $ssTreeFile ) = @_;

  $cmd = "rm -f $ssTreeFile; $nw_rename $treeFile $annFile | $nw_order -  > $ssTreeFile";
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
  map {print OUT $_ . "\t" . $tbl{$_} . "\n"} @{$r};
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
  my $treeLeavesFile = "$outDir" . "/family_sppSeqIDs.leaves";
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
	my $ogCladeTreeFile = "$outDir/family" . "_clade.tree";
	$cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
	#print "\tcmd=$cmd\n" if $dryRun || $debug;
	system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

	#print "\t--- Extracting leaves of the OG clade\n";
	my $ogCladeTreeLeavesFile = "$outDir/family" . "_clade.leaves";
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
    my $ogCladeTreeFile = "$outDir/family" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$outDir/family" . "_OG_clade.leaves";
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
    my $ogCladeTreeFile = "$outDir/family" . "_OG_clade.tree";
    $cmd = "rm -f $ogCladeTreeFile; $nw_clade $treeFile @og > $ogCladeTreeFile";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

    #print "\t--- Extracting leaves of the OG clade\n";
    my $ogCladeTreeLeavesFile = "$outDir/family" . "_OG_clade.leaves";
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

# extract unique elements from an array
sub unique{

  my $a = shift;
  my %saw;
  my @out = grep(!$saw{$_}++, @{$a});

  return @out;
}

exit 0;
