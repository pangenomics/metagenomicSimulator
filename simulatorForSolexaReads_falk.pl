#!/bin/usr/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw/ tempdir /;
use Bio::SeqIO;
use Bio::Seq::Quality;

sub createAutoInfo($ $ $);
sub readFasta($);


##please set these variables here or using the commandline options
#the folder the binaries are stored in:
my $binfolder = "/g/bork6/mende/scripts/dev/simulator/bin/";
#the folder the genomes are stored in:
my $genomeFolder = "./genomes/";


##variables from input
my $outfolder = "";
#abundanceFile = file with readnumbers for each genomeFile used (Fastafiles are in $genomeFolder)
my $abundanceFile = "";
my $insertSize = 180;
my $insertSD = 0;
my $readLength = 75;
my $insertNumber = 0;
#my $genomeFolder = "/g/bork8/mende/Metagenomics/cMESS/Genomes/All/";

my $help = 0;
my $outprefix = "IlluminaSimReads";
my $qualityFiles = "";
my $genomeInfo = "";


#variables from input for special circumstances -> add new reads to old simulation
my $laneNumber = 1;
my $tileNumber = 1;
my $xCoord = 1;
my $yCoord = 1;
my $multiplexIdentifier = "0";

#variables from input for special circumstances -> folder of insert generator changed
#my $insertGenerator = "/g/bork6/mende/MGSimulation/solexa/version2/programs/generateReads";
#my $insertGenerator = "/g/bork6/mende/MGSimulation/solexa/version2/programs/generateReads";

#internal variables (used by subroutines)
my $currentInsertName = "";
my $currentInsert = "";
my $leftRead = "";
my $rightRead = "";
my $leftReadName = "";
my $rightReadName = "";

my %lengthHash = ();
my %copyNumHash = ();
my %abundanceHash = ();
my %probHash = ();
my %probSumHash = ();
my %newAbundanceHash = ();
my %readAbundanceHash = ();

my $randomFolder = 1;

my @probArray = ();


GetOptions (
     "help|?" => \$help,
     "outfolder=s" => \$outfolder,
     "outprefix=s" => \$outprefix,
     "abundancefile=s" => \$abundanceFile,
     "genomeInfo=s" => \$genomeInfo,
     "genomefolder=s" => \$genomeFolder,
     "readlength=i" => \$readLength,
     "insertSize=i" => \$insertSize,
     "insertSD=i" => \$insertSD,
     "insertNumber=i" => \$insertNumber,
     "binfolder=s" => \$binfolder,
     "lanenumber=i" => \$laneNumber,
     "tilenumber=i" => \$tileNumber,
     "xcoord=i" => \$xCoord,
     "ycoord=i" => \$yCoord,
     "multiplexidentifier=s" => \$multiplexIdentifier,
     "qualityfiles=s" => \$qualityFiles,
     "randomFolder=i" => \$randomFolder
     );

pod2usage(-verbose => 2) if $help;

if ($outfolder eq ""){
print("outfolder not specified, please run simulatorForSolexaReads.pl --help for more information\n");
exit;
}
if ($abundanceFile eq ""){
print("abundanceFile not specified, please run simulatorForSolexaReads.pl --help for more information\n");
exit;
}
if ($genomeFolder eq ""){
print("genomeFolder not specified, please run simulatorForSolexaReads.pl --help for more information\n");
exit
}
if ($qualityFiles eq ""){
print("qualityFiles not specified, please run simulatorForSolexaReads.pl --help for more information\n");
exit
}
if ($genomeInfo eq ""){
$genomeInfo  = $abundanceFile ;
$genomeInfo  =~ s/\.[^\.]+$/_ginf.txt/;
print("genomeInfo not specified, creating best guess file (1 genome per file)..\nFile (reuseable): $genomeInfo  \n");
createAutoInfo($abundanceFile,$genomeInfo,$genomeFolder);
}


my $rawreadsname = "rawReadsXXXX";

my $rawreadsfolder = $outfolder;
if ($randomFolder ==1){
     $rawreadsfolder = tempdir( $rawreadsname, DIR => $outfolder);
}

#my $rawreads = tempdir ( $rawreadsname, DIR => $outfolder, CLEANUP => 1);

open(GENOMEINFO, "<$genomeInfo") or die "Cannot open the genome information file: $genomeInfo\n";

while(my $currentLine = <GENOMEINFO>){
     chomp($currentLine);
     my @currentLineArray = split(/\t/,$currentLine);
     my $genomeID = $currentLineArray[0];
     my $genomeLength = $currentLineArray[1];
     my $genomeCopyNumber = $currentLineArray[2];

     $lengthHash{$genomeID} = $genomeLength;
     $copyNumHash{$genomeID} = $genomeCopyNumber;
}
#read in organism abundance values
open(ABUNDANCES, "<$abundanceFile") or die "Cannot open the genome abundance file: $abundanceFile\n";

while(my $currentLine = <ABUNDANCES>){
     chomp($currentLine);
     my @currentLineArray = split(/\t/,$currentLine);
     my $genomeID = $currentLineArray[0];
     my $abundance = $currentLineArray[1];
     my $genomeLength = $lengthHash{$genomeID};
     my $genomeCopyNumber = $copyNumHash{$genomeID};

     $abundanceHash{$genomeID} = $abundance*$genomeLength*$genomeCopyNumber;
}

#calculate total organism abundance
my $totalAbundance = 0;
for my $currentGenome (keys %abundanceHash ) {
     $totalAbundance += $abundanceHash{$currentGenome};
}

#find right multiplicator
my $multiplicator = 1;
my $newTotalAbundance = $totalAbundance;
while ($newTotalAbundance < $insertNumber){
     $newTotalAbundance = $newTotalAbundance * $multiplicator;
     $multiplicator ++;
}

my $currentGenome ="";

#apply multiplicator to all organism abundances
for $currentGenome (keys %abundanceHash ) {
     $newAbundanceHash{$currentGenome} = int($abundanceHash{$currentGenome}*$multiplicator+0.5);
}

#fill probArray: an array to randomly pick a which genome a read comes from according to the genomes abundances and sequencelengths.
$newTotalAbundance = 0;
my @probArrayH; my @probArrayL; #high & low bond of probability array.. anyway a sorted list and the old vec approach is killing mem
my @probArrayID; my $cnt=0;
for $currentGenome (keys %abundanceHash ) {

     my $currentNewAbundance = $newAbundanceHash{$currentGenome};
     $probArrayID[$cnt] = $currentGenome;
         $probArrayL[$cnt] = $newTotalAbundance;
         $probArrayH[$cnt] = $newTotalAbundance+$currentNewAbundance-1;

#    for my $j ($newTotalAbundance .. $newTotalAbundance+$currentNewAbundance-1){
#            $probArray[$j] = $currentGenome;
#        }
     $newTotalAbundance += $currentNewAbundance;
         $cnt++;
}

#calculate how many reads come from which genome.
for (my $i = 1; $i <= $insertNumber; $i++){
     my $currentIndex = int(rand($newTotalAbundance));
     #$currentGenome = $probArray[$currentIndex];
     #slower but much more memory efficient..
     my $cnt=0;
     while ($currentIndex > $probArrayH[$cnt]){$cnt++;}
     $currentGenome = $probArrayID[$cnt];

     $readAbundanceHash{$currentGenome} += 1;
}

#report how many reads are generated from which species
#should be disabled later on maybe
for my $currentGenome (keys %readAbundanceHash ) {
     print $currentGenome."\t".$readAbundanceHash{$currentGenome}."\n";
}


my $leftends = "$rawreadsfolder\/$outprefix\.1";
my $rightends = "$rawreadsfolder\/$outprefix\.2";

open(LEFTENDS, ">$leftends") or die "Cannot open the data file for data output: $leftends\n";
open(RIGHTENDS, ">$rightends") or die "Cannot open the data file for data output: $rightends\n";

#while(my $currentLine = <ABUNDANCES>){
#    chomp($currentLine);
#    my @currentLineArray = split(/\t/,$currentLine);
#    my $genomeID = $currentLineArray[0];
#    my $currentReadNumber = $currentLineArray[1];
#    my $genomeFile = $genomeFolder."\/".$genomeID;

#    my $insertnumber = int(($currentReadNumber+1)/2);

for my $currentGenome (keys %readAbundanceHash ) {
     my $insertnumber = $readAbundanceHash{$currentGenome};
     my $genomeFile = $genomeFolder."\/".$currentGenome;

     #rightends and leftends file
     #print STDERR "$binfolder\/generateReads_d -i $genomeFile -l $insertSize -c $insertnumber|\n";
     open PIPE, "$binfolder\/generateReads -i $genomeFile -l $insertSize -c $insertnumber -s $insertSD |" or die "Cannot open pipe from $binfolder\/generateReads: $!";

     $currentInsertName = <PIPE>;
     $currentInsert = "";
     $leftRead = "";
     $rightRead = "";
     $leftReadName = "";
     $rightReadName = "";

     while ( my $currentFastaLine=<PIPE> ){
         chomp($currentFastaLine);
         if (substr($currentFastaLine,0,1) eq "\>") {
             &getReadsFromInsert($currentInsert);
             &getReadNames($currentGenome);
             print LEFTENDS "$leftReadName\n";
             print LEFTENDS "$leftRead\n";

             print RIGHTENDS "$rightReadName\n";
             print RIGHTENDS "$rightRead\n";

             #cleanup and start next insert
             $currentInsertName = $currentFastaLine;
             $currentInsert = "";
             $leftRead = "";
             $rightRead = "";
         }

         else{
             $currentInsert = $currentInsert.$currentFastaLine;
         }
     }
     &getReadsFromInsert($currentInsert);
     &getReadNames($currentGenome);
     print LEFTENDS "$leftReadName\n";
     print LEFTENDS "$leftRead\n";

     print RIGHTENDS "$rightReadName\n";
     print RIGHTENDS "$rightRead\n";

}
close(RIGHTENDS);
close(LEFTENDS);
close(ABUNDANCES);

my $qualityname = "qualitiesXXXX";
my $qualityfolder = tempdir( $qualityname, DIR => $outfolder);
#my $qualityfolder = tempdir( $qualityname, DIR => $outfolder, CLEANUP => 1);


system "$binfolder\/simulateSequencingErrors --outfasta=$qualityfolder\/$outprefix\.1.fasta --outqual=$qualityfolder\/$outprefix\.1.fasta.qual -q $qualityFiles $leftends";
system "$binfolder\/simulateSequencingErrors --outfasta=$qualityfolder\/$outprefix\.2.fasta --outqual=$qualityfolder\/$outprefix\.2.fasta.qual -q $qualityFiles $rightends";



#convert fasta + qual to fastq in illumina scale
my $fastqname = "fastqXXXX";
my $fastqfolder = tempdir( $fastqname, DIR => $outfolder);
#my $fastqfolder = tempdir( $fastqname, DIR => $outfolder, CLEANUP => 1);

open(FASTQFILE, ">$fastqfolder\/$outprefix\.1.fq") or die "Cannot open the data file for data output: $fastqfolder\/$outprefix\.1.fq";
close(FASTQFILE);
open(FASTQFILE, ">$fastqfolder\/$outprefix\.2.fq") or die "Cannot open the data file for data output: $fastqfolder\/$outprefix\.2.fq";
close(FASTQFILE);

my $fasta_1_obj = Bio::SeqIO->new( -file      => "$qualityfolder\/$outprefix\.1.fasta",
                                    -format    => 'fasta',
                                    -variant   => 'sanger' );

my $qual_1_obj = Bio::SeqIO->new(  -file      => "$qualityfolder\/$outprefix\.1.fasta.qual",
                                    -format    => 'qual',
                                    -variant   => 'sanger' );


my $fastq_1_obj = Bio::SeqIO->new(   -file    => ">$fastqfolder\/$outprefix\.1.fq",
                                      -format  => 'fastq',
                                      -variant => 'illumina' );

while (my $fasta_1_obj_entry = $fasta_1_obj->next_seq){
   ## create objects for both a seq and its associated qual
   my $qual_1_obj_entry = $qual_1_obj->next_seq;

   die "Fasta and Quality file are not synchronized!\n" unless $fasta_1_obj_entry->id eq $qual_1_obj_entry->id;

   ## Here we use seq and qual object methods feed info for new BSQ
   ## object.
   my $bioSeqQual_1_obj = Bio::Seq::Quality-> new( -id   => $fasta_1_obj_entry->id,
                                                   -seq  => $fasta_1_obj_entry->seq,
                                                   -qual => $qual_1_obj_entry->qual );

   ## and print it out.
   $fastq_1_obj->write_fastq($bioSeqQual_1_obj);
}

  my $fasta_2_obj = Bio::SeqIO->new( -file     => "$qualityfolder\/$outprefix\.2.fasta",
                                    -format    => 'fasta',
                                    -variant   => 'sanger' );

my $qual_2_obj = Bio::SeqIO->new(  -file      => "$qualityfolder\/$outprefix\.2.fasta.qual",
                                    -format    => 'qual',
                                    -variant   => 'sanger' );

my $fastq_2_obj = Bio::SeqIO->new(   -file    => ">$fastqfolder\/$outprefix\.2.fq",
                                      -format  => 'fastq',
                                      -variant => 'illumina' );

while (my $fasta_2_obj_entry  = $fasta_2_obj->next_seq){
   ## create objects for both a seq and its associated qual
   my $qual_2_obj_entry = $qual_2_obj->next_seq;

   die "Fasta and Quality file are not synchronized!\n" unless 
$fasta_2_obj_entry->id eq $qual_2_obj_entry->id;

   ## Here we use seq and qual object methods feed info for new BSQ
   ## object.
   my $bioSeqQual_2_obj = Bio::Seq::Quality-> new( -id   => $fasta_2_obj_entry->id,
                                                   -seq  => $fasta_2_obj_entry->seq,
                                                   -qual => $qual_2_obj_entry->qual );

   ## and print it out.
   $fastq_2_obj->write_fastq($bioSeqQual_2_obj);
}

#system "/g/bork2/sunagawa/bin/convert_project -f fasta -t fastq $qualityfolder\/$outprefix\.1.fasta $fastqfolder\/$outprefix\.1";
#system "/g/bork2/sunagawa/bin/convert_project -f fasta -t fastq $qualityfolder\/$outprefix\.2.fasta $fastqfolder\/$outprefix\.2";

#system "perl /g/bork6/mende/MGSimulation/solexa/version2/programs/fastq_qual_converter.pl $fastqfolder\/$outprefix\.1.fastq";
#system "perl /g/bork6/mende/MGSimulation/solexa/version2/programs/fastq_qual_converter.pl $fastqfolder\/$outprefix\.2.fastq";


#this needs to be done outside right now, needs rescripting later
#system "perl -i.orig -lane '$line = ($. / 2); $line2 = $. / 4; if ($line !~ /\D/ && $line2 =~ /\D/){$_ =~ tr/[a-z]/[A-Z]/; $_ =~ s/[^ATGC]/N/g; print $_}else{print $_}' $fastqfolder\/$outprefix\.1\.fq;";
#system "perl -i.orig -lane '$line = ($. / 2); $line2 = $. / 4; if ($line !~ /\D/ && $line2 =~ /\D/){$_ =~ tr/[a-z]/[A-Z]/; $_ =~ s/[^ATGC]/N/g; print $_}else{print $_}' $fastqfolder\/$outprefix\.2\.fq;";




#generates 2 reads from an insert, decides which strand the reads are from, and stores the sequence in global variables $leftRead, $rightRead and $strand;
sub getReadsFromInsert{
     my $insertSequence = shift;
     my $leftEnd = "";
     my $rightEnd = "";
     $leftEnd = substr($insertSequence, 0, $readLength);
     $rightEnd = 
substr($insertSequence,length($insertSequence)-$readLength,$readLength);
     $rightEnd = reverse($rightEnd);

     # if probabilites for each strand are different then 0.5 to 0.5 change next line
     my $strand = int(rand(2));

     if ($strand == 0) {
         $rightEnd =~ tr/ACGTacgt/TGCAtgca/;
         $leftRead = $leftEnd;
         $rightRead = $rightEnd;
     }
     elsif ($strand == 1){
         $rightEnd =~ tr/ACGTacgt/TGCAtgca/;
         $leftRead = $rightEnd;
         $rightRead = $leftEnd;
     }
     else{
         print "Random number for strand buggy";
     }
}

#generate names for the reads according to the Illumina fastq naming convention
sub getReadNames{
     my $genomeID = shift;
     $leftReadName = ">$genomeID\:$laneNumber\:$tileNumber\:$xCoord\:$yCoord\#$multiplexIdentifier\/1";
     $rightReadName = ">$genomeID\:$laneNumber\:$tileNumber\:$xCoord\:$yCoord\#$multiplexIdentifier\/2"; 


     if ($yCoord < 10000){
         $yCoord++;
     }
     elsif ($xCoord < 10000){
         $yCoord = 1;
         $xCoord++;
     }
     elsif ($tileNumber < 100){
         $xCoord = 1;
         $tileNumber++;
     }
     elsif ($laneNumber <= 8){
         $tileNumber = 1;
         $laneNumber++;
     }
}



sub createAutoInfo($ $ $){
   my ($abF,$outF,$genoP) = @_;
   open I,"<",$abF or die $!;
   open O,">",$outF or die $!;

   while (my $line = <I>){
     chomp($line); next if ($line =~ m/^#/); next if (length($line) < 2); my @spl = split("\t",$line);
     my $tmpFile = $genoP."/".$spl[0];
     my $href = readFasta($tmpFile); my %fash = %{$href};
     my $len = 0;
     foreach my $ke (keys %fash){
       $len += length($fash{$ke});
     }
     print O $spl[0]."\t".$len."\t1\n";
   }
   close O; close I;
}

sub readFasta($){
   my ($fil) = @_;
   open(FAS,"<","$fil") || die("Couldn't open FASTA file $fil .");
      my %Hseq;
      my $temp;
      my $line; my $hea=<FAS>; chomp ($hea);
       my $trHe = ($hea);
       # get sequence
     while($line = <FAS>){
       #next if($line =~ m/^;/);
       if ($line =~ m/^>/){
         chomp($line);
         $Hseq{$trHe} = $temp;
         $trHe = ($line);
         $temp = "";
         next;
       }
     #chomp($line);
     #$line =~ s/\s//g;
     $temp .= ($line);
     }
     $Hseq{$trHe} = $temp;
     return \%Hseq;
}




=head1 NAME

iMESSi - Simulator for Illumina metagenomic sequences

=head1 SYNOPSIS

simulatorForSolexaReads.pl [options]

=head1 OPTIONS

=over 4

=item B<--help>

     Prints this help message

=item B<--outfolder>

     Folder for output

=item B<--outprefix>

     Prefix of output filenames (default: SolexaSimReads)

=item B<--abundanceFile>

     input file with genomes and their organism number

=item B<--genomeInfo>

     file with information about the genomes: genome length, copy number in an organism.

=item B<--genomeFolder>

     folder the genomes mentioned in abundanceFile are at.

=item B<--readlength>

     Read length of simulated reads.

=item B<--insertNumber>

     Number of inserts to simulate, from each insert 2 paired-end reads are generated.

=item B<--insertSize>

     Mean insert size of simulated reads.

=item B<--insertSD>

     Standard deviation from the mean insert size.

=item B<--qualityfiles>

     location of quality files for sequencing error and quality generation. Note that the first file can be given by just specifying it's location and multiple files can be given using folling scheme --qualityfiles="file1.qual -q file2.qual -q file3.qual".

=item B<--binfolder>

     Location of generateReads and simulateSequencingErrors, has to be set here or inside the script

=item B<Options to add more data to old simulations>

     --multiplexIdentifier    Identifier for multiplex Libraries (default is 0)

     --laneNumber        Starting lane number

     --tileNumber        Starting tile number

     --xCoord        Starting xCoord number

     --yCoord        Starting yCoord number

=back

=head1 DESCRIPTION

This program is used to simulate Illumina reads from a specified genome or metagenome

=head1 AUTHOR

Daniel Mende, E<lt>mende@embl.deE<gt> et al. Please cite us: .

=cut

