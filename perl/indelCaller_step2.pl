#!/usr/bin/env perl

########## LICENCE ##########
# Copyright (c) 2020-2021 Genome Research Ltd
# 
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
# 
# This file is part of NanoSeq.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
###########################
# * Perl script
#     * Read all the indels, group by read bundle, mark sites in read bundles where we would believe an indel
#         * Add for each site if it is masked or not - when a final event overlaps a masked site, remove
#         * Reformat read bundle names to: RB:Z:1,11373,11754,GCT,TTA
#     * Extract read bundles from bam to create mini bam
#         * Remove duplicate flag
#     * samtools + bcftools to call indels
#     * read vcf output and decide which are ok and if no indel is called, warn!

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Which;
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use Capture::Tiny qw(capture);

my %opts;
GetOptions('b|bam=s'    => \my $botseq_bam_file,
           'r|ref=s'    => \my $ref_genome,
           'o|out=s'    => \$opts{'o'},
           's|sample=s' => \my $sample_name,
           'k|keep'     => \$opts{'k'},
           't|sort'     => \$opts{'t'},
           'i|index'    => \$opts{'i'},
           'h|help'     => \$opts{'h'}
) or pod2usage(2);
pod2usage(-verbose => 1,  -exitval => 0) if(defined $opts{'h'});
pod2usage(2) if( @ARGV != 1 );
die ("samtools not found in path\n") unless ( which 'samtools' );
die ("bcftools not found in path\n") unless ( which 'bcftools' );

die ("\nOutput prefix not defined\n") unless(defined $opts{'o'});
my $out_name = basename($opts{'o'});
my $out_dir = dirname($opts{'o'});
die ("\nMust define the reference\n") unless ( $ref_genome);
die ("\nReference $ref_genome not found\n") unless ( -e $ref_genome );

die ("\nMust a BAM file\n") unless ( $botseq_bam_file );
die ("\nBAM / CRAM $botseq_bam_file not found\n") unless ( -e $botseq_bam_file );
my ($ext) = $botseq_bam_file =~ /(\.[^.]+)$/;
$ext =~ s/.$/i/;
die ("\nBAM / CRAM index for $botseq_bam_file not found\n") unless ( -e $botseq_bam_file . $ext);

$sample_name = "sample_1" unless ( $sample_name);

my $FILE = $ARGV[0];
die ("\nInput file $FILE not found\n") unless ( -e $FILE );
open(IN, "zcat $FILE |") or die( "\nProblem with gunzip $FILE\n" );

my $tempdir;
if ( defined $opts{'k'} ){
  $tempdir = tempdir("tmp.XXXXXXXX", DIR => $out_dir, CLEANUP=> 0);
} else {
  $tempdir = tempdir("tmp.XXXXXXXX", DIR => $out_dir, CLEANUP=> 1);
}

#chromosomes in biological order
sub cmp_chr {
  my $ca = $a;
  my $cb = $b;
  $ca=~s/^chr//;
  $cb=~s/^chr//;
  if ($ca !~ /\D/ ) {
    if ($cb !~ /\D/ ) {
      return ( $ca <=> $cb);
    }
    else { return -1; }
  }
  elsif ($cb !~ /\D/) { return 1;}
  elsif (($ca eq "M" or $ca eq "MT") and ( $cb eq "X" or $cb eq "Y") ) { return 1;}
  elsif (($ca eq "X" or $ca eq "Y" ) and ( $cb eq "M" or $cb eq "MT")) { return -1;}
  elsif ( $ca eq "X" ) {return -1;}
  elsif ( $cb eq "X" ) {return 1;}
  elsif ( $ca eq "Y" ) {return -1;}
  elsif ( $cb eq "Y" ) {return 1;}
  elsif ( $ca eq "M" or $ca eq "MT") {return -1;}
  elsif ( $cb eq "M" or $cb eq "MT") {return 1;}
  else { return $ca cmp $cb;}
}

sub runCmd {
  my $cmd = $_[0];
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };
  die "Error calling $cmd, $stderr\n" if ( $exit != 0 );
}

# Load indels:
my %indels;
my $counts = 0;
while(<IN>) { # input bed files are generated by indelCaller.pl [step 1]

  next if ( /BULK_SEEN/ ); #filter if site comes from the bulk sample
  chomp;
  my($chr,$pos0,$pos1,$info)             = (split(/\t/,$_  ))[0,1,2,3]    ;
  my($rb_id,$dp,$qpos,$context,$sw,$snp) = (split(/;/,$info))[0,1,2,3,4,5];
  my @tmp = split(/[:\-\|]/,$rb_id);
  $"      = ",";
  $rb_id  = "@tmp";
  $sw     =~ s/SW=//;
  $snp    =~ s/cSNP=//;
  $indels{$rb_id}->{$pos1}->{"dp"     }  = $dp     ;
  $indels{$rb_id}->{$pos1}->{"qpos"   }  = $qpos   ;
  $indels{$rb_id}->{$pos1}->{"context"}  = $context;
  $indels{$rb_id}->{$pos1}->{"sw"     }  = $sw     ;
  $indels{$rb_id}->{$pos1}->{"snp"    }  = $snp    ;
  $counts++;
}
print STDOUT "$counts indel sites seen\n";
print STDOUT scalar(keys(%indels)), " readbundles with indels\n";

open(OUT, "| bcftools sort -Ov -T $tempdir -o $tempdir/$out_name.vcf -") || die "Failed to write to $out_name.final.vcf\n";
my $header = "";
my $vcf_content = "";
my $once = 1;
print STDOUT "\ntmp dir : $tempdir\n\n";
foreach my $rb_id (keys %indels) {
  my($chr,$start,$end) = (split(/,/,$rb_id))[0,1,2];
  my(%good_sites,%bad_sites);
  foreach my $pos ( keys %{$indels{$rb_id}} ) {
    if($indels{$rb_id}->{$pos}->{"sw"} == 0 && $indels{$rb_id}->{$pos}->{"snp"} == 0) {
      $good_sites{$pos} = 1;
    } else {
      $bad_sites{$pos} = 0;
    }  
  }
  if(scalar(keys(%good_sites)) == 0) {
    next; # all indels overlap noisy or common SNP sites
  }
  $"=",";
  my @tmp = keys %{$indels{$rb_id}};
  my $positions = "@tmp";
  print STDOUT "processing read bundle $rb_id ...\n"; 
  print STDOUT "Step 1...\n";
  open(BAM_OUT," | samtools view -bo $tempdir/$out_name.tmp.bam -") or die("\nError running: samtools view -bo $tempdir/$out_name.tmp.bam - :$!\n");
  open(BAM,"samtools view -h $botseq_bam_file $chr:$start-$end | ") or die( "\nError running : samtools view -h $botseq_bam_file $chr:$start-$end : $!\n");
  while(<BAM>) {
    next unless ( /(RB:Z:$rb_id)|(^@)/);
    if ( /^@/ ) {
      print BAM_OUT $_;
    } else {
      my @fields = (split(/\t/,$_));
      $fields[1] = $fields[1] - 1024 if($fields[1] > 1024); #remove duplicate  flag
      print BAM_OUT join "\t",@fields; 
    }
  }
  close( BAM_OUT);

  print STDOUT "Step 2...\n";
  &runCmd( "samtools index $tempdir/$out_name.tmp.bam" );

  print STDOUT "Step 3...\n";
  #`samtools mpileup --no-BAQ  -d 250 -m 2 -F 0.5 -r $chr:$start-$end --BCF --output-tags DP,DV,DP4,SP -f $ref_genome -o $tempdir/$out_name.bcf $tempdir/$out_name.tmp.bam`;
  &runCmd("bcftools mpileup --no-BAQ  -L 250 -m 2 -F 0.5 -r $chr:$start-$end -O b -a DP,DV,DP4,SP -f $ref_genome -o $tempdir/$out_name.bcf $tempdir/$out_name.tmp.bam");
  
  print STDOUT "Step 4...\n";
  &runCmd( "bcftools index -f $tempdir/$out_name.bcf $tempdir/$out_name.indexed.bcf");

  print STDOUT "Step 5...\n";
  #`bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v $tempdir/$out_name.bcf -o $tempdir/$out_name.tmp.vcf`;
  &runCmd("bcftools call --ploidy 1 --skip-variants snps --multiallelic-caller --variants-only  -O v $tempdir/$out_name.bcf -o $tempdir/$out_name.tmp.vcf"); 

  print STDOUT "Step 6...\n";
  &runCmd("bcftools norm -f $ref_genome $tempdir/$out_name.tmp.vcf > $tempdir/$out_name.tmp2.vcf"); 

  #print "$rb_id:$positions:\n";
  my $get_header = 0;
  if($header eq "") {
    $get_header = 1;
  }
  print STDOUT "\nChecking result\n";
  open(IN_TMP, "<$tempdir/$out_name.tmp2.vcf") or die "Couldn't open $tempdir/$out_name.tmp2.vcf\n";
  my @info = ();
  my $count = 0;
  my $best_qual = -10;
  my $best_count = 0;
  while(<IN_TMP>) {
    if(/^#/) {
      if($get_header == 1) {
        #add info and filter fields to header
        if ($once and /^##FORMAT/){
          $_ = "##FILTER=<ID=MASKED,Description=\"Site overlaps with SW or SNP site\">\n" . $_;
          $_ = "##INFO=<ID=RB,Number=1,Type=String,Description=\"Readbundle ID\">\n" . $_ ;
          $once = 0;
        }
        if ( /^#CHROM/ ){ #use a simple sample name, without path
          chomp;
          my @fields = split(/\t/,$_);
          pop(@fields);
          push(@fields, "$sample_name\n");
          $_ = join "\t", @fields;
        }
        #clean up these lines so no tmp files appear
        if (/^##bcftoolsCommand=mpileup/) {
          $_ ="##bcftoolsCommand=mpileup --no-BAQ -L 250 -m 2 -F 0.5 -O b -a DP,DV,DP4,SP\n";
        }
        if (/^##bcftools_callCommand=call/) {
          $_ ="##bcftools_callCommand=call --ploidy 1 --skip-variants snps --multiallelic-caller --variants-only -O v\n";
        }
        if (/^##bcftools_normCommand=norm/) {
          $_ ="##bcftools_normCommand=norm\n";
        }
        print OUT $_;
        $header .= $_;
      }
      next;
    }
    print "RESULT: $_\n";
    my @fields = split(/\t/,$_);
    $fields[7] = "$fields[7];RB=$rb_id";
    $fields[6] = "PASS";
    # If there is overlap with SW or cSNP, flag it as MASKED
    for(my $i=$fields[1]; $i<= ($fields[1] + length($fields[3])); $i++) {
      if(exists($indels{$rb_id}->{$i}) && ($indels{$rb_id}->{$i}->{"sw"} == 1 || $indels{$rb_id}->{$i}->{"snp"} == 1)) {
        $fields[6] = "MASKED";
      }
    }
    $info[$count] = join "\t",@fields;
    if($fields[5] > $best_qual) {
      $best_qual = $fields[5];
      $best_count = $count;
    }
    $count++;
  }
  close(IN_TMP);
  print OUT $info[$best_count] if (  @info > 0 );
}
if (  $header eq "") { #minimal header for an empty vcf
  print OUT "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n"; 
  OUT->flush();
}
close OUT or die("Error closing output pipe to $tempdir/$out_name.vcf.gz\n");

if ( defined $opts{'t'}) {
  #generate a new header by reordering the reference contigs
  open(REFI,"<".$ref_genome.".fai") or die ("Couldn't open $ref_genome.fai\n");
  my %chrdic = ();
  while(<REFI>){
    my ($ichr,$ilength) = split;
    $chrdic{ $ichr} = $ilength;
  }
  close(REFI);
  my @sortedChr = sort cmp_chr (keys %chrdic);
  my $newHead = "";
  foreach my $ichr ( @sortedChr ) {
    $newHead.="##contig=<ID=".$ichr.",length=".$chrdic{$ichr}.">\n";
  } 
  open(VCFI,"<$tempdir/$out_name.vcf") or die ("Couldn't open $tempdir/$out_name.vcf\n");
  open(VCFO,">$tempdir/$out_name.sorted.vcf") or die ("Couldn't write to $out_dir/$out_name.sorted.vcf\n");
  my $writeNewHead = 1;
  while(<VCFI>){
    if ( /^##contig=/) {
      if ( $writeNewHead == 1 ) {
        $writeNewHead = 0;
        $_ = $newHead;
      }
      else {
        next;
      }
    }
    print VCFO $_;
  }
  close(VCFO);
  close(VCFI);
  &runCmd("bcftools sort -Oz -T $tempdir -o $out_dir/$out_name.vcf.gz $tempdir/$out_name.sorted.vcf");
} else {
  &runCmd("bcftools view -Oz -o $out_dir/$out_name.vcf.gz $tempdir/$out_name.vcf");
}

if ( defined $opts{'i'} ){
  &runCmd("bcftools index -f -t $out_dir/$out_name.vcf.gz");
}

__END__

=head1 NAME

indelCaller_step2.pl - Filter sites of interest to look for indels.

=head1 SYNOPSIS

indelCaller_step2.pl  [options] -r ref -b BAM -o prefix input.bed.gz

    -out               -o   Output prefix
    -ref               -r   reference
    -bam               -b   bam / cram file

  Optional parameters:
    -keep              -k   Keep intermediate files
    -index             -i   Create vcf index
    -sample            -s   Sample name ( sample_1 )
    -sort              -t   Put contigs in normal order
    -help              -h

=head1 OPTIONS

=over 8

=item B<-out>

Output prefix. 
Final ouptput is a fileterd bed file.

=item B<-ref>

Reference file

=item B<-bam>

Original NanoSeq BAM (CRAM) for sample.

=item B<-keep>

Do not remove teporary intermediate directory.

=item B<-index>

Create index for output VCF.

=item B<-sample>

Sample name for the VCF output. (sample_1)

=item B<-sort>

Resort cotigs from vcf to be in normal biological order.

=back

=cut
