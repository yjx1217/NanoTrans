#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: plot_isoform_differential_expression.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.09
#  description: plot isoform differential expression
#  example: perl plot_isoform_differential_expression.pl -i Master_Sample_Table.txt -threads 1 -b $batch_id -lfc 1 -p 0.05 
##############################################################


my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $ucsc_dir=$ENV{ucsc_dir};
my $blast_dir=$ENV{blast_dir};
my $bedtools_dir=$ENV{bedtools_dir};
my $bowtie2_dir=$ENV{bowtie2_dir};
my $threads = 1;
my $output_prefix;
my $refseq_fasta_file;
my $refseq_gtf_file;
my $debug = "no";

GetOptions('refseq_fasta_file|f:s' => \$refseq_fasta_file,
	   'refseq_gtf_file|g:s' => \$refseq_gtf_file,
	   'prefix|p:s' => \$output_prefix,
	   'threads|t:i' => \$threads,
	   'debug|d:s' => \$debug);

my $local_time = localtime();
print "\n\n[$local_time] Preparing reference genome for JAFFAL ..\n";

my $base_dir = cwd();
print "base_dir=$base_dir\n";
chdir("$base_dir") or die "cannot change directory to: $base_dir!\n";
system("ln -s $refseq_fasta_file genome.fa");

# add "EN" prefix to genome annotation gtf file
my $refseq_gtf_mod_file = "ref.genome.mod.gtf";
my $refseq_gtf_mod_file_fh = write_file($refseq_gtf_mod_file);

my $refseq_gtf_fh = read_file($refseq_gtf_file);
while (<$refseq_gtf_fh>) {
    chomp;
    if (/^#/) {
	print $refseq_gtf_mod_file_fh "$_\n";
    } elsif (/^\s*$/) {
	print $refseq_gtf_mod_file_fh "$_\n";
    } else {
	my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $_;
	my $mod_attribute = "";
	my @attribute = split ";", $attribute;
	foreach my $a (@attribute) {
	    $a =~ s/^\s+//gi;
	    # print "attribute_a: $a\n";
	    my $feature_id;
	    my $feature_type;
	    if ($a =~ /(gene_id|transcript_id|exon_id) \"([^\"]+)\"/) {
		$feature_type = $1;
		$feature_id = $2;
		# print "feature_id = $feature_id\n";
		if ($feature_id !~ /^EN/) {
		    $feature_id = "EN.${feature_id}";
		}
		if ($feature_id =~ /\_/) {
		    $feature_id =~ s/\_/\./gi;
		}
		# print "mod_feature_id = $feature_id\n";
		$a = "$feature_type \"$feature_id\"";
	    }
	    $mod_attribute .= "$a; ";
	}
	$mod_attribute =~ s/\s+$//gi;
	print $refseq_gtf_mod_file_fh "$chr\t$source\t$feature\t$start\t$end\t$score\t$strand\t$frame\t$mod_attribute\n";
    }
}
close $refseq_gtf_mod_file_fh;

# prepare genome annotation tab file
system("$ucsc_dir/gtfToGenePred -genePredExt  -ignoreGroupsWithoutExons $refseq_gtf_mod_file $output_prefix.ucsc_tab.raw.tab");
my $ucsc_tab_raw = "$output_prefix.ucsc_tab.raw.tab";
my $ucsc_tab_raw_fh = read_file($ucsc_tab_raw);
my $ucsc_tab_tidy = "$output_prefix.tab";
my $ucsc_tab_tidy_fh = write_file($ucsc_tab_tidy);
print $ucsc_tab_tidy_fh "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames\n";
while (<$ucsc_tab_raw_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    print $ucsc_tab_tidy_fh "0\t$_\n";
}

# prepare genome annotation sequence file
my $refseq_fasta_fh = read_file($refseq_fasta_file);
my %refseq_fasta = ();
my @refseq_fasta = ();
parse_fasta_file($refseq_fasta_fh, \%refseq_fasta, \@refseq_fasta);

my $refseq_annotation_fasta_file = "$output_prefix.fa";
my $refseq_annotation_fasta_fh = write_file($refseq_annotation_fasta_file);
my $refseq_annotation_bed_file = "$output_prefix.bed";
my $refseq_annotation_bed_fh = write_file($refseq_annotation_bed_file);
$refseq_gtf_mod_file_fh = read_file($refseq_gtf_mod_file);
my %combined_exon_seq = ();
my @transcript_id = ();
my $combined_exon_seq_id;
my $transcript_id;

while (<$refseq_gtf_mod_file_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $_;
    if ($feature eq "transcript") {
	($transcript_id) = ($attribute =~ /transcript_id \"([^\"]+)\"/);
	push @transcript_id, $transcript_id;
	$combined_exon_seq_id = "genome_annotation_${transcript_id}__range=$chr:${start}-${end}__5'pad=0__3'pad=0__strand=${strand}__repeatMasking=none";
	$combined_exon_seq{$transcript_id}{'exon_seq_id'} = $combined_exon_seq_id;
	$combined_exon_seq{$transcript_id}{'exon_seq'} = "";
    } elsif ($feature eq "exon") {
	my $new_start = $start -1;
	print $refseq_annotation_bed_fh "$chr\t$new_start\t$end\n";

	# ($exon_number) = ($attribute =~ /exon_number \"([^\"]+)\"/);
	my $exon_seq = substr $refseq_fasta{$chr}, $start - 1, $end - $start + 1;
	if ($strand eq "+") {
	    $combined_exon_seq{$transcript_id}{'exon_seq'} =  $combined_exon_seq{$transcript_id}{'exon_seq'}.$exon_seq;
	} else {
	    $exon_seq = revcom($exon_seq);
	    $combined_exon_seq{$transcript_id}{'exon_seq'} =  $combined_exon_seq{$transcript_id}{'exon_seq'}.$exon_seq;
	}
    }
}

foreach my $transcript_id (@transcript_id) {
    print $refseq_annotation_fasta_fh ">$combined_exon_seq{$transcript_id}{'exon_seq_id'}\n$combined_exon_seq{$transcript_id}{'exon_seq'}\n";
}


system("$bedtools_dir/bedtools maskfasta -fi $refseq_fasta_file -fo Masked_genome.fa -bed $output_prefix.bed");
system("$bowtie2_dir/bowtie2-build genome_annotation.fa genome_annotation");
system("$bowtie2_dir/bowtie2-build Masked_genome.fa Masked_genome");
system("$blast_dir/makeblastdb -in $output_prefix.fa -dbtype nucl -out ${output_prefix}_blast");


$local_time = localtime();
print "\n[$local_time] Done\n";



sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    }
    else {
        open($fh, $file) or die "can't open $file";
    }
    return $fh;
}

sub write_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub parse_sample_table {
    my ($fh, $sample_table_hashref, $sample_table_arrayref) = @_;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	    /^\s*$/ and next;
        my ($sample_id, $comparison_group, $replicate_id, $basecalled_fastq_file, $basecalled_fast5_dir, $note) = split /\s+/, $_;
        push @$sample_table_arrayref, $sample_id;
        $$sample_table_hashref{$sample_id}{'comparison_group'} = $comparison_group;
        $$sample_table_hashref{$sample_id}{'replicate_id'} = $replicate_id;
        $$sample_table_hashref{$sample_id}{'basecalled_fastq_file'} = $basecalled_fastq_file;
        $$sample_table_hashref{$sample_id}{'basecalled_fast5_dir'} = $basecalled_fast5_dir;
        $$sample_table_hashref{$sample_id}{'note'} = $note;
    }
}

sub parse_fasta_file {
    my ($fh, $seq_hashref, $seq_arraryref) = @_;
    my $seq_id = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(\S+)/) {
            $seq_id = $1;
            $$seq_hashref{$seq_id} = "";
            push @$seq_arraryref, $seq_id;
        } else {
            $$seq_hashref{$seq_id} .= (uc $_);
        }
    }
}

sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}

sub ensembl_gtf2exon_seq {
}

sub ensembl_gtf2exon_bed {
    my $refseq_gtf_fh =shift @_;
    my @exon_bed = ();
    while (<$refseq_gtf_fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $_;
	if ($feature eq "exon") {
	}
    }
    return @exon_bed;
}
