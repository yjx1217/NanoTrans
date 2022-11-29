#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(first);

##############################################################
#  script: tidy_jaffal_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.12
#  description: tidy up JAFFAL's output files
#  example: perl tidy_jaffal_output.pl -icsv input.csv(.gz) -ifa input.fasta(.gz) -p output_prefix -x transcript2gene_map
##############################################################

my ($input_csv, $input_fasta, $output_prefix, $u2t_conversion, $transcript2gene_map);
GetOptions('input_csv|icsv:s' => \$input_csv, 
	   'input_fasta|ifa:s' => \$input_fasta,
	   'output_prefix|p:s' => \$output_prefix,
	   'u2t_conversion|u2t:s' => \$u2t_conversion,
           'transcript2gene_map|x:s' => \$transcript2gene_map);



my $transcript2gene_map_fh = read_file($transcript2gene_map);
my %transcript2gene_map = ();
my %gene_id2name = ();
while (<$transcript2gene_map_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    /^transcript_id\tgene_id\tgene_name/ and next;
    my ($transcript_id, $gene_id, $gene_name) = split /\t/, $_;
    $transcript2gene_map{$transcript_id}{'gene_id'} = $gene_id;
    $transcript2gene_map{$transcript_id}{'gene_name'} = $gene_name;
    if (not exists $gene_id2name{$gene_id}) {
        $gene_id2name{$gene_id} = $gene_name;
    }
}

my $input_csv_fh = read_file($input_csv);
my $output_tsv = "$output_prefix.gene_fusion.transcripts.txt";
my $output_tsv_fh = write_file($output_tsv);

print $output_tsv_fh "sample_id\tfusion_gene_by_name\tfusion_gene_by_id\tpart1_chr\tpart1_breakpoint\tpart1_strand\tpart2_chr\tpart2_breakpoint\tpart2_strand\tbreakpoint_distance_by_kb(for_intrachromosomal_breakpoints_only)\tspanning_read_pairs\tspanning_reads\tinframe\taligns\trearrangement\tread_id\tread_break\tclassification\tknown\n";

while (<$input_csv_fh>) {
    chomp;
    /^\s*$/ and next;
    /^\s*#/ and next;
    /^sample,fusion genes/ and next;
    my @line = split /,/, $_;
    my $sample_id;
    if ($line[0]=~ /basecalled_reads.Q/) {
	($sample_id) = ($line[0] =~ /(\S+)(?:\.basecalled_reads.Q.*.pass.fastq|\.fq)$/);
    } else {
	$sample_id = $line[0];
    }
    # print "line[0]=$line[0], sample_id = $sample_id\n";
    my $fusion_gene_by_id = $line[1];
    my @fusion_gene_by_id = split /:/, $fusion_gene_by_id;
    my @fusion_gene_by_id_tidy = ();
    my @fusion_gene_by_name = ();
    foreach my $gene_id (@fusion_gene_by_id) {
	$gene_id =~ s/^EN\.//gi;
	my $gene_name = $gene_id;
	if (exists $gene_id2name{$gene_id}) {
	    $gene_name = $gene_id2name{$gene_id};
	}
	push @fusion_gene_by_id_tidy, $gene_id;
	push @fusion_gene_by_name, $gene_name;
    }
    my $fusion_gene_by_id_tidy = join ":", @fusion_gene_by_id_tidy;
    my $fusion_gene_by_name = join ":", @fusion_gene_by_name;
    my @other_info = splice @line, 2;
    my $other_info = join "\t", @other_info;
    print $output_tsv_fh "$sample_id\t$fusion_gene_by_name\t$fusion_gene_by_id_tidy\t$other_info\n";
}

my $input_fasta_fh = read_file($input_fasta);
my %input_fasta = ();
my @input_fasta = ();
parse_fasta_file($input_fasta_fh, \%input_fasta, \@input_fasta);

my $output_fasta = "$output_prefix.gene_fusion.transcripts.fa";
my $output_fasta_fh = write_file($output_fasta);
foreach my $id (@input_fasta) {
    my ($sample_id_raw, $fusion_gene_by_id, $read_id) = split /\|/, $id;
    my $sample_id;
    if ($sample_id_raw =~ /basecalled_reads.Q/) {
	($sample_id) = ($sample_id_raw =~ /(\S+)(?:\.basecalled_reads.Q.*.pass.fastq|\.fq)$/);
    } else {
	$sample_id = $sample_id_raw;
    }
    # print "line[0]=$line[0], sample_id = $sample_id\n";

    my @fusion_gene_by_id = split /:/, $fusion_gene_by_id;
    my @fusion_gene_by_id_tidy = ();
    my @fusion_gene_by_name = ();
    foreach my $gene_id (@fusion_gene_by_id) {
	$gene_id =~ s/^EN\.//gi;
	my $gene_name = $gene_id;
        if (exists $gene_id2name{$gene_id}) {
            $gene_name = $gene_id2name{$gene_id};
	}
	push @fusion_gene_by_id_tidy, $gene_id;
	push @fusion_gene_by_name, $gene_name;
    }
    my $fusion_gene_by_id_tidy = join ":", @fusion_gene_by_id_tidy;
    my $fusion_gene_by_name = join ":",@fusion_gene_by_name;
    print $output_fasta_fh ">$sample_id|$fusion_gene_by_name|$fusion_gene_by_id_tidy|$read_id\n";
    if ($u2t_conversion eq "yes") {
	$input_fasta{$id} =~ s/U/T/gi;
    }
    print $output_fasta_fh "$input_fasta{$id}\n";
}


sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    } else {
        open($fh, $file) or die "can't open $file";
    }
    return $fh;
}

sub write_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
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
