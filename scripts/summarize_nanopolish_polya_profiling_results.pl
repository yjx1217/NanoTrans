#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Descriptive::Discrete;
use Getopt::Long;

##############################################################                                                                                            #  script: summarize_nanopolish_polya_profiling_results.pl                                                                                        
#  author: Jia-Xing Yue (GitHub ID: yjx1217)                                                                                               
#  last edited: 2021.10.06                                                                                                                   
#  description: summarize nanopolish's polya results
#  example: perl summarize_nanopolish_polya_profiling_results.pl -i input.nanopolish.raw.txt -m transcript2gene_map.txt  -o output.nanopolish.filtered.txt
##############################################################                                                                               

my $input;
my $output;
my $transcript2gene_map;

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'transcript2gene_map|m:s' => \$transcript2gene_map);

my $input_fh = read_file($input);
my $output_fh = write_file($output);
my %data = ();
my @data = ();

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    /readname\tcontig/ and next;
    my ($readname, $contig, $position, $leader_start, $adapter_start, $polya_start, $transcript_start, $read_rate, $polya_length, $qc_tag) = split /\t/, $_;
    if ($qc_tag eq "PASS") {
	if (not exists $data{$contig}) {
	    @{$data{$contig}{'polya_length'}} = ($polya_length);
	} else {
	    push @{$data{$contig}{'polya_length'}}, $polya_length;
	}
    }
}

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


print $output_fh "isoform_id\tgene_id\tgene_name\tmedian\tmean\tmin\tmax\n";
foreach my $contig (sort keys %data) {
    my $polya_length_stat = new Statistics::Descriptive::Discrete;
    $polya_length_stat ->add_data(@{$data{$contig}{'polya_length'}});
    my $mean = $polya_length_stat->mean();
    my $median = $polya_length_stat->median();
    my $min = $polya_length_stat->min();
    my $max = $polya_length_stat->max();
    #my $transcript_id = $contig;
    #my $gene_id =  $transcript2gene_map{$transcript_id}{'gene_id'};
    my ($transcript_id, $gene_id) = split /_/, $contig;
    my $gene_name = $gene_id;
    if (exists $transcript2gene_map{$transcript_id}) {
	$gene_name = $transcript2gene_map{$transcript_id}{'gene_name'};
    } elsif (exists $gene_id2name{$gene_id}) {
	$gene_name = $gene_id2name{$gene_id};
    }
    print $output_fh "$transcript_id\t$gene_id\t$gene_name\t$median\t$mean\t$min\t$max\n";
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
    my ($fh, $input_hashref, $input_arrayref) = @_;
    my $seq_name = "";
    while (<$fh>) {
	chomp;
	if (/^\s*$/) {
	    next;
	} elsif (/^\s*#/) {
	    next;
	} elsif (/^>(.*)/) {
	    $seq_name = $1;
	    push @$input_arrayref, $seq_name;
	    $$input_hashref{$seq_name} = "";
	} else {
	    $$input_hashref{$seq_name} .= $_;
	}
    }
}

