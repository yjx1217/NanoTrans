#! /usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use List::Util qw(first);

##############################################################
#  script: tidy_count_matrix_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.30
#  description: tidy up flair's count matrix output file
#  example: perl tidy_count_matrix_output.pl -s Master_Sample_Table.txt -i input.txt(.gz) -o output.txt(.gz) -x transcript2gene_map
##############################################################

my ($sample_table, $input, $output, $filtered_count_matrix, $transcript2gene_map);
GetOptions('sample_table|sample_table:s' => \$sample_table,
	   'input|i:s' => \$input, 
	   'output|o:s' => \$output,
           'transcript2gene_map|x:s' => \$transcript2gene_map);

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my %comparison_groups = ();
foreach my $sample_id (@sample_table) {
    my $g = $sample_table{$sample_id}{'comparison_group'};
    if (not exists $comparison_groups{$g}) {
        @{$comparison_groups{$g}} = ($sample_id);
    } else {
        push @{$comparison_groups{$g}}, $sample_id;
    }
}

my $sample_count = scalar @sample_table;


# print "sample_tags: @sample_tags\n";

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


my $input_fh = read_file($input);
my %input = ();
parse_input_file($input_fh, \%input);

my $output_fh = write_file($output);
print $output_fh "isoform_id\tgene_id\tgene_name";
foreach my $g (sort keys %comparison_groups) {
    foreach my $s (@{$comparison_groups{$g}}) {
	print $output_fh "\t$s";
    }
}
print $output_fh "\n";

foreach my $contig_id (sort keys %input) {
    my ($isoform_id, $gene_id) = split "_", $contig_id;
    my $gene_name = $gene_id;
    if (exists $gene_id2name{$gene_id}) {
	$gene_name = $gene_id2name{$gene_id};
    }
    print $output_fh "$isoform_id\t$gene_id\t$gene_name";
    foreach my $g (sort keys %comparison_groups) {
	foreach my $s (@{$comparison_groups{$g}}) {
	    print $output_fh "\t$input{$contig_id}{$s}";
	}
    }
    print $output_fh "\n";
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

sub parse_input_file {
    my ($fh, $input_hashref) = @_;
    my @raw_sample_tags = ();
    my $sample_count;
    while (<$fh>) {
	chomp;
	/^\s*$/ and next;
        /^\s*#/ and next;
	if ($. == 1) {
	    my @header = split "\t", $_;
	    $sample_count = (scalar @header) - 1;
	    @raw_sample_tags = splice @header, 1, $sample_count;
	} 
	else {
	    my @count = split "\t", $_;
	    my $contig = $count[0];
	    for (my $i = 0; $i < $sample_count; $i++) {
		my ($sample_id, $comparison_group, $replicate_id) = split /\_/, $raw_sample_tags[$i]; 
		$$input_hashref{$contig}{$sample_id} = $count[$i+1];
	    }
	}
    }
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
