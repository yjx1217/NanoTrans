#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: pool_polyA_length_distribution_summary.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.08
#  description: pool together polyA length distribution summary files
#  example: perl pool_polyA_length_distribution_summary.pl  -i Master_Sample_Table.txt(.gz) -o output.txt(.gz)
##############################################################

my ($sample_table, $output);
GetOptions('sample_table|i:s' => \$sample_table, # input table filename
           'output|o:s' => \$output); # output table filename

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);

my $output_fh = write_file($output);
print $output_fh "sample_id\tcomparison_group\treplicate_id\tisoform_id\tgene_id\tgene_name\tmedian_length\tmean_length\tmin_length\tmax_length\n";
my %comparison_groups = ();
foreach my $sample_id (@sample_table) {
    my $input = "./../$sample_id/$sample_id.nanopolish.polya_profiling.filtered.summary.txt";
    my $input_fh = read_file($input);
    while (<$input_fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/^isoform_id\tgene_id/ and next;
	print $output_fh "$sample_id\t$sample_table{$sample_id}{'comparison_group'}\t$sample_table{$sample_id}{'replicate_id'}\t$_\n";
    }
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
