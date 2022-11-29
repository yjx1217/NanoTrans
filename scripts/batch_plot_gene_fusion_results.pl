#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_plot_gene_fusion_results.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.11.02
#  description: plot gene fusion results
#  example: perl batch_plot_gene_fusion_results.pl  -i Master_Sample_Table.txt -b $batch_id -g $gtf 
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $gtf;
my $threads = 1;
my $debug = "no";
GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch|b:s' => \$batch_id,
	   'gtf|g:s' => \$gtf,
	   'debug|d:s' => \$debug);

my $local_time = localtime();
print "\n\n[$local_time] Starting gene fusion plotting for batch $batch_id ..\n";

my $base_dir = cwd();
#print "base_dir=$base_dir\n";

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);


my $combined_output_dir = "$base_dir/$batch_id/all_samples_combined";
chdir("$combined_output_dir") or die "cannot change directory to: $combined_output_dir!\n";

system("mkdir gene_fusion_plots");
foreach my $sample_id (@sample_table) {
    chdir("$combined_output_dir/gene_fusion_plots") or die "cannot change directory to: $combined_output_dir/gene_fusion_plots!\n";
    system("Rscript --vanilla $NANOTRANS_HOME/scripts/plot_gene_fusion_results.R --input $combined_output_dir/$sample_id.gene_fusion.transcripts.txt --gtf $base_dir/$gtf --prefix $batch_id.$sample_id");
    if (-e "Rplots.pdf") {
	system("rm Rplots.pdf");
    }
    chdir("./../");
}

$local_time = localtime();
print "\n[$local_time] Finishing fusion gene plotting for batch $batch_id.\n";
$local_time = localtime();
print "\n[$local_time] Done\n";

chdir("$base_dir") or die "cannot change directory to: $base_dir\n";

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
