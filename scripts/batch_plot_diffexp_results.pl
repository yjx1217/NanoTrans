#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_plot_diffexp_results.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.09
#  description: plot differential expression results
#  example: perl batch_plot_diffexp_results.pl -i Master_Sample_Table.txt -threads 1 -b $batch_id -lfc 1 -p 0.05 
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $log2foldchange_cutoff = 1;
my $adj_p_value_cutoff = 0.05;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch|b:s' => \$batch_id,
	   'log2foldchange_cutoff|lfc:s' => \$log2foldchange_cutoff,
	   'adj_p_value_cutoff|p:s' => \$adj_p_value_cutoff,
	   'debug|d:s' => \$debug);

my $local_time = localtime();
print "\n\n[$local_time] Starting differential expression plotting for batch $batch_id ..\n";

my $base_dir = cwd();
#print "base_dir=$base_dir\n";

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);


my $combined_output_dir = "$base_dir/$batch_id/all_samples_combined";
my %comparison_groups = ();
foreach my $sample_id (@sample_table) {
    my $g = $sample_table{$sample_id}{'comparison_group'};
    if (not exists $comparison_groups{$g}) {
        @{$comparison_groups{$g}} = ($sample_id);
    } else {
        push @{$comparison_groups{$g}}, $sample_id;
    }
}

$local_time = localtime();
print "\n[$local_time] Check numbers of comparison groups and in-group samples for plotting ..\n";

my @sample_comparison_pair = ();
my @comparison_groups = sort keys  %comparison_groups;
my $comparison_groups_count = scalar @comparison_groups;

my $min_replicate_count = 99999;
foreach my $g (@comparison_groups) {
    my $in_group_sample_count = scalar @{$comparison_groups{$g}};
    if ($in_group_sample_count < $min_replicate_count) {
	$min_replicate_count = $in_group_sample_count;
    }
}

$local_time = localtime();
# print "\n[$local_time] Detect $comparison_groups_count groups defined in the $sample_table.\n";
# print "\n[$local_time] The minimal number of in-group replicates is $min_replicate_count. \n";

chdir("$combined_output_dir") or die "cannot change directory to: $!\n";
my $between_group_comparison_count = 0;
if ($comparison_groups_count <= 1) {
    print "\n!!!!!!!!!!\n";
    print "Not enough comparison groups for making comparison! Only $comparison_groups_count group was defined in the $sample_table.\n";
    print "Exit\n";
    print "!!!!!!!!!!\n";
} else {
    $local_time = localtime();
    print "\n[$local_time] Perform between-group comparison plotting ..\n";
    for (my $i = 0; $i < $comparison_groups_count - 1; $i++) {
	for (my $j = 1; $j < $comparison_groups_count; $j++) {
	    my $comparison_group_pair = "$comparison_groups[$i]-$comparison_groups[$j]";
	    $local_time = localtime();
	    print "\n[$local_time] Making comparison between groups $comparison_groups[$i] and $comparison_groups[$j] ..\n";
	    $between_group_comparison_count++;
	    if ($min_replicate_count >= 2) {
		$local_time = localtime();
		print "\n[$local_time] Minimal replicate count = $min_replicate_count (>= 2). Perform between group comparison with replicate-based batch effect removal ..\n";
		$local_time = localtime();
                chdir("${batch_id}_differential_expression_output") or die "cannot change directory to: $!\n";
		print "\n[$local_time] Plotting differential isoform and gene expression results ..\n";

		system("Rscript --vanilla $NANOTRANS_HOME/scripts/plot_diffexp_results.R --input_data $batch_id.${comparison_group_pair}_comparison.differential_isoform_expression.tidy.txt --output_prefix $batch_id.${comparison_group_pair}_comparison.differential_isoform_expression.tidy.for_plot --comparison_pair_by_file yes --log2foldchange_cutoff $log2foldchange_cutoff --adj_p_value_cutoff $adj_p_value_cutoff");
		system("Rscript --vanilla $NANOTRANS_HOME/scripts/plot_diffexp_results.R --input_data $batch_id.${comparison_group_pair}_comparison.differential_gene_expression.tidy.txt --output_prefix $batch_id.${comparison_group_pair}_comparison.differential_gene_expression.tidy.for_plot --comparison_pair_by_file yes --log2foldchange_cutoff $log2foldchange_cutoff --adj_p_value_cutoff $adj_p_value_cutoff");
		$local_time = localtime();
	    } else {
		$local_time = localtime();
		print "\n[$local_time] Minimal replicate count = $min_replicate_count (< 3). Perform between group comparison without replicate-based batch effect removal ..\n";
		$local_time = localtime();
		$between_group_comparison_count = 0;
		print "\n[$local_time] Skipping plotting ..\n";
                # print "\n[$local_time] Plotting isoforms with differential expression and splicing respectively ..\n";
		
		# foreach my $sample_in_group_i (@{$comparison_groups{$comparison_groups[$i]}}) {
		#     foreach my $sample_in_group_j (@{$comparison_groups{$comparison_groups[$j]}}){
		# 	my $comparison_sample_pair = "${sample_in_group_i}-${sample_in_group_j}";
		# 	$local_time = localtime();
		# 	print "\n[$local_time] Plotting the comparison between samples ${sample_in_group_i} and ${sample_in_group_j} ..\n";
		# 	my $sample_in_group_i_tag = "${sample_in_group_i}_$sample_table{$sample_in_group_i}{'comparison_group'}_$sample_table{$sample_in_group_i}{'replicate_id'}";
		# 	my $sample_in_group_j_tag = "${sample_in_group_j}_$sample_table{$sample_in_group_j}{'comparison_group'}_$sample_table{$sample_in_group_j}{'replicate_id'}";
		# 	# system("Rscript --vanilla $NANOTRANS_HOME/scripts/plot_diffexp_results.R --input_data $batch_id.$comparison_sample_pair.differential_expressed_isoforms.statistically_tested.tidy.txt --output_prefix $batch_id.$comparison_sample_pair.differential_expressed_isoforms.statistically_tested.tidy --sample_a $sample_in_group_i_tag --sample_b $sample_in_group_j_tag --comparison_pair_by_file no --log2foldchange_cutoff $log2foldchange_cutoff --adj_p_value_cutoff $adj_p_value_cutoff");
		# 	# my @diff_splicing_types = qw(ir es alt3 alt5);
		# 	# foreach my $diff_splicing_type (@diff_splicing_types) {
		# 	#     system("perl $NANOTRANS_HOME/scripts/tidy_diffsplice_fishers_exact_output.pl -i $batch_id.$comparison_sample_pair.differential_spliced_isoforms.statistically_tested.$diff_splicing_type.events.quant.raw.txt -o $batch_id.$comparison_sample_pair.differential_spliced_isoforms.statistically_tested.$diff_splicing_type.events.quant.tidy.txt -x $base_dir/$transcript2gene_map");
		# 	# }
		#     }
                # }
	    }
	}
    }
}

$local_time = localtime();
print "\n[$local_time] Finishing differential expression plotting for batch $batch_id.\n";
print "\nA total of $between_group_comparison_count comparison between $comparison_groups_count groups were plotted in total.\n";
$local_time = localtime();
print "\n[$local_time] Done\n";

chdir("$base_dir") or die "cannot change directory to: $!\n";

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
