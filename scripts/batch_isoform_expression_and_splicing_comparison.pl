#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_isoform_differential_expression_and_splicing_identification.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.09.21
#  description: isoform differential expression and splicing identification
#  example: perl batch_isoform_differential_expression_and_splicing_identification.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -cq_dir ./../
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $minimap2_dir = $ENV{minimap2_dir};
my $bedtools_dir = $ENV{bedtools_dir};
my $samtools_dir = $ENV{samtools_dir};
my $ucsc_dir = $ENV{ucsc_dir};
my $flair_dir = $ENV{flair_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $isoform_cq_dir = "./../02.Isoform_Clustering_and_Quantification";
my $transcript2gene_map;
my $read_counts_cutoff = 5;
my $log2foldchange_cutoff = 1;
my $adj_p_value_cutoff = 0.05;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'isoform_cq_dir|isoform_clustering_and_quantification_dir:s' => \$isoform_cq_dir,
	   'transcript2gene_map|x:s' => \$transcript2gene_map,
	   'read_counts_cutoff|e:i' => \$read_counts_cutoff,
	   'log2foldchange_cutoff|lfc:s' => \$log2foldchange_cutoff,
	   'adj_p_value_cutoff|p:s' => \$adj_p_value_cutoff,
	   'debug|d:s' => \$debug);


my $local_time = localtime();
print "\n\n[$local_time] Starting isoform differential expression and splicing identification for batch $batch_id ..\n";
my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
my $output_dir = "$batch_id";
system("mkdir $output_dir");
chdir("$output_dir") or die "cannot change directory to: $!\n";

my $isoform_clustering_bed_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.flair_all_collapsed.isoforms.bed";
my $isoform_quantification_tsv_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.counts_matrix.tsv";

if (-e $isoform_clustering_bed_file) {
    print "\nSuccessfully located the isoform clustering BED file: $isoform_clustering_bed_file\n";
} else {
    print "\nCannot find the isoform clustering BED file: $isoform_clustering_bed_file \n";
    print "Exit!\n";
    exit;
}

if (-e $isoform_quantification_tsv_file) {
    print "\nSuccessfully located the isoform quantification TSV file: $isoform_quantification_tsv_file\n"; 
} else {
    print "\nCannot find the isoform quantification TSV file: $isoform_quantification_tsv_file\n";
    print "Exit!\n";
    exit;
}


$local_time = localtime();
print "\n[$local_time] Prepair yml file for experimental design ..\n";

my $combined_output_dir = "$base_dir/$output_dir/all_samples_combined";
system("mkdir -p $combined_output_dir");

my $design_yml_file = "$batch_id.experimental_design.yml";
my $design_yml_fh = write_file($design_yml_file);
my %comparison_groups = ();
foreach my $sample_id (@sample_table) {
    my $g = $sample_table{$sample_id}{'comparison_group'};
    if (not exists $comparison_groups{$g}) {
        @{$comparison_groups{$g}} = ($sample_id);
    } else {
        push @{$comparison_groups{$g}}, $sample_id;
    }
}

print $design_yml_fh "data:\n";
foreach my $g (sort keys %comparison_groups) {
    print $design_yml_fh "  $g:\n";
    foreach my $s (@{$comparison_groups{$g}}) {
        my $rep_id = $sample_table{$s}{'replicate_id'};
        print $design_yml_fh "    ${rep_id}: $base_dir/$output_dir/$s\n";
        # print $design_yml_fh "    rep${rep_id}: $base_dir/$output_dir/$s\n";
    }
}
print $design_yml_fh "out: $combined_output_dir\n";

$local_time = localtime();
print "\n[$local_time] Check numbers of comparison groups and in-group samples ..\n";

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
print "\n[$local_time] Detect $comparison_groups_count groups defined in the $sample_table.\n";
print "\n[$local_time] The minimal number of in-group replicates is $min_replicate_count. \n";

chdir("$combined_output_dir") or die "cannot change directory to: $!\n";
my $between_group_comparison_count = 0;
if ($comparison_groups_count <= 1) {
    print "\n!!!!!!!!!!\n";
    print "Not enough comparison groups for making comparison! Only $comparison_groups_count group was defined in the $sample_table.\n";
    print "Exit\n";
    print "!!!!!!!!!!\n";
} else {
    $local_time = localtime();
    print "\n[$local_time] Perform between group comparison ..\n";
    for (my $i = 0; $i < $comparison_groups_count - 1; $i++) {
	for (my $j = 1; $j < $comparison_groups_count; $j++) {
	    my $comparison_group_pair = "$comparison_groups[$i]-$comparison_groups[$j]";
	    $local_time = localtime();
	    print "\n[$local_time] Making comparison between groups: $comparison_groups[$i] and $comparison_groups[$j] ..\n";
	    $between_group_comparison_count++;
	    if ($min_replicate_count >= 3) {
		$local_time = localtime();
		print "\n[$local_time] Minimal replicate count = $min_replicate_count (>= 3). Perform between group comparison with replicate-based batch effect removal ..\n";
		$local_time = localtime();
		print "\n[$local_time] Detecting isoforms with differential expression ..\n";
		system("$flair_dir/flair diffExp -t $threads -q $isoform_quantification_tsv_file -e ${read_counts_cutoff} -o ${batch_id}_differential_expression_output");
		chdir("${batch_id}_differential_expression_output") or die "cannot change directory to: $!\n";
		system("mv workdir/dge_stderr.txt workdir/log.txt");
		# generate tidy result tables
		my $comparison_group_pair_tag;
		$comparison_group_pair_tag = $comparison_groups[$i] ."_v_". $comparison_groups[$j];
		if (! -e "die_${comparison_group_pair_tag}_deseq2_results_shrinkage.tsv") {
		    $comparison_group_pair_tag = $comparison_groups[$j] ."_v_". $comparison_groups[$i];
		}
		# print "comparison_group_pair_tag=$comparison_group_pair_tag\n";
		system("$NANOTRANS_HOME/scripts/tidy_diffExp_die_deseq2_output.pl  -sample_table $base_dir/$sample_table -i workdir/isoforms_deseq2_${comparison_group_pair_tag}_results_shrinkage.tsv -o ${batch_id}.${comparison_group_pair}_comparison.differential_isoform_expression.tidy.txt  -m workdir/filtered_iso_counts_ds2.tsv  -x $base_dir/$transcript2gene_map");
		system("$NANOTRANS_HOME/scripts/tidy_diffExp_dge_deseq2_output.pl  -sample_table $base_dir/$sample_table -i workdir/genes_deseq2_${comparison_group_pair_tag}_results_shrinkage.tsv    -o ${batch_id}.${comparison_group_pair}_comparison.differential_gene_expression.tidy.txt     -m workdir/filtered_gene_counts_ds2.tsv -x $base_dir/$transcript2gene_map");
		system("$NANOTRANS_HOME/scripts/tidy_diffExp_diu_drimseq_output.pl -sample_table $base_dir/$sample_table -i workdir/isoforms_drimseq_${comparison_group_pair_tag}_results.tsv          -o ${batch_id}.${comparison_group_pair}_comparison.differential_isoform_usage.tidy.txt       -m workdir/filtered_iso_counts_drim.tsv -x $base_dir/$transcript2gene_map");
		#system("mkdir intermediate_files");
		#system("mv filtered_iso_counts_ds2.tsv intermediate_files");
		#system("mv filtered_gene_counts_ds2.tsv intermediate_files");
		#system("mv filtered_iso_counts_drim.tsv intermediate_files");
		#system("mv formula_matrix.tsv intermediate_files");
		#system("mv dge_QCplots_${comparison_group_pair_tag}.pdf intermediate_files");
		#system("mv die_QCplots_${comparison_group_pair_tag}.pdf intermediate_files");
		#system("mv log.txt intermediate_files");
		#system("mv die_${comparison_group_pair_tag}_deseq2_results_shrinkage.tsv intermediate_files");
		#system("mv dge_${comparison_group_pair_tag}_deseq2_results_shrinkage.tsv intermediate_files");
		#system("mv die_${comparison_group_pair_tag}_deseq2_results.tsv intermediate_files");
		#system("mv dge_${comparison_group_pair_tag}_deseq2_results.tsv intermediate_files");
		#system("mv diu_${comparison_group_pair_tag}_drimseq2_results.tsv intermediate_files");
		chdir("./../") or die "cannot change directory to: $!\n";
    
		$local_time = localtime();
		print "\n[$local_time] Detecting isoforms with differential isoform splicing ..\n";
		system("mkdir ${batch_id}_differential_splicing_output");
		system("$flair_dir/flair diffSplice -t $threads -i $isoform_clustering_bed_file -q $isoform_quantification_tsv_file --test --batch --conditionA $comparison_groups[$i] --conditionB $comparison_groups[$j] --drim1 6 --drim2 3 --drim3 15 --drim4 5 -o $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing");
		#system("mv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.* ${batch_id}_differential_splicing_output");
		chdir("${batch_id}_differential_splicing_output") or die "cannot change directory to: $!\n";

                $comparison_group_pair_tag = $comparison_groups[$i] ."_v_". $comparison_groups[$j];
                if (! -e "./../$batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing/workdir/es_${comparison_group_pair_tag}_drimseq_results.tsv") {
                    $comparison_group_pair_tag = $comparison_groups[$j] ."_v_". $comparison_groups[$i];
                }

		my @diff_splicing_types = qw(ir es alt3 alt5);
		foreach my $diff_splicing_type (@diff_splicing_types) {
		    system("$NANOTRANS_HOME/scripts/tidy_diffsplice_drimseq_output.pl -i ./../$batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing/workdir/${diff_splicing_type}_${comparison_group_pair_tag}_drimseq_results.tsv -a $comparison_groups[$i] -b $comparison_groups[$j] -q ./../$batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing/diffsplice.${diff_splicing_type}.events.quant.tsv -o $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.tidy.txt -x $base_dir/$transcript2gene_map");
		}
		#system("mkdir intermediate_files");
		#system("mv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.stderr.txt log.txt");
		#system("mv log.txt intermediate_files");
		#system("mv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.*.events.quant.tsv intermediate_files");
		#system("mv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.*drimseq_results.tsv intermediate_files");
		chdir("./../") or die "cannot change directory to: $!\n";
	    } else {
		$local_time = localtime();
		print "\n[$local_time] Minimal replicate count = $min_replicate_count (< 3). Perform between group comparison without replicate-based batch effect removal ..\n";
		$local_time = localtime();
		print "\n[$local_time] Detecting isoforms with differential expression and splicing respectively ..\n";
		
		system("mkdir ${batch_id}_differential_expression_output");
		system("mkdir ${batch_id}_differential_splicing_output");
		system("mkdir -p ${batch_id}_differential_expression_output/intermediate_files");
		system("mkdir -p ${batch_id}_differential_splicing_output/intermediate_files");


		system("$flair_dir/flair diffSplice -t $threads -i $isoform_clustering_bed_file -q $isoform_quantification_tsv_file -o $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw");
		system("cp $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw.es.events.quant.tsv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw.es.events.quant.tsv.bk");
		system("perl $NANOTRANS_HOME/scripts/fix_flair_diffSpice_id_parsing_issue.pl -i $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw.es.events.quant.tsv.bk -o $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw.es.events.quant.tsv");

		foreach my $sample_in_group_i (@{$comparison_groups{$comparison_groups[$i]}}) {
		    foreach my $sample_in_group_j (@{$comparison_groups{$comparison_groups[$j]}}){
			my $comparison_sample_pair = "${sample_in_group_i}-${sample_in_group_j}";
			$local_time = localtime();
			print "\n[$local_time] Making comparison between samples: ${sample_in_group_i} and ${sample_in_group_j} ..\n";
			my $sample_in_group_i_tag = "${sample_in_group_i}_$sample_table{$sample_in_group_i}{'comparison_group'}_$sample_table{$sample_in_group_i}{'replicate_id'}";
			my $sample_in_group_j_tag = "${sample_in_group_j}_$sample_table{$sample_in_group_j}{'comparison_group'}_$sample_table{$sample_in_group_j}{'replicate_id'}";
			system("$flair_dir/diff_iso_usage $isoform_quantification_tsv_file $sample_in_group_i_tag $sample_in_group_j_tag $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.raw.txt");
			system("perl $NANOTRANS_HOME/scripts/tidy_diff_iso_usage_output.pl -i $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.raw.txt -o $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.tidy.txt.tmp -a $sample_in_group_i_tag -b $sample_in_group_j_tag -x $base_dir/$transcript2gene_map");
			system("Rscript --vanilla $NANOTRANS_HOME/scripts/apply_p_value_fdr_adjustment.R $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.tidy.txt.tmp $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.tidy.txt");
			system("rm $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.tidy.txt.tmp");
			system("mv $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.raw.txt ${batch_id}_differential_expression_output/intermediate_files");
			system("mv $batch_id.${comparison_sample_pair}_comparison.differential_isoform_usage.tidy.txt ${batch_id}_differential_expression_output");
			my @diff_splicing_types = qw(ir es alt3 alt5);
			foreach my $diff_splicing_type (@diff_splicing_types) {
			    system("$flair_dir/diffsplice_fishers_exact $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw.$diff_splicing_type.events.quant.tsv  $sample_in_group_i_tag $sample_in_group_j_tag $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.raw.txt");
			    system("perl $NANOTRANS_HOME/scripts/tidy_diffsplice_fishers_exact_output.pl -i $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.raw.txt -o $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.tidy.txt.tmp -x $base_dir/$transcript2gene_map");
			    system("Rscript --vanilla $NANOTRANS_HOME/scripts/apply_p_value_fdr_adjustment.R $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.tidy.txt.tmp $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.tidy.txt");
			    system("rm $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.tidy.txt.tmp");
			    system("mv $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.raw.txt ${batch_id}_differential_splicing_output/intermediate_files");
			    system("mv $batch_id.${comparison_sample_pair}_comparison.differential_isoform_splicing.$diff_splicing_type.events.quant.tidy.txt ${batch_id}_differential_splicing_output");
			}
		    }
		}
		system("mv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw*.tsv ${batch_id}_differential_splicing_output/intermediate_files");
		system("mv $batch_id.${comparison_group_pair}_comparison.differential_isoform_splicing.raw*.tsv.bk ${batch_id}_differential_splicing_output/intermediate_files");
	    }
	}
    }
}

$local_time = localtime();
print "\n[$local_time] Finishing the isoform differential expression and splicing identification for batch $batch_id.\n";
print "\nA total of $between_group_comparison_count comparison between $comparison_groups_count groups were made in total.\n";
$local_time = localtime();
print "\n[$local_time] Done\n";

# chdir("$base_dir") or die "cannot change directory to: $!\n";

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
