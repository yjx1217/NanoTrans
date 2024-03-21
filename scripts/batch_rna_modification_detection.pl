#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_rna_modification_detection.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.09.21
#  description: RNA modification detection
#  example: perl batch_rna_modification_detection.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -isoform_cq_dir ./../02.Isoform_Clustering_and_Quantification -long_reads_dir ./../00.Long_Reads
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $minimap2_dir = $ENV{minimap2_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $java_dir = $ENV{java_dir};
my $xpore_dir = $ENV{xpore_dir};
my $nanopolish_dir = $ENV{nanopolish_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $long_reads_dir = "./../00.Long_Reads";
my $isoform_cq_dir = "./../02.Isoform_Clustering_and_Quantification";
my $transcript2gene_map = "./../00.Reference_Genome/ref.transcript2gene_map.txt";
my $min_mapping_quality = 0;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
           'long_reads_dir|long_reads_dir:s' => \$long_reads_dir,
	   'isoform_cq_dir|isoform_clustering_and_quantification_dir:s' => \$isoform_cq_dir,
	   'transcript2gene_map|x:s' => \$transcript2gene_map,
           'min_mapping_quality|q:i' => \$min_mapping_quality,
	   'debug|d:s' => \$debug);


my $local_time = localtime();
print "\n[$local_time] Starting RNA modification detection for batch $batch_id ..\n";
my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
my $output_dir = "$batch_id";
system("mkdir $output_dir");

my $isoform_clustering_fasta_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.flair_all_collapsed.isoforms.fa";
my $isoform_clustering_gtf_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.flair_all_collapsed.isoforms.gtf";
my $isoform_quantification_tsv_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.counts_matrix.tsv";

if (-e $isoform_clustering_fasta_file) {
    print "Successfully located the isoform clustering fasta file: $isoform_clustering_fasta_file\n";
} else {
    print "Cannot find the isoform clustering fasta file: $isoform_clustering_fasta_file \n";
    print "Exit!\n";
    exit;
}

if (-e $isoform_clustering_gtf_file) {
    print "Successfully located the isoform clustering gtf file: $isoform_clustering_gtf_file\n";
} else {
    print "Cannot find the isoform clustering gtf file: $isoform_clustering_gtf_file \n";
    print "Exit!\n";
    exit;
}

if (-e $isoform_quantification_tsv_file) {
    print "Successfully located the isoform quantification tsv file: $isoform_quantification_tsv_file\n"; 
} else {
    print "Cannot find the isoform quantification tsv file: $isoform_quantification_tsv_file\n";
    print "Exit!\n";
    exit;
}


# long-read mapping to the isoform_clustering_file

my $sample_count = 0;
foreach my $sample_id (@sample_table) {
    $sample_count++;
    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    $local_time = localtime();
    print "\n[$local_time] Processing isoform-based mapping for sample $sample_id ..\n";

    my $basecalled_fastq_file = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fastq_file'}";
    my $basecalled_fast5_dir = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}";
    my $basecalled_sequencing_summary = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}/sequencing_summary.txt";
    print "Check the specified long read file:\n";
    if (-e $basecalled_fastq_file) {
        print "Successfully located the specified long read file: $basecalled_fastq_file.\n";
    } else {
        print "Cannot find the specified long read file: $basecalled_fastq_file!\n";
        print "Exit!\n";
        exit;
    }
    if (-e $basecalled_fast5_dir) {
        print "Successfully located the specified long read file: $basecalled_fast5_dir.\n";
    } else {
        print "Cannot find the specified long read file: $basecalled_fast5_dir!\n";
        print "Exit!\n";
        exit;
    }
    if (-e $basecalled_sequencing_summary) {
        print "Successfully located the specified long read file: $basecalled_sequencing_summary.\n";
    } else {
        print "Cannot find the specified long read file: $basecalled_sequencing_summary!\n";
        print "Exit!\n";
        exit;
    }

    my $local_time = localtime();
    print "\n[$local_time] processing sample $sample_id for clustered isoform based mapping ..\n";
    system("$minimap2_dir/minimap2 -t $threads -ax map-ont -p 0.99 -uf --secondary=no $isoform_clustering_fasta_file $basecalled_fastq_file | $samtools_dir/samtools view --threads $threads -F260 -hbS -T $isoform_clustering_fasta_file -q $min_mapping_quality - >$sample_id.bam");
    # -F 260  output primary aligned mapped reads
    # read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
    # bit 3 + bit 9 = 4 + 256 = 260
    
    ## sort bam file by picard-tools: SortSam
    $local_time = localtime();
    print "\n[$local_time] Sorting the mapping BAM file for sample $sample_id ..\n";
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam -INPUT $sample_id.bam -OUTPUT $sample_id.sort.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT -MAX_RECORDS_IN_RAM 50000");
    if ($debug eq "no") {
	system("rm $sample_id.bam");
    }
    # index final bam file
    $local_time = localtime();
    print "\n[$local_time] Indexing the mapping BAM file for sample $sample_id ..\n";
    system("$samtools_dir/samtools index $sample_id.sort.bam");
    
    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with nanopolish pre-processing ..\n";
    if ((not -e "$basecalled_fastq_file.index") or (not -e "$basecalled_fastq_file.index.gzi") or (not -e "$basecalled_fastq_file.index.fai") or (not -e "$basecalled_fastq_file.index.readdb")) {
	system("$nanopolish_dir/nanopolish index -s $basecalled_sequencing_summary -d $basecalled_fast5_dir $basecalled_fastq_file");
    }
    system("$nanopolish_dir/nanopolish eventalign --threads $threads --reads $basecalled_fastq_file --bam $sample_id.sort.bam --genome $isoform_clustering_fasta_file --signal-index --scale-events --summary $sample_id.isoform_based_mapping.eventalign.summary.txt > $sample_id.isoform_based_mapping.eventalign.txt");

    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with xpore dataprep ..\n";
    system("$xpore_dir/xpore dataprep --n_processes $threads --eventalign $sample_id.isoform_based_mapping.eventalign.txt --gtf_or_gff $isoform_clustering_gtf_file --transcript_fasta $isoform_clustering_fasta_file --out_dir dataprep");
    system("gzip $sample_id.isoform_based_mapping.eventalign.txt");
    chdir("./../") or die "cannot change directory to: $!\n";
}

$local_time = localtime();
print "\n[$local_time] Prepair yml file for experimental design ..\n";
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
	print $design_yml_fh "    ${rep_id}: $base_dir/$output_dir/$s/dataprep\n";
    }
}

my $combined_output_dir = "$base_dir/$output_dir/all_samples_combined";
system("mkdir -p $combined_output_dir");
print $design_yml_fh "out: $combined_output_dir\n";


$local_time = localtime();
print "\n[$local_time] Perform differential modification profiling ..\n";
system("$xpore_dir/xpore diffmod --n_processes $threads --config $batch_id.experimental_design.yml");
system("$xpore_dir/xpore postprocessing --diffmod_dir $combined_output_dir");

system("perl $NANOTRANS_HOME/scripts/tidy_xpore_output.pl -i $combined_output_dir/diffmod.table -o $combined_output_dir/$batch_id.rna_modification.diffmod.table.tidy.txt -x $base_dir/$transcript2gene_map");
system("perl $NANOTRANS_HOME/scripts/tidy_xpore_output.pl -i $combined_output_dir/majority_direction_kmer_diffmod.table -o $combined_output_dir/$batch_id.rna_modification.majority_direction_kmer_diffmod.table.tidy.txt -x $base_dir/$transcript2gene_map");

$local_time = localtime();
print "\n[$local_time] Clean up un-needed files for the batch $batch_id ..\n";

if ($debug eq "no") {
    system("rm $combined_output_dir/diffmod.table");
    system("rm $combined_output_dir/majority_direction_kmer_diffmod.table");
    system("rm $combined_output_dir/diffmod.log");
    system("rm -r $combined_output_dir/models");
}
$local_time = localtime();
print "\n[$local_time] Finishing the RNA modification identification for the batch $batch_id.\n";
$local_time = localtime();
print "\n[$local_time] Done!\n";

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
