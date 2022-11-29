#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_isoform_clustering_and_quantification.pl 
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.11.14
#  description: isoform clustering and quantification
#  example: perl batch_flair_analysis.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -ref_dir ./../00.Reference_Genome  -mapping_dir ./../01.Reference_based_Mapping -long_reads_dir ./../00.Long_Reads -method nanopolish -x ./../00.Reference_Genome/ref.transcript2gene_map.txt
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $minimap2_dir = $ENV{minimap2_dir};
my $bedtools_dir = $ENV{bedtools_dir};
my $samtools_dir = $ENV{samtools_dir};
my $flair_dir = $ENV{flair_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $transcript2gene_map;
my $ref_dir = "./../00.Reference_Genome";
my $long_reads_dir = "./../00.Long_Reads";
my $mapping_dir = "./../01.Reference_Genome_based_Read_Mapping";
my $min_mapping_quality = 1;
my $method = "nanopolish"; 
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_dir|ref_dir:s' => \$ref_dir,
	   'long_reads_dir|long_reads_dir:s' => \$long_reads_dir,
	   'mapping_dir|read_mapping_dir:s' => \$mapping_dir,
	   'min_mapping_quality|mq:i' => \$min_mapping_quality,
	   'method|m:s' => \$method,
	   'transcript2gene_map|x:s' => \$transcript2gene_map,
	   'debug|d:s' => \$debug);


my $local_time = localtime();
print "\n\n[$local_time] Starting isoform clustering and quantification for batch $batch_id ..\n";
my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
my $output_dir = "$batch_id";
system("mkdir $output_dir");

my $sample_count = 0;
my @basecalled_fastq_file_array = ();
my @basecalled_fast5_dir_array = ();
my @aln_correction_output_array = ();

my $ref_genome_file = "$base_dir/$ref_dir/ref.genome.fa";
my $ref_genome_annotation_file = "$base_dir/$ref_dir/ref.genome.gtf";
my $ref_transcriptome_file = "$base_dir/$ref_dir/ref.transcriptome.fa";

foreach my $sample_id (@sample_table) {
    my $basecalled_fastq_file = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fastq_file'}";
    my $basecalled_fast5_dir = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}";;
    $sample_count++;
    my $bam_file = "$base_dir/$mapping_dir/$batch_id/$sample_id/$sample_id.sort.bam";
    print "Check the specified reference genome file:\n";
    if (-e $ref_genome_file) {
        print "Successfully located the specified reference genome file: $ref_genome_file.\n";
    } else {
        print "Cannot find the specified reference genome file: $ref_genome_file!\n";
        print "Exit!\n";
        exit;
    }
    print "Check the specified reference transcriptome file:\n";
    if (-e $ref_transcriptome_file) {
        print "Successfully located the specified reference transcriptome file: $ref_transcriptome_file.\n";
    } else {
        print "Cannot find the specified reference transcriptome file: $ref_transcriptome_file!\n";
        print "Exit!\n";
        exit;
    }
    print "Check the specified long-read mapping bam file:\n";
    if (-e $bam_file) {
        print "Successfully located the specified long-read mapping bam file: $bam_file.\n";
    } else {
        print "Cannot find the specified long-read mapping bam file: $bam_file!\n";
        print "Exit!\n";
        exit;
    }
    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    $local_time = localtime();
    print "\n[$local_time] Processing BAM correction for sample $sample_id ..\n";
    # alignment correction 
    if (-e $ref_genome_file) {
	system("$flair_dir/bam2Bed12 -i $bam_file > $sample_id.bam2bed.bed12");
	system("$flair_dir/flair correct -t $threads -q $sample_id.bam2bed.bed12 -g $ref_genome_file --nvrna -f $base_dir/$ref_dir/ref.genome.gtf -o $sample_id.flair");
	push @basecalled_fastq_file_array, $basecalled_fastq_file;
	push @basecalled_fast5_dir_array, $basecalled_fast5_dir;
	push @aln_correction_output_array, "$sample_output_dir/$sample_id.flair_all_corrected.bed";
    } else {
	print "Error! Cannot find $ref_genome_file or $ref_transcriptome_file\n";
    }
    system("rm -r tmp");
    chdir("./../") or die "cannot change directory to: $!\n";
}

# isoform collapsing
$local_time = localtime();
print "\n[$local_time] Combining all samples from the batch $batch_id for isoform clustering ..\n";

my $combined_output_dir = "$base_dir/$output_dir/all_samples_combined";
system("mkdir -p $combined_output_dir");
chdir("$combined_output_dir") or die "cannot change directory to: $!\n";
system("mkdir tmp");
# combine alignment correction results
my $basecalled_fastq_file_array_list = join " ", @basecalled_fastq_file_array;
my $basecalled_fast5_dir_array_list = join " ", @basecalled_fast5_dir_array;
my $aln_correction_output_array_list = join " ", @aln_correction_output_array;

system("cat $aln_correction_output_array_list > $batch_id.all_samples_combined.flair_all_corrected.bed");
# generate reads_manifest.tsv 
my $manifest_file = "$batch_id.all_samples_combined.reads_manifest.tsv";
my $manifest_fh = write_file($manifest_file);
foreach my $sample_id (@sample_table) {
    # my $replicate_batch_id = "batch_" . $sample_table{$sample_id}{'replicate_id'};
    print $manifest_fh "$sample_id\t$sample_table{$sample_id}{'comparison_group'}\t$sample_table{$sample_id}{'replicate_id'}\t$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fastq_file'}\n";
}

if ($debug eq "no") {
    system("$flair_dir/flair collapse -t $threads -q $batch_id.all_samples_combined.flair_all_corrected.bed -g $ref_genome_file -f $base_dir/$ref_dir/ref.genome.gtf -m $minimap2_dir/minimap2 -b $bedtools_dir/bedtools -sam $samtools_dir/samtools -r $basecalled_fastq_file_array_list --isoformtss --no_redundant none --generate_map --annotation_reliant generate --quality $min_mapping_quality -o $batch_id.all_samples_combined.flair_all_collapsed --temp_dir ./tmp");
} else {
    system("$flair_dir/flair collapse -t $threads -q $batch_id.all_samples_combined.flair_all_corrected.bed -g $ref_genome_file -f $base_dir/$ref_dir/ref.genome.gtf -m $minimap2_dir/minimap2 -b $bedtools_dir/bedtools -sam $samtools_dir/samtools -r $basecalled_fastq_file_array_list --isoformtss --no_redundant none --generate_map --annotation_reliant generate --quality $min_mapping_quality -o $batch_id.all_samples_combined.flair_all_collapsed --keep_intermediate --temp_dir ./tmp");
}

print "\n[$local_time] processing all samples from the Batch $batch_id for isoform quantification ..\n";
system("$flair_dir/flair quantify -m $minimap2_dir/minimap2 -sam $samtools_dir/samtools -t $threads --quality $min_mapping_quality -r $manifest_file -i $batch_id.all_samples_combined.flair_all_collapsed.isoforms.fa -o $batch_id.all_samples_combined.counts_matrix.tsv"); 

print "\n[$local_time] processing all samples from the Batch $batch_id for productivity prediction ..\n";
system("$flair_dir/predictProductivity --append_column --longestORF -i $batch_id.all_samples_combined.flair_all_collapsed.isoforms.bed -g $ref_genome_annotation_file -f $ref_genome_file > $batch_id.all_samples_combined.flair_all_collapsed.isoforms.with_productivity.bed");

system("$NANOTRANS_HOME/scripts/tidy_count_matrix_output.pl -s $base_dir/$sample_table -i $batch_id.all_samples_combined.counts_matrix.tsv -o $batch_id.all_samples_combined.counts_matrix.tidy.txt -x $base_dir/$transcript2gene_map");

system("rm -r tmp");
system("mkdir intermediate_files");
system("mv $batch_id.all_samples_combined.flair_all_collapsed.annotated_transcripts.* intermediate_files");
system("mv $batch_id.all_samples_combined.flair_all_collapsed.combined.isoform.read.map.txt intermediate_files");
system("mv $batch_id.all_samples_combined.flair_all_collapsed.isoform.read.map.txt intermediate_files");
system("mv $batch_id.all_samples_combined.flair_all_corrected.bed intermediate_files");
system("mv $batch_id.all_samples_combined.reads_manifest.tsv intermediate_files");

$local_time = localtime();
print "\n[$local_time] finishing the isoform clustering and quantification for batch $batch_id.\n";
print "A total of $sample_count samples were processed!\n";
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
