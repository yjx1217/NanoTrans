#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_gene_fusion_detection.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.09
#  description: gene fusion detection by long reads
#  example: perl batch_gene_fusion_detection.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -long_reads_dir ./../00.Long_Reads  -x ./../00.Reference_Genome/ref.transcript2gene_map.txt
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $java_dir = $ENV{java_dir};
my $jaffal_dir = $ENV{jaffal_dir};
my $bpipe_dir = $ENV{bpipe_dir};
# my $ucsc_dir = $ENV{ucsc_dir};
# my $minimap2_dir = $ENV{minimap2_dir};
# my $samtools_dir = $ENV{samtools_dir};
# my $picard_dir = $ENV{picard_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $transcript2gene_map;
my $for_human = "no";
my $long_reads_dir = "./../00.Long_Reads";
my $ref_dir = "./../00.Reference_Genome";
my $min_mapping_quality = 1;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
           'for_human|h:s' => \$for_human,
           'ref_dir|ref_dir:s' => \$ref_dir,
	   'long_reads_dir|long_reads_dir:s' => \$long_reads_dir,
	   'min_mapping_quality|mq:i' => \$min_mapping_quality,
	   'transcript2gene_map|x:s' => \$transcript2gene_map,
	   'debug|d:s' => \$debug);

my $local_time = localtime();
print "\n\n[$local_time] Starting gene fusion detection for the batch $batch_id ..\n";

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
my $output_dir = "$batch_id";
system("mkdir $output_dir");

my $sample_count = 0;
foreach my $sample_id (@sample_table) {
    my $basecalled_fastq_file = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fastq_file'}";
    # my $basecalled_fast5_dir = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}";
    # my $basecalled_sequencing_summary = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}/sequencing_summary.txt";

    $local_time = localtime();
    print "\n[$local_time] Check the specified long read file ..:\n";
    print "Check the specified long read file:\n";
    if (-e $basecalled_fastq_file) {
        print "Successfully located the specified long read file: $basecalled_fastq_file.\n";
    } else {
        print "Cannot find the specified long read file: $basecalled_fastq_file!\n";
        print "Exit!\n";
        exit;
    }

    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    $local_time = localtime();
    print "\n[$local_time] Processing sample $sample_id for read mapping\n";
    my $local_time = localtime();
    print "\n[$local_time] processing sample $sample_id for clustered isoform based mapping ..\n";
    if ($for_human eq "yes") {
	system("$bpipe_dir/bpipe run --threads $threads -p refBase=$base_dir/$ref_dir -p genome=genome -p annotation=annotation -p knownTable=$NANOTRANS_HOME/misc/human.known_fusions.txt  -p jaffa_output=$sample_output_dir/ $jaffal_dir/JAFFAL.groovy $basecalled_fastq_file");
    } else {
	system("$bpipe_dir/bpipe run --threads $threads -p refBase=$base_dir/$ref_dir -p genome=genome -p annotation=annotation -p jaffa_output=$sample_output_dir/ $jaffal_dir/JAFFAL.groovy $basecalled_fastq_file");
    }
    system("perl $NANOTRANS_HOME/scripts/tidy_jaffal_output.pl -icsv jaffa_results.csv -ifa jaffa_results.fasta -p $sample_id -x $base_dir/$transcript2gene_map -u2t yes");
    $local_time = localtime();
    print "\n[$local_time] Finishing JAFFAL fusion gene detection for sample $sample_id \n";
    system("rm -r tmp");
    $sample_count++;
    chdir("./../") or die "cannot change directory to: $!\n";
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
    }
}
print $design_yml_fh "out: $combined_output_dir\n";

$local_time = localtime();

chdir("$combined_output_dir") or die "cannot change directory to: $!\n";
print "\n[$local_time] Copying final results into $combined_output_dir ..\n";
foreach my $sample_id (@sample_table) {
    system("cp ./../$sample_id/$sample_id.gene_fusion.transcripts.fa .");
    system("cp ./../$sample_id/$sample_id.gene_fusion.transcripts.txt .");
}

$local_time = localtime();
print "\n[$local_time] A total of $sample_count samples were processed!\n";


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
