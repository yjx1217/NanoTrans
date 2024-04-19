#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_polya_tail_length_profiling.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2024.04.17
#  description: profiling polyA tail length for gene and transcripts
#  example: perl batch_polya_tail_length_profiling.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -isoform_cq_dir ./../02.Isoform_Clustering_and_Quantification -long_reads_dir ./../00.Long_Reads -method nanopolish -x ./../00.Reference_Genome/ref.transcript2gene_map.txt
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $java_dir = $ENV{java_dir};
my $minimap2_dir = $ENV{minimap2_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $nanopolish_dir = $ENV{nanopolish_dir};
# my $tailfindr_dir = $ENV{tailfindr_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $transcript2gene_map;
my $isoform_cq_dir = "./../02.Isoform_Clustering_and_Quantification";
my $long_reads_dir = "./../00.Long_Reads";
my $min_mapping_quality = 1;
my $method = "nanopolish"; 
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'long_reads_dir|long_reads_dir:s' => \$long_reads_dir,
	   'isoform_clustering_and_quantification_dir|isoform_cq_dir:s' => \$isoform_cq_dir,
	   'min_mapping_quality|mq:i' => \$min_mapping_quality,
	   'method|m:s' => \$method,
	   'transcript2gene_map|x:s' => \$transcript2gene_map,
	   'debug|d:s' => \$debug);

my $local_time = localtime();
print "\n\n[$local_time] Starting polyA length profiling for batch $batch_id ..\n";

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
    my $basecalled_fast5_dir = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}";
    my $basecalled_sequencing_summary = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fast5_dir'}/sequencing_summary.txt";

    my $isoform_clustering_fasta_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.flair_all_collapsed.isoforms.fa";
    my $isoform_clustering_gtf_file = "$base_dir/$isoform_cq_dir/$batch_id/all_samples_combined/$batch_id.all_samples_combined.flair_all_collapsed.isoforms.gtf";

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
    if (-e $basecalled_fast5_dir) {
        print "Successfully located the specified long read file: $basecalled_fast5_dir.\n";
        if (-e $basecalled_sequencing_summary) {
            print "Successfully located the specified long read file: $basecalled_sequencing_summary.\n";
        } else {
            print "Cannot find the specified long read file: $basecalled_sequencing_summary!\n";
            print "INFO: Likely Dorado was chosen for basecalling, please adjust the basecalled_fast5_dir column in Master_Sample_table to raw_fast5_dir; \n";
            print "INFO: or double-check it is correct for $basecalled_fast5_dir!\n";
        }
    } else {
        print "Cannot find the specified long read file: $basecalled_fast5_dir!\n";
        print "Exit!\n";
        exit;
    }


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

    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    $local_time = localtime();
    print "\n[$local_time] Processing sample $sample_id for read mapping\n";
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
    print "\n[$local_time] Processing sample $sample_id for polyA tail profiling\n";
    # polyA profiling  
    if ($method eq "nanopolish") {
	if ((not -e "$basecalled_fastq_file.index") or (not -e "$basecalled_fastq_file.index.gzi") or (not -e "$basecalled_fastq_file.index.fai") or (not -e "$basecalled_fastq_file.index.readdb")) {
	    if (( not -e $basecalled_sequencing_summary)) {
            system("$nanopolish_dir/nanopolish index -d $basecalled_fast5_dir $basecalled_fastq_file"); 
        } else {
            system("$nanopolish_dir/nanopolish index -s $basecalled_sequencing_summary -d $basecalled_fast5_dir $basecalled_fastq_file"); 
        }

	}
	system("$nanopolish_dir/nanopolish polya --threads $threads --reads $basecalled_fastq_file --bam $sample_id.sort.bam --genome $isoform_clustering_fasta_file  > $sample_id.nanopolish.polya_profiling.raw.txt");
	system("head -1  $sample_id.nanopolish.polya_profiling.raw.txt > $sample_id.nanopolish.polya_profiling.header");
	system("cat $sample_id.nanopolish.polya_profiling.raw.txt | egrep  \"PASS\" > $sample_id.nanopolish.polya_profiling.filtered.content");
	system("cat $sample_id.nanopolish.polya_profiling.header $sample_id.nanopolish.polya_profiling.filtered.content > $sample_id.nanopolish.polya_profiling.filtered.txt");
	system("rm $sample_id.nanopolish.polya_profiling.header $sample_id.nanopolish.polya_profiling.filtered.content");
	system("perl $NANOTRANS_HOME/scripts/summarize_nanopolish_polya_profiling_results.pl -i $sample_id.nanopolish.polya_profiling.filtered.txt -m ./../../$transcript2gene_map -o $sample_id.nanopolish.polya_profiling.filtered.summary.txt");
    # } elsif ($method eq "tailfindr") {
    #   Note: tailfindr require turning off adaptor trimming during basecalling, which makes it incovenient to packed together with the nanopolish polyA analysis protocol.
    # 	my $run_tailfindr_find_tails_fh = write_file("$sample_id.run_tailfindr.find_tails.R");
    # 	print $run_tailfindr_find_tails_fh "#!/usr/bin/env Rscript\nlibrary(tailfindr); find_tails(fast5_dir = \"$basecalled_fast5_dir\" , save_dir = \".\", csv_filename = \"$sample_id.tailfindr.polya_profiling.find_tails.raw.csv\", num_cores = $threads);";
    # 	close $run_tailfindr_find_tails_fh;
    # 	system("$tailfindr_dir/R --slave --no-save < $sample_id.run_tailfindr.find_tails.R");
    # 	system("$samtools_dir/samtools view -h $sample_id.filtered.bam  > $sample_id.filtered.sam");
    # 	my $run_tailfindr_annotate_tails_fh = write_file("$sample_id.run_tailfindr.annotate_tails.R");
    # 	print $run_tailfindr_annotate_tails_fh "#!/usr/bin/env Rscript\nlibrary(tailfindr); annotate_tails(\"$sample_output_dir/$sample_id.filtered.sam\", \"$sample_id.tailfindr.polya_profiling.find_tails.raw.csv\", \"$sample_id.tailfindr.polya_profiling.annotate_tails.raw.csv\")\n";
    # 	close $run_tailfindr_annotate_tails_fh;
    # 	system("$tailfindr_dir/R --slave --no-save < $sample_id.run_tailfindr.annotate_tails.R");
    # 	system("perl $NANOTRANS_HOME/scripts/summarize_tailfindr_polya_profiling_results.pl -i  $sample_id.tailfindr.polya_profiling.annotate_tails.raw.csv -o  $sample_id.tailfindr.polya_profiling.filtered.summary.txt -m ./../../$transcript2gene_map -q $min_mapping_quality");
    } else {
	# die "Unknown method: $method! Please use \"nanopolish\" or \"tailfindr\"!\n";
	die "Unknown method: $method! Please use \"nanopolish\"!\n";
    }
    
    $local_time = localtime();
    print "\n[$local_time] Finishing Nanopolish polyA profiling for sample $sample_id \n";
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
        print $design_yml_fh "    ${rep_id}: $base_dir/$output_dir/$s/dataprep\n";
    }
}
print $design_yml_fh "out: $combined_output_dir\n";

$local_time = localtime();



chdir("$combined_output_dir") or die "cannot change directory to: $!\n";

print "\n[$local_time] Examing polyA length difference between compared groups ..\n";
system("perl $NANOTRANS_HOME/scripts/pool_polyA_length_distribution_summary.pl -i $base_dir/$sample_table -o $batch_id.all_samples_combined.polya_profiling.summary.txt");
system("Rscript --vanilla --slave $NANOTRANS_HOME/scripts/plot_polyA_length_distribution.R --input $batch_id.all_samples_combined.polya_profiling.summary.txt --prefix $batch_id.all_samples_combined.polya_profiling.comparison");

$local_time = localtime();
print "\n[$local_time] A total of $sample_count samples were processed!\n";
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
