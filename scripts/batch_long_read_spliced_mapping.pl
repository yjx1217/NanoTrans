#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_long_read_spliced_mapping.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.09.26
#  description: run long-read mapping for a batch of samples
#  example: perl batch_long_read_mapping_spliced_mapping.pl -i Master_Sample_Table.txt -b $batch_id -threads 4  -ref_dir ./../00.Reference_Genome -long_reads_dir ./../00.Long_Reads -long_read_technology nanopore_drs -long_read_mapper minimap2-mapping_mode map2genome
##############################################################

my $NANOTRANS_HOME = $ENV{NANOTRANS_HOME};
my $java_dir = $ENV{java_dir};
my $minimap2_dir = $ENV{minimap2_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $sample_table = "Master_Sample_Table.txt";
my $threads = 1;
my $batch_id = "Batch_TEST";
my $long_read_technology = "nanopore_drs"; # "nanopore_drs" or "pacbio_isoseq"
my $long_read_mapper = "minimap2";
my $mapping_mode = "map2genome"; # or "map2genome";
my $min_mapping_quality = 0;
my $ref_dir = "./../00.Reference_Genome";
my $long_reads_dir = "./../00.Long_Reads";
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_dir|ref_dir:s' => \$ref_dir,
	   'mapping_mode|mapping_mode:s' => \$mapping_mode,
	   'long_reads_dir|long_reads_dir:s' => \$long_reads_dir,
	   'long_read_technology|long_read_technology:s' => \$long_read_technology,
	   'long_read_mapper|m:s' => \$long_read_mapper,
	   'min_mapping_quality|q:i' => \$min_mapping_quality,
	   'debug|d:s' => \$debug);


my $local_time = localtime();
print "\n\n[$local_time] Starting long-read-based spliced mapping for batch $batch_id\n";

my $output_dir = "$batch_id";
my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
system("mkdir $output_dir");

my $sample_count = 0;
my $mapping_count = 0;

my $ref_genome_file = "$base_dir/$ref_dir/ref.genome.fa";
my $ref_transcriptome_file = "$base_dir/$ref_dir/ref.transcriptome.fa";

if ($mapping_mode eq "map2genome") {
    print "Check the specified reference genome file:\n";
    if (-e $ref_genome_file) {
	print "Successfully located the specified reference genome file: $ref_genome_file.\n";
    } else {
	print "Cannot find the specified reference genome file: $ref_genome_file!\n";
	print "Exit!\n";
	exit;
    }
} elsif ($mapping_mode eq "map2transcriptome") {
    print "Check the specified reference transcriptome file:\n";
    if (-e $ref_transcriptome_file) {
	print "Successfully located the specified reference transcriptome file: $ref_transcriptome_file.\n";
    } else {
	print "Cannot find the specified reference transcriptome file: $ref_transcriptome_file!\n";
	print "Exit!\n";
	exit;
    }
} else {
    print "Error! Unsupported mapping mode: $mapping_mode!\n";
    print "please specify mapping mode as \"map2genome\" or \"map2transcritpome\"!\n";
}


foreach my $sample_id (@sample_table) {
    $local_time = localtime();
    print "\n[$local_time] Processing mapping for sample $sample_id ..\n";
    $sample_count++;
    
    my $basecalled_fastq_file = "$base_dir/$long_reads_dir/$sample_table{$sample_id}{'basecalled_fastq_file'}";
    print "Check the specified long read file:\n";
    if (-e $basecalled_fastq_file) {
	print "Successfully located the specified long read file: $basecalled_fastq_file.\n";
    } else {
	print "Cannot find the specified long read file: $basecalled_fastq_file!\n";
	print "Exit!\n";
	exit;
    }
    print "processing sample $sample_id\n";
    print "basecalled_fastq_file = $basecalled_fastq_file\n";
    my $sample_output_dir;
    $mapping_count++;
    $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    print("mapping the reads by $long_read_mapper\n");
    if ($mapping_mode eq "map2genome") {
	print "mapping $basecalled_fastq_file to reference genome: $ref_genome_file\n";
    } elsif ($mapping_mode eq "map2transcriptome") {
	print "mapping $basecalled_fastq_file to reference transcriptome: $ref_transcriptome_file\n";
    } else {
	die "unsupported mapping mode: $mapping_mode\n";
    }
    if ($long_read_mapper eq "minimap2") {
	if ($long_read_technology eq "pacbio_isoseq") {
	    if ($mapping_mode eq "map2genome") {
		system("$minimap2_dir/minimap2 -t $threads -ax splice:hq -uf --secondary=no $ref_genome_file $basecalled_fastq_file | $samtools_dir/samtools view --threads $threads -F260 -hbS -T $ref_genome_file -q $min_mapping_quality - >$sample_id.bam");
		# -F 260  output primary aligned mapped reads
                # read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
		# bit 3 + bit 9 = 4 + 256 = 260
	    } elsif ($mapping_mode eq "map2transcriptome") {
		system("$minimap2_dir/minimap2 -t $threads -ax map-hifi -p 0.99 -uf --secondary=no $ref_transcriptome_file $basecalled_fastq_file | $samtools_dir/samtools view --threads $threads -F260 -hbS -T $ref_transcriptome_file -q $min_mapping_quality - >$sample_id.bam");
		# -F 260  output primary aligned mapped reads
                # read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
		# bit 3 + bit 9 = 4 + 256 = 260
	    }
	} else {
	    # nanopore_drs
	    if ($mapping_mode eq "map2genome") {
		system("$minimap2_dir/minimap2 -t $threads -ax splice -uf -k14 --secondary=no $ref_genome_file $basecalled_fastq_file | $samtools_dir/samtools view --threads $threads -F260 -hbS -T $ref_genome_file -q $min_mapping_quality - >$sample_id.bam");
		# -F 260  output primary aligned mapped reads
                # read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
		# bit 3 + bit 9 = 4 + 256 = 260
	    } elsif ($mapping_mode eq "map2transcriptome") {
		system("$minimap2_dir/minimap2 -t $threads -ax map-ont -p 0.99 -uf --secondary=no $ref_transcriptome_file $basecalled_fastq_file | $samtools_dir/samtools view --threads $threads -F260 -hbS -T $ref_transcriptome_file -q $min_mapping_quality - >$sample_id.bam");
		# -F 260  output primary aligned mapped reads
                # read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
		# bit 3 + bit 9 = 4 + 256 = 260
            }
	}
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
	my $bam_index_statu = system("$samtools_dir/samtools index $sample_id.sort.bam");
	if ($bam_index_statu != 0) {
		print("\n Failed to create a bai index. Try using a csi index. \n");
		system("$samtools_dir/samtools index -c $sample_id.sort.bam");
	}
	# generate samtools mpileup 
	# system("$samtools_dir/samtools mpileup -Q 0 -C 50 -q $min_mapping_quality -f $base_dir/$ref_dir/ref.genome.fa $sample_id.sort.bam |gzip -c >${sample_id}.mpileup.gz");	
	# compute basic alignment statistics by samtools
	$local_time = localtime();
	print "\n[$local_time] Computing the BAM statistics for sample $sample_id ..\n";

	system("$samtools_dir/samtools flagstat $sample_id.sort.bam >$sample_id.samstat");
	# # calculate per-base depth
	# system("$samtools_dir/samtools depth -aa $sample_id.sort.bam |gzip -c >$sample_id.depth.txt.gz");
	# calculate read mapping coverage statistics
	if ($mapping_mode eq "map2genome") {
	    system("perl $NANOTRANS_HOME/scripts/summarize_mapping_coverage.pl -r $base_dir/$ref_dir/ref.genome.fa -s $sample_id.samstat -m lite -t $sample_id -o $sample_id.coverage_summary.txt");
	} elsif ($mapping_mode eq "map2transcriptome") {
            system("perl $NANOTRANS_HOME/scripts/summarize_mapping_coverage.pl -r $base_dir/$ref_dir/ref.transcriptome.fa -s $sample_id.samstat -m lite -t $sample_id -o $sample_id.coverage_summary.txt");
	}
	if ($mapping_count == 1) {
	    system("cp $sample_id.coverage_summary.txt ./../all_samples.coverage_summary.txt");
	} else {
	    system("tail -1 $sample_id.coverage_summary.txt  >> ./../all_samples.coverage_summary.txt");
	}
	system("rm -r tmp");
	$local_time = localtime();
	print "\n[$local_time] Finishing long-read-based spliced mapping for sample $sample_id.\n";	
    }
    chdir("./../") or die "cannot change directory to: $!\n";
}

$local_time = localtime();
print "\n[$local_time] Finishing long-read-based spliced mapping for batch $batch_id.\n";
print "\nA total of $sample_count samples were processed!\n";
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

