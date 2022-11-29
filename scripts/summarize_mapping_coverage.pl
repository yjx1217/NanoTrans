#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Descriptive::Discrete;
use Getopt::Long;

##############################################################
#  script: summarize_mapping_coverage.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.08.26
#  description: summarize read mapping coverage statistics based on the samtools samstat and depth outputs. 
#  example: perl summarize_mapping_coverage.pl -t sample_tag -r refseq.fa(.gz) -s samstat.out(.gz) -m normal  -d depth.out(.gz) -c 5 -o sample.coverage_summary.txt
##############################################################

my $sample_tag;
my $refseq;
my $samstat;
my $depth;
my $min_depth_cutoff = 5;
my $mode = "normal"; # "normal" or "lite"
my $output;

GetOptions('sample_tag|t:s' => \$sample_tag,
	   'refseq|r:s' => \$refseq,
	   'samstat|s:s' => \$samstat,
	   'depth|d:s' => \$depth,
	   'mode|m:s' => \$mode,
	   'min_depth_cutoff|c:f' => \$min_depth_cutoff,
	   'output|o:s' => \$output);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);
close $refseq_fh;

my $samstat_fh = read_file($samstat);
my %samstat_stat = parse_samstat_file($samstat_fh);
close $samstat_fh;

my $depth_fh;
my %depth_stat;

my $output_fh = write_file($output);

if ($mode eq "normal") {
    $depth_fh = read_file($depth);
    %depth_stat = parse_depth_file_Nfiltered($depth_fh, $min_depth_cutoff, \%refseq);
    close $depth_fh;
}

if ($mode eq "normal") {
    print $output_fh "sample\tref_genome\ttotal_reads\tmapping_rate\tpairing_rate\t(>=$min_depth_cutoff)_depth_proportion\tzero_depth_proportion\tmean_depth\tmedian_depth";
    foreach my $chr (@refseq) {
	print $output_fh "\t${chr}_mean_depth\t${chr}_median_depth";
    }
    print $output_fh "\n";
} else {
 print $output_fh "sample\tref_genome\ttotal_reads\tmapping_rate\tpairing_rate\n";
}   

if ($mode eq "normal") {
    print $output_fh "$sample_tag\t$refseq\t$samstat_stat{'total_reads'}\t$samstat_stat{'mapping_rate'}\t$samstat_stat{'pairing_rate'}";
    print $output_fh "\t$depth_stat{'min_depth_proportion'}\t$depth_stat{'zero_depth_proportion'}";
    print $output_fh "\t$depth_stat{'mean_depth'}\t$depth_stat{'median_depth'}";
    foreach my $chr (@refseq) {
	print $output_fh "\t$depth_stat{'chr_mean_depth'}{$chr}\t$depth_stat{'chr_median_depth'}{$chr}";
    }
    print $output_fh "\n";
} else {
    print $output_fh "$sample_tag\t$refseq\t$samstat_stat{'total_reads'}\t$samstat_stat{'mapping_rate'}\t$samstat_stat{'pairing_rate'}";
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

sub parse_fasta_file {
    my ($fh, $input_hashref, $input_arrayref) = @_;
    my $seq_name = "";
    while (<$fh>) {
	chomp;
	if (/^\s*$/) {
	    next;
	} elsif (/^\s*#/) {
	    next;
	} elsif (/^>(.*)/) {
	    $seq_name = $1;
	    push @$input_arrayref, $seq_name;
	    $$input_hashref{$seq_name} = "";
	} else {
	    $$input_hashref{$seq_name} .= $_;
	}
    }
}

sub parse_samstat_file {
    my $fh = shift @_;
    my $total_reads = 0;
    my $qc_passed_reads = 0;
    my $qc_failed_reads = 0;
    my $qc_passed_mapped_reads = 0;
    my $qc_passed_mapped_paired_reads = 0;
    my $mapping_rate = 0;
    my $pairing_rate = 0;
    my %stat = ();
    while (<$fh>) {
        chomp;
	if (/in total\s+\(/) {
	    ($qc_passed_reads, $qc_failed_reads) = ($_ =~ /^(\S+) \+ (\S+) in total \(QC\-passed reads \+ QC\-failed reads\)/);
	    $total_reads = $qc_passed_reads + $qc_failed_reads;
	} elsif (/mapped\s+\(/) {
	    ($qc_passed_mapped_reads) = ($_ =~ /^(\S+) \+ \S+/);
	} elsif (/properly paired\s+\(/) {
	    ($qc_passed_mapped_paired_reads) = ($_ =~ /^(\S+) \+ \S+/);
	}
    }
    if ($total_reads != 0) {
	$mapping_rate = $qc_passed_mapped_reads/$total_reads;
	$pairing_rate = $qc_passed_mapped_paired_reads/$total_reads;
    }
    $stat{'total_reads'} = $total_reads;
    $stat{'mapping_rate'} = sprintf("%.2f", $mapping_rate);
    $stat{'pairing_rate'} = sprintf("%.2f", $pairing_rate);
    return %stat;
}

sub parse_depth_file {
    my ($fh, $min_depth_cutoff) = @_;
    my $zero_depth_count = 0;
    my $min_depth_count = 0;
    my $total_base_count = 0;
    my %depth = ();
    my %stat = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $depth) = split /\t/, $_;
	if (not exists $depth{$chr}{$depth}) {
	    $depth{$chr}{$depth} = 1;
	} else {
	    $depth{$chr}{$depth}++;
	}
	$total_base_count++;
	if ($depth == 0) {
	    $zero_depth_count++;
	} elsif ($depth >= $min_depth_cutoff) {
	    $min_depth_count++;
	}
    }
    my $zero_depth_proportion = 0;
    my $min_depth_proportion = 0;
    if ($total_base_count != 0) {
	$zero_depth_proportion = $zero_depth_count/$total_base_count;
	$min_depth_proportion = $min_depth_count/$total_base_count;
    }

    $stat{"zero_depth_proportion"} = sprintf("%.2f", $zero_depth_proportion);
    $stat{"min_depth_proportion"} = sprintf("%.2f", $min_depth_proportion);

    my @total_depth_data_tuple = ();
    my %chr_depth_data_tuple = ();
    foreach my $chr (sort keys %depth) {
	foreach my $depth (keys %{$depth{$chr}}) {
	    my $freq = $depth{$chr}{$depth};
	    push @total_depth_data_tuple, $depth;
	    push @total_depth_data_tuple, $freq;
	    if (exists $chr_depth_data_tuple{$chr}) {
		push @{$chr_depth_data_tuple{$chr}}, $depth;
		push @{$chr_depth_data_tuple{$chr}}, $freq;
	    } else {
		@{$chr_depth_data_tuple{$chr}} = ($depth, $freq);
	    }
	}
    }

    my $total_depth_stat = Statistics::Descriptive::Discrete->new();
    $total_depth_stat->add_data_tuple(@total_depth_data_tuple);
    my $mean_depth = $total_depth_stat->mean();
    my $median_depth = $total_depth_stat->median();
    $stat{"mean_depth"} = sprintf("%.2f", $mean_depth);
    $stat{"median_depth"} = sprintf("%.2f", $median_depth);
    @total_depth_data_tuple = ();
    foreach my $chr (sort keys %depth) {
	my $chr_depth_stat = Statistics::Descriptive::Discrete->new();
	$chr_depth_stat->add_data_tuple(@{$chr_depth_data_tuple{$chr}});
	my $chr_mean_depth = $chr_depth_stat->mean();
	my $chr_median_depth = $chr_depth_stat->median();
	$stat{"chr_mean_depth"}{$chr} = sprintf("%.2f", $chr_mean_depth);
	$stat{"chr_median_depth"}{$chr} = sprintf("%.2f", $chr_median_depth);
	@{$chr_depth_data_tuple{$chr}} = ();
    }
    return %stat;
}

sub parse_depth_file_Nfiltered {
    my ($fh, $min_depth_cutoff, $genome_hashref) = @_;
    my $zero_depth_count = 0;
    my $min_depth_count = 0;
    my $total_base_count = 0;
    my %depth = ();
    my %stat = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $depth) = split /\t/, $_;
	# print "chr=$chr, pos=$pos\n";
	my $base = substr $$genome_hashref{$chr}, $pos - 1, 1;
	if ($base !~ /(N|n)/) {
	    if (not exists $depth{$chr}{$depth}) {
		$depth{$chr}{$depth} = 1;
	    } else {
		$depth{$chr}{$depth}++;
	    }
	    $total_base_count++;
	    if ($depth == 0) {
		$zero_depth_count++;
	    } elsif ($depth >= $min_depth_cutoff) {
		$min_depth_count++;
	    }
	}
    }
    my $zero_depth_proportion = 0;
    my $min_depth_proportion = 0;
    if ($total_base_count != 0) {
	$zero_depth_proportion = $zero_depth_count/$total_base_count;
	$min_depth_proportion = $min_depth_count/$total_base_count;
    }

    $stat{"zero_depth_proportion"} = sprintf("%.2f", $zero_depth_proportion);
    $stat{"min_depth_proportion"} = sprintf("%.2f", $min_depth_proportion);

    my @total_depth_data_tuple = ();
    my %chr_depth_data_tuple = ();
    foreach my $chr (sort keys %depth) {
	foreach my $depth (keys %{$depth{$chr}}) {
	    my $freq = $depth{$chr}{$depth};
	    push @total_depth_data_tuple, $depth;
	    push @total_depth_data_tuple, $freq;
	    if (exists $chr_depth_data_tuple{$chr}) {
		push @{$chr_depth_data_tuple{$chr}}, $depth;
		push @{$chr_depth_data_tuple{$chr}}, $freq;
	    } else {
		@{$chr_depth_data_tuple{$chr}} = ($depth, $freq);
	    }
	}
    }

    my $total_depth_stat = Statistics::Descriptive::Discrete->new();
    $total_depth_stat->add_data_tuple(@total_depth_data_tuple);
    my $mean_depth = $total_depth_stat->mean();
    my $median_depth = $total_depth_stat->median();
    $stat{"mean_depth"} = sprintf("%.2f", $mean_depth);
    $stat{"median_depth"} = sprintf("%.2f", $median_depth);
    @total_depth_data_tuple = ();
    foreach my $chr (sort keys %depth) {
	my $chr_depth_stat = Statistics::Descriptive::Discrete->new();
	$chr_depth_stat->add_data_tuple(@{$chr_depth_data_tuple{$chr}});
	my $chr_mean_depth = $chr_depth_stat->mean();
	my $chr_median_depth = $chr_depth_stat->median();
	$stat{"chr_mean_depth"}{$chr} = sprintf("%.2f", $chr_mean_depth);
	$stat{"chr_median_depth"}{$chr} = sprintf("%.2f", $chr_median_depth);
	@{$chr_depth_data_tuple{$chr}} = ();
    }
    return %stat;
}



    
