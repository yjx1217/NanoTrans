#! /usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
#use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);

##############################################################
#  script: tidy_diff_iso_usage_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.09.26
#  description: tidy up flair's tidy_diff_iso_usage.py output file
#  example: perl tidy_diff_iso_usage_output.pl -i input.txt(.gz) -o output.txt(.gz) -m count_matrix.tsv -a sample_a -b sample_b -x transcript2gene_map
##############################################################

my ($input, $output, $sample_a, $sample_b, $transcript2gene_map);
GetOptions('input|i:s' => \$input, # input genome fasta file
	   'output|o:s' => \$output,
	   'sample_a|a:s' => \$sample_a,
	   'sample_b|b:s' => \$sample_b,
           'transcript2gene_map|x:s' => \$transcript2gene_map);

my $transcript2gene_map_fh = read_file($transcript2gene_map);
my %transcript2gene_map = ();
my %gene_id2name = ();
while (<$transcript2gene_map_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    /^transcript_id\tgene_id\tgene_name/ and next;
    my ($transcript_id, $gene_id, $gene_name) = split /\t/, $_;
    $transcript2gene_map{$transcript_id}{'gene_id'} = $gene_id;
    $transcript2gene_map{$transcript_id}{'gene_name'} = $gene_name;
    if (not exists $gene_id2name{$gene_id}) {
        $gene_id2name{$gene_id} = $gene_name;
    }
}

my $input_fh = read_file($input);
my $output_fh = write_file($output);
my %input = ();

my $i = 0;
my @p_values = ();
while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    /^\s*#/ and next;
    my ($gene_id, $isoform_id, $p_value, $sample_a_isoform_count, $sample_b_isoform_count, $sample_a_other_isoforms_count, $sample_b_other_isoforms_count) = split /\t/, $_;
    if ($p_value eq "NA") {
	$p_value = 9999;
    }
    $i++;
    $input{$p_value}{$i}{'gene_id'} = $gene_id;
    $input{$p_value}{$i}{'isoform_id'} = $isoform_id;
    $input{$p_value}{$i}{'sample_a_isoform_count'} = $sample_a_isoform_count;
    $input{$p_value}{$i}{'sample_b_isoform_count'} = $sample_b_isoform_count;
    $input{$p_value}{$i}{'sample_a_other_isoforms_count'} = $sample_a_other_isoforms_count;
    $input{$p_value}{$i}{'sample_b_other_isoforms_count'} = $sample_b_other_isoforms_count;
}

my $flag = 0;
foreach my $p_value (sort {$a <=> $b} keys %input) {
    foreach my $i (sort {$a <=> $b} keys %{$input{$p_value}}) {
	my $gene_id = $input{$p_value}{$i}{'gene_id'};
	my $isoform_id = $input{$p_value}{$i}{'isoform_id'};
	my $gene_name = $gene_id;
	if (exists $transcript2gene_map{$isoform_id}) {
	    $gene_name = $transcript2gene_map{$isoform_id}{'gene_name'};
	} elsif (exists $gene_id2name{$gene_id}) {
	    $gene_name = $gene_id2name{$gene_id};
	}
	my $sample_a_isoform_count = $input{$p_value}{$i}{'sample_a_isoform_count'};
	my $sample_b_isoform_count = $input{$p_value}{$i}{'sample_b_isoform_count'};
	my $sample_a_other_isoforms_count = $input{$p_value}{$i}{'sample_a_other_isoforms_count'};
	my $sample_b_other_isoforms_count = $input{$p_value}{$i}{'sample_b_other_isoforms_count'};
	my ($sample_a_lite) = ($sample_a =~ /^([^\_]+)_/);
	my ($sample_b_lite) = ($sample_b =~ /^([^\_]+)_/);

	my $odds_ratio;
	if ($sample_b_other_isoforms_count * $sample_a_isoform_count == 0) {
	    $odds_ratio = "NA";
	} else {
	    $odds_ratio = ($sample_b_isoform_count * $sample_a_other_isoforms_count)/($sample_b_other_isoforms_count * $sample_a_isoform_count);
	}
	if ($flag == 0) {
	    print $output_fh "isoform_id\tgene_id\tgene_name\tcomparison_sample_pair\todds_ratio\tp_value\tadj_p_value\t${sample_a_lite}_isoform_count\t${sample_b_lite}_isoform_count\t${sample_a_lite}_other_isoforms_count\t${sample_b_lite}_other_isoforms_count\n";
	    $flag = 1;
	}
	if ($flag == 1) {
	    if ($p_value == 9999) {
		print $output_fh "$isoform_id\t$gene_id\t$gene_name\t$sample_b_lite/$sample_a_lite\tNA\tNA\tNA\t$sample_a_isoform_count\t$sample_a_isoform_count\t$sample_a_other_isoforms_count\t$sample_b_other_isoforms_count\n";
	    } else {
		my $adj_p_value = "NA";
		print $output_fh "$isoform_id\t$gene_id\t$gene_name\t$sample_b_lite/$sample_a_lite\t$odds_ratio\t$p_value\t$adj_p_value\t$sample_a_isoform_count\t$sample_a_isoform_count\t$sample_a_other_isoforms_count\t$sample_b_other_isoforms_count\n";
	    }
	}
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

sub parse_count_matrix_file {
    my $fh = shift @_;
    my @sample_tags = ();
    my $flag = 0;
    while (<$fh>) {
	chomp;
	/^\s*$/ and next;
        /^\s*#/ and next;
	if ((/^ids\s/) and ($flag == 0)) {
	   my @header = split "\t", $_;
	   my $sample_count = (scalar @header) - 1;
	   @sample_tags = splice @header, 1, $sample_count;
	   $flag = 1;
	} 
    }
    return @sample_tags;
}


sub parse_input_file {
    my $fh = shift @_;
    my %input = ();
    my @header = ();
    my @sample_tags = ();
    my $sample_count;
    my $comparison_sample_pair;
    my $line_count = 0;
    while (<$fh>) {
        chomp;
	/^\s*$/ and next;
	/^\s*#/ and next;
	if (/^feature_id\s+coordinate/) {
	    @header = split /\t/, $_;
	    $sample_count = (scalar @header) -4;
	    $comparison_sample_pair = $header[-1];
	    $comparison_sample_pair =~ s/_pval//gi;
	    @sample_tags = splice @header, 2, $sample_count;	    
	    print "sample_tags: @sample_tags\n";
	} else {
	    my @line = split /\t/, $_;
	    my $feature_id = $line[0];
	    my $coordinate = $line[1];
	    my $isoform_ids = $line[-2];
	    my $p_value = $line[-1];
	    my @sample_values = splice @line, 2, $sample_count;
	    $line_count++;
	    $input{$p_value}{$line_count}{$feature_id}{'coordinate'} = $coordinate;
	    $input{$p_value}{$line_count}{$feature_id}{'isoform_ids'} = $isoform_ids;
	    $input{$p_value}{$line_count}{$feature_id}{'comparison_sample_pair'} = $comparison_sample_pair;
	    $input{$p_value}{$line_count}{$feature_id}{'sample_tags'} = join "\t", @sample_tags;
	    $input{$p_value}{$line_count}{$feature_id}{'sample_values'} = join "\t", @sample_values;
	}
    }
    return %input;
}

