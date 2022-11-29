#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(first);

##############################################################
#  script: tidy_diff_iso_usage_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.09.26
#  description: tidy up flair's tidy_diff_iso_usage.py output file
#  example: perl tidy_diff_iso_usage_output.pl -i input.txt(.gz) -o output.txt(.gz) -m count_matrix.tsv -a sample_a -b sample_b -x transcript2gene_map
##############################################################

my ($input, $output, $transcript2gene_map);
GetOptions('input|i:s' => \$input, # input genome fasta file
	   'output|o:s' => \$output,
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
while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    /^\s*#/ and next;
    if (/^id,position,kmer/) {
	my @header = split /,/, $_;
	my @header_to_keep = splice @header, 1;
	my $header_to_keep = join "\t", @header_to_keep;
	print $output_fh "isoform_id\tgene_id\tgene_name\t$header_to_keep\n";
    } else {
	my @content = split /,/,$_;
	my @content_to_keep = splice @content, 1;
	my $content_to_keep = join "\t",@content_to_keep;
	my ($isoform_id, $gene_id) = split "_", $content[0];
	my $p_value = $content_to_keep[3];
	my $z_score = $content_to_keep[4];
	my $z_score_abs = abs($z_score);
	if ($p_value eq "NA") {
	    $p_value = 9999;
	}
	$i++;
	$input{$z_score_abs}{$p_value}{$i}{'gene_id'} = $gene_id;
	$input{$z_score_abs}{$p_value}{$i}{'isoform_id'} = $isoform_id;
	$input{$z_score_abs}{$p_value}{$i}{'content_to_keep'} = $content_to_keep;
    }
}


my $flag = 0;
foreach my $z_score_abs (sort {$b <=> $a} keys %input) {
    foreach my $p_value (sort {$a <=> $b} keys %{$input{$z_score_abs}}) {
	foreach my $i (sort {$a <=> $b} keys %{$input{$z_score_abs}{$p_value}}) {
	    my $gene_id = $input{$z_score_abs}{$p_value}{$i}{'gene_id'};
	    my $isoform_id = $input{$z_score_abs}{$p_value}{$i}{'isoform_id'};
	    my $gene_name = $gene_id;
	    if (exists $transcript2gene_map{$isoform_id}) {
		$gene_name = $transcript2gene_map{$isoform_id}{'gene_name'};
	    } elsif (exists $gene_id2name{$gene_id}) {
		$gene_name = $gene_id2name{$gene_id};
	    }
	    my $content_to_keep = $input{$z_score_abs}{$p_value}{$i}{'content_to_keep'};
	    print $output_fh "$isoform_id\t$gene_id\t$gene_name\t$content_to_keep\n";
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

