#! /usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use List::Util qw(first);

##############################################################
#  script: tidy_diffExp_die_drimseq_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.10.30
#  description: tidy up flair's diffExp diu drimseq output file
#  example: perl tidy_diffExp_die_drimseq_output.pl -sample_table Master_Sample_Table.txt -i input.txt(.gz) -o output.txt(.gz) -m count_matrix.tsv  -x transcript2gene_map
##############################################################

my ($sample_table, $input, $output, $filtered_count_matrix, $transcript2gene_map, $master_sample_table);
GetOptions('sample_table|sample_table:s' => \$sample_table,
	   'input|i:s' => \$input, 
	   'output|o:s' => \$output,
	   'count_matrix|m:s' => \$filtered_count_matrix,
           'transcript2gene_map|x:s' => \$transcript2gene_map);

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my %comparison_groups = ();
foreach my $sample_id (@sample_table) {
    my $g = $sample_table{$sample_id}{'comparison_group'};
    if (not exists $comparison_groups{$g}) {
        @{$comparison_groups{$g}} = ($sample_id);
    } else {
        push @{$comparison_groups{$g}}, $sample_id;
    }
}

my $filtered_count_matrix_fh = read_file($filtered_count_matrix);
my %filtered_count = ();
parse_filtered_count_matrix_file($filtered_count_matrix_fh, \%filtered_count);

#print "input=$input\n";
my ($comparison_group_a, $comparison_group_b) = ($input =~ /(?:genes_deseq2|isoforms_deseq2|isoforms_drimseq)\_([^\-]+)_v_([^\-]+)\_(?:results)/);
# print "comparison_group_a: $comparison_group_a\n";
# print "comparison_group_b: $comparison_group_b\n";

my $sample_count = scalar @sample_table;
my @sample_tags = ();
foreach my $s (@{$comparison_groups{$comparison_group_a}}) {
    push @sample_tags, $s;
}
foreach my $s (@{$comparison_groups{$comparison_group_b}}) {
    push @sample_tags, $s;
}

# print "sample_tags: @sample_tags\n";

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
    /feature_id\s+gene_id/ and next;
    my @line = split /\t/, $_;
    my $line_length = scalar @line;
    my $contig_id = $line[1];
    my ($isoform_id, $gene_id) = split /_/, $contig_id; 
    my $lr = "NA";
    my $df = "NA";
    my $p_value = 9999;
    my $adj_p_value = 9999;
    # print "line_length = $line_length\n";
    if ($line_length == 15) {
	if ($line[11] ne "") {
	    $lr = $line[11];
	}
	if ($line[12] ne "") {
	    $df = int $line[12];
	}
	if ($line[13] ne "") {
            $p_value = $line[13];
        } 
	if ($line[14] ne "") {
	    $adj_p_value = $line[14];
	}
    }
    my @sample_values = ();
    foreach my $sample_id (@sample_tags) {
	# print "sample_id=$sample_id, contig_id=$contig_id\n";
	# print "count=$filtered_count{$contig_id}{$sample_id}\n";
	push @sample_values, $filtered_count{$contig_id}{$sample_id};
    }

    my $sample_values = join "\t", @sample_values;
    $i++;
    $input{$p_value}{$i}{'gene_id'} = $gene_id;
    $input{$p_value}{$i}{'isoform_id'} = $isoform_id;
    $input{$p_value}{$i}{'lr'} = $lr,
    $input{$p_value}{$i}{'df'} = $df,
    $input{$p_value}{$i}{'adj_p_value'} = $adj_p_value;
    $input{$p_value}{$i}{'sample_values'} = $sample_values;
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
	my $lr = $input{$p_value}{$i}{'lr'};
	my $df = $input{$p_value}{$i}{'df'};
	my $adj_p_value = $input{$p_value}{$i}{'adj_p_value'};
	my $sample_tags = join "\t", @sample_tags;
	my $sample_values = join "\t", $input{$p_value}{$i}{'sample_values'};
	if ($flag == 0) {
	    print $output_fh "isoform_id\tgene_id\tgene_name\tcomparison_group_pair\tlog_likelihood_ratio\tdf\tp_value\tadj_p_value\t$sample_tags\n";
	    $flag = 1;
	}
	if ($flag == 1) {
	    if ($p_value == 9999) {
		print $output_fh "$isoform_id\t$gene_id\t$gene_name\t${comparison_group_b}\/${comparison_group_a}\tNA\tNA\tNA\tNA\t$sample_values\n";
	    } else {
		print $output_fh "$isoform_id\t$gene_id\t$gene_name\t${comparison_group_b}\/${comparison_group_a}\t$lr\t$df\t$p_value\t$adj_p_value\t$sample_values\n";
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

sub parse_filtered_count_matrix_file {
    my ($fh, $filtered_count_matrix_hashref) = @_;
    my @raw_sample_tags = ();
    my $sample_count;
    while (<$fh>) {
	chomp;
	/^\s*$/ and next;
        /^\s*#/ and next;
	if ($. == 1) {
	    my @header = split "\t", $_;
	    $sample_count = (scalar @header) - 3;
	    @raw_sample_tags = splice @header, 3, $sample_count;
	    # print "sample_count: $sample_count, raw_sample_tags: @raw_sample_tags\n";
	} 
	else {
	    my @count = split "\t", $_;
	    my $contig = $count[2];
	    for (my $i = 0; $i < $sample_count; $i++) {
		my ($sample_id, $comparison_group, $replicate_id, $subfix) = split /\_/, $raw_sample_tags[$i]; 
		# print "i=$i, sample_id=$sample_id, count=$count[$i+3]\n";
		$$filtered_count_matrix_hashref{$contig}{$sample_id} = $count[$i+3];
	    }
	}
    }
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
