#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

##############################################################
#  script: tidy_diffsplice_drimseq_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.11.14
#  description: tidy up flair's diffsplice drimseq output file
#  example: perl tidy_diffsplice_drimseq_output.pl -i input.tsv(.gz) -q quant.tsv(.gz) -o output.txt(.gz) -a $group_a -b $group_b
##############################################################

my ($input, $output, $quant, $group_a, $group_b, $transcript2gene_map);
GetOptions('input|i:s' => \$input, 
	   'quant|q:s' => \$quant,
	   'output|o:s' => \$output,
	   'group_a|a:s' => \$group_a,
	   'group_b|b:s' => \$group_b,
	   'transcript2gene_map|x:s' => \$transcript2gene_map);

my $input_fh = read_file($input);
my %input = parse_input_file($input_fh);

my $quant_fh = read_file($quant);
my %quant = parse_quant_file($quant_fh);

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


my $output_fh = write_file($output);

my $flag = 0;
foreach my $p_value (sort {$a <=> $b} keys %input) {
    foreach my $i (sort {$a <=> $b} keys %{$input{$p_value}}) {
	foreach my $feature_id (sort keys %{$input{$p_value}{$i}}) {
	    my $coordinate = $input{$p_value}{$i}{$feature_id}{'coordinate'};
	    my $isoform_ids = $quant{$feature_id}{'isoform_ids'};
	    my $comparison_group_pair = "${group_a}-${group_b}";
	    my $sample_tags = $input{$p_value}{$i}{$feature_id}{'sample_tags'};
	    my $sample_values = $input{$p_value}{$i}{$feature_id}{'sample_values'};
	    my $lr = $input{$p_value}{$i}{$feature_id}{'lr'};
	    my $df = $input{$p_value}{$i}{$feature_id}{'df'};
	    my $adj_p_value = $input{$p_value}{$i}{$feature_id}{'adj_p_value'};
	    if ($flag == 0) {
		print $output_fh "isoform_id\tgene_id\tgene_name\tfeature_id\tcoordinate\tcomparison_group_pair\tlog_likelihood_ratio\tdf\tp_value\tadj_p_value\t$sample_tags\n";
		$flag = 1;
	    }
	    if ($flag == 1) {
		my @isoform_id_array = ();
		if ($isoform_ids =~ /,/) {
		    @isoform_id_array = split /,/, $isoform_ids;
		} else {
		    @isoform_id_array = ($isoform_ids);
		}
		my @gene_id_array = ();
		my @gene_name_array = ();
		foreach my $isoform_id (@isoform_id_array) {
		    my $gene_id = "NA";
		    my $gene_name = "NA";
		    if (exists $transcript2gene_map{$isoform_id}) {
			$gene_id = $transcript2gene_map{$isoform_id}{'gene_id'};
			$gene_name = $transcript2gene_map{$isoform_id}{'gene_name'};
		    } 
		    if (($gene_name eq "NA") and (exists $gene_id2name{$gene_id})) {
			$gene_name = $gene_id2name{$gene_id};
		    }
		    push @gene_id_array, $gene_id;
		    push @gene_name_array, $gene_name;
		}
		my $isoform_id_list = join ",", @isoform_id_array;
		my $gene_id_list = join ",", @gene_id_array;
		my $gene_name_list = join ",", @gene_name_array;
		if ($p_value ne "9999") {
		    print $output_fh "$isoform_id_list\t$gene_id_list\t$gene_name_list\t$feature_id\t$coordinate\t$comparison_group_pair\t$lr\t$df\t$p_value\t$adj_p_value\t$sample_values\n";
		} else {
		    print $output_fh "$isoform_id_list\t$gene_id_list\t$gene_name_list\t$feature_id\t$coordinate\t$comparison_group_pair\t$lr\t$df\tNA\tNA\t$sample_values\n";
		}
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
	if (/feature_id\s+gene_id/) {
	    @header = split /\t/, $_;
	    $sample_count = (scalar @header) - 7;
	    @sample_tags = splice @header, 3, $sample_count;
	    foreach my $s (@sample_tags) {
		($s) = ($s =~ /^([^\_]+)\_/);
	    }
	    print "sample_tags: @sample_tags\n";
	} else {
	    my @line = split /\t/, $_;
	    my $feature_id = $line[1];
	    my $coordinate = $line[2];
	    my $lr = $line[-4];
	    my $df = $line[-3];
	    my $p_value = $line[-2];
	    my $adj_p_value = $line[-1];
	    my @sample_values = splice @line, 3, $sample_count;
	    $line_count++;
	    if (($p_value eq "") or ($p_value eq "NA")) {
		$p_value = 9999;
	    }
	    $input{$p_value}{$line_count}{$feature_id}{'coordinate'} = $coordinate;
	    $input{$p_value}{$line_count}{$feature_id}{'lr'} = $lr;
	    $input{$p_value}{$line_count}{$feature_id}{'df'} = $df;
	    $input{$p_value}{$line_count}{$feature_id}{'p_value'} = $p_value;
	    $input{$p_value}{$line_count}{$feature_id}{'adj_p_value'} = $adj_p_value;
	    $input{$p_value}{$line_count}{$feature_id}{'sample_tags'} = join "\t", @sample_tags;
	    $input{$p_value}{$line_count}{$feature_id}{'sample_values'} = join "\t", @sample_values;
	}
    }
    return %input;
}

sub parse_quant_file {
    my $fh = shift @_;
    my %quant = ();
    my @header = ();
    my @sample_tags = ();
    my $sample_count;
    my $line_count = 0;
    while (<$fh>) {
        chomp;
	/^\s*$/ and next;
	/^\s*#/ and next;
	if (/feature_id\s+coordinate/) {
	    @header = split /\t/, $_;
	    $sample_count = (scalar @header) -3;
	    @sample_tags = splice @header, 2, $sample_count;
	    foreach my $s (@sample_tags) {
		($s) = ($s =~ /^([^\_]+)\_/);
	    }
	    # print "sample_tags: @sample_tags\n";
	} else {
	    my @line = split /\t/, $_;
	    my $feature_id = $line[0];
	    my $coordinate = $line[1];
	    my $isoform_ids = $line[-1];
	    my @sample_values = splice @line, 2, $sample_count;
	    $line_count++;
	    $quant{$feature_id}{'coordinate'} = $coordinate;
	    $quant{$feature_id}{'isoform_ids'} = $isoform_ids;
	    $quant{$feature_id}{'sample_tags'} = join "\t", @sample_tags;
	    $quant{$feature_id}{'sample_values'} = join "\t", @sample_values;
	}
    }
    return %quant;
}

