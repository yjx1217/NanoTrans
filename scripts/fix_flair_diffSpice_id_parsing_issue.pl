#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

##############################################################
#  script: tidy_diffsplice_fishers_exact_output.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.09.26
#  description: tidy up flair's diffsplice_fishers_exact.py output file
#  example: perl tidy_diffsplice_fishers_exact_output.pl -i input.txt(.gz) -o output.txt(.gz)
##############################################################

my ($input, $output);
GetOptions('input|i:s' => \$input, 
	   'output|o:s' => \$output);

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    if (/^feature_id\tcoordinate/) {
	print $output_fh "$_\n";
    } else {
	my @line = split /\t/, $_;
	my @isoform_ids = split /,/, $line[-1];
	my @new_isoform_ids = ();
	foreach my $isoform_ids (@isoform_ids) {
	    my $new_isoform_ids = $isoform_ids;
	    if ($isoform_ids =~ /\_/) {
		($new_isoform_ids) = ($isoform_ids =~ /([^_]+)\_/);
	    }
	    push @new_isoform_ids, $new_isoform_ids;
	}
	$line[-1] = join ",", @new_isoform_ids;
	my $new_line = join "\t", @line;
	print $output_fh "$new_line\n";
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

