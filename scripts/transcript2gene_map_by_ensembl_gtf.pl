#!/usr/bin/perl
#use warnings FATAL => 'all';
use warnings;
use strict;
use Getopt::Long;

my ($input, $output, $ignore_version_number);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'ignore_version_number|v:s' => \$ignore_version_number);

my $input_fh = read_file($input);
my $output_fh = write_file($output);
my %seen = ();

print $output_fh "transcript_id\tgene_id\tgene_name\n";

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $_;
    if (($attribute =~ /transcript_id/) and ($attribute =~ /gene_id/)) {
	my ($transcript_id) = ($attribute =~ /transcript_id \"([^\"]+)\"/);
	if (not exists $seen{$transcript_id}) {
	    my ($gene_id) = ($attribute =~ /gene_id \"([^;]+)\"/);
	    if ($gene_id =~ /(\S+)\.\d+$/) {
		$gene_id = $1;
	    }
	    if ($attribute =~ /gene_name/) {
		# print "line=$_\n";
		my ($gene_name) = ($attribute =~ /gene_name \"([^\"]+)\"/);
		print $output_fh "$transcript_id\t$gene_id\t$gene_name\n";
	    } else {
		print $output_fh "$transcript_id\t$gene_id\t$gene_id\n";
	    }
	    $seen{$transcript_id} = 1;
	}
    } 
}

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
