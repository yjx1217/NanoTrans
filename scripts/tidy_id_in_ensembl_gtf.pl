#!/usr/bin/perl
#use warnings FATAL => 'all';
use warnings;
use strict;
use Getopt::Long;

my ($input, $output);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output);

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    if (/^#/) {
	print $output_fh "$_\n";
    } elsif (/^\s*$/) {
	print $output_fh "$_\n";
    } else {
	my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $_;
	my $mod_attribute = "";
        my @attribute = split ";", $attribute;
        foreach my $a (@attribute) {
            $a =~ s/^\s+//gi;
            # print "attribute_a: $a\n";
            my $feature_id;
            my $feature_type;
            if ($a =~ /(gene_id|transcript_id|exon_id) \"([^\"]+)\"/) {
                $feature_type = $1;
                $feature_id = $2;
                # print "feature_id = $feature_id\n";
		$feature_id =~ s/\(/\_/gi;
		$feature_id =~ s/\)/\_/gi;
                # print "mod_feature_id = $feature_id\n";
                $a = "$feature_type \"$feature_id\"";
            }
            $mod_attribute .= "$a; ";
        }
        $mod_attribute =~ s/\s+$//gi;
        print $output_fh "$chr\t$source\t$feature\t$start\t$end\t$score\t$strand\t$frame\t$mod_attribute\n";
    }
}
close $input_fh;
close $output_fh;


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
