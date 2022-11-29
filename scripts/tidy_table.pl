#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: tidy_table.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.04.14
#  description: tidy up a space/tab/comma-separated multi-column text file and to remove extra spaces that might be introduced during copy and paste
#  example: perl tidy_table.pl -i input.txt(.gz) -o output.txt(.gz)
##############################################################

my ($input, $output, $deliminator);
$deliminator = "space"; # "space" or "tab" or "comma"
GetOptions('input|i:s' => \$input, # input table filename
           'output|o:s' => \$output, # output table filename
           'deliminator|d:s' => \$deliminator); 

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my @line = ();
    if ($deliminator eq "space") {
	@line = split /\s+/, $_;
    } elsif ($deliminator eq "tab") {
	@line = split /\t/, $_;
    } elsif ($deliminator eq "comma") {
	@line = split /,/, $_;
    }
    my @tidy_line = ();
    foreach my $c (@line) {
	$c =~ s/^\s+//gi;
	$c =~ s/\s+$//gi;
	push @tidy_line, $c;
    }
    my $tidy_line;
    if ($deliminator eq "space") {
	$tidy_line = join " ", @tidy_line;
    } elsif ($deliminator eq "tab") {
	$tidy_line = join "\t", @tidy_line;
    } elsif ($deliminator eq "comma") {
	$tidy_line = join ",", @tidy_line;
    }
    print $output_fh "$tidy_line\n";
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

