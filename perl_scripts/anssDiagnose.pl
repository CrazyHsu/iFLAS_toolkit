#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
    Output: ID, left spliced supported,
            supported junctions count in SEs, junctions count in SEs,
            right spliced supported, skipped supported
Options:
    -j --junc   FILE    The junction file
    -a --anss   5/3     
    -h --help           Print this help information
HELP
    exit(-1);
}

my ($juncFile, $anss);
GetOptions(
            'j|junc=s'  => \$juncFile,
            'a|anss=i'  => \$anss,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open JUNC, "$juncFile" or die "Can't read file ($juncFile): $!";

my %juncHash;
while(<JUNC>){
    chomp;
    my @fields = split "\t";
    my ($chr, $strand) = @fields[0, 5];
    my @blockSizes = split ',', $fields[10];
    my ($juncStart, $juncEnd) = ($fields[1]+$blockSizes[0], $fields[2]-$blockSizes[1]);
    $juncHash{left}{"$chr:$strand:$juncStart"} = '';
    $juncHash{right}{"$chr:$strand:$juncEnd"} = '';
}

while(<IN>){
    chomp;
    my ($chr, $spL, $spR, $strand) = $_ =~ /(.+):.+:(\d+)-(\d+):([+-])$/;
    my ($incSup, $excSup);
    if( ($anss == 5 && $strand eq '+') ||
        $anss == 3 && $strand eq '-' ){
        $excSup = exists $juncHash{left}{"$chr:$strand:$spL"} ? 'Y' : 'N';
        $incSup = exists $juncHash{left}{"$chr:$strand:$spR"} ? 'Y' : 'N';
    }else{
        $incSup = exists $juncHash{right}{"$chr:$strand:$spL"} ? 'Y' : 'N';
        $excSup = exists $juncHash{right}{"$chr:$strand:$spR"} ? 'Y' : 'N';
    }
    say join "\t", ($_, $incSup, $excSup);
}

