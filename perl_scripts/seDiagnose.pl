#!/usr/bin/env perl

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
    -h --help           Print this help information
HELP
    exit(-1);
}

my $juncFile;
GetOptions(
            'j|junc=s'  => \$juncFile,
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
    $juncHash{"$chr:$strand:$juncStart-$juncEnd"} = '';
}

while(<IN>){
    chomp;
    my ($chr, $spL, $ses, $spR, $strand) = $_ =~ /(.+):.+:(\d+)@(.+)@(\d+):([+-])$/;
    my @ses = split ';', $ses;
    my (@seStarts, @seEnds);
    for my $se(@ses){
        my ($seStart, $seEnd) = split '-', $se;
        push @seStarts, $seStart;
        push @seEnds, $seEnd;
    }
    my $spLSup = exists $juncHash{"$chr:$strand:$spL-$seStarts[0]"} ? 'Y' : 'N';
    my $sesSup = 0;
    for(my $i = 1; $i <= $#ses; $i++){
        $sesSup++ if(exists $juncHash{"$chr:$strand:" . $seEnds[$i-1] . "-$seStarts[$i]"});
    }
    my $spRSup = exists $juncHash{"$chr:$strand:$seEnds[-1]-$spR"} ? 'Y' : 'N';
    my $skipSup = exists $juncHash{"$chr:$strand:$spL-$spR"} ? 'Y' : 'N';
    say join "\t", ($_, $spLSup, $sesSup, $#ses, $spRSup, $skipSup);
}

