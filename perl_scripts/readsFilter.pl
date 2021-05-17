#!/usr/bin/perl -w

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use strict;
use Getopt::Long;
use File::Basename;

my ($minCov, $minIde, $maxRatio) = (67, 75, 0.9);
GetOptions(
            'c|coverage=s'  => \$minCov,
            'i|identity=s'  => \$minIde,
            'r|ratio=s'     => \$maxRatio,
            'h|help'        => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";

my %reads;
while(<IN>){
    chomp;
    my ($readName, $exonCount, $coverage, $identity) = (split "\t")[3, 9, 12, 13];
    if($coverage < $minCov && $identity < $minIde){
        print STDERR join "\t", ($_, "Low_CV_and_ID\n");
        next;
    }
    if($coverage < $minCov){
        print STDERR join "\t", ($_, "Low_CV\n");
        next;
    }
    if($identity < $minIde){
        print STDERR join "\t", ($_, "Low_ID\n");
        next;
    }
    my $score = $coverage * $identity;
    if(exists $reads{$readName}){
        my @reads = @{$reads{$readName}};
        if(@reads == 1){
            if($exonCount == 1 && $reads{$readName}->[0]->[2] == 1 ||
               $exonCount > 1 && $reads{$readName}->[0]->[2] > 1){
                if($score > $reads{$readName}->[0]->[0]){
                    unshift @{$reads{$readName}}, [$score, $_, $exonCount];
                }else{
                    push @{$reads{$readName}}, [$score, $_, $exonCount];
                }
            }elsif($exonCount > 1){
                print STDERR join "\t", ($reads{$readName}->[0]->[1], "Filtered_by_spliced\n");
                $reads{$readName}->[0] = [$score, $_, $exonCount];
            }else{
                print STDERR join "\t", ($_, "Filtered_by_spliced\n");
            }
        }else{ # @reads >=2
            if($exonCount == 1 && $reads{$readName}->[0]->[2] == 1 ||
               $exonCount > 1 && $reads{$readName}->[0]->[2] > 1){
                if($score >= $reads{$readName}->[0]->[0]){
                    print STDERR join "\t", ($reads{$readName}->[1]->[1], "Low_score\n");
                    $reads{$readName}->[1] = $reads{$readName}->[0];
                    $reads{$readName}->[0] = [$score, $_, $exonCount];
                }elsif($score > $reads{$readName}->[1]->[0]){
                    print STDERR join "\t", ($reads{$readName}->[1]->[1], "Low_score\n");
                    $reads{$readName}->[1] = [$score, $_, $exonCount];
                }else{
                    print STDERR join "\t", ($_, "Low_score\n");
                }
            }else{
                if($reads{$readName}->[0]->[2] == 1){
                    print STDERR join "\t", ($reads{$readName}->[0]->[1], "Filtered_by_spliced\n");
                    print STDERR join "\t", ($reads{$readName}->[1]->[1], "Filtered_by_spliced\n");
                    $reads{$readName}->[0] = [$score, $_, $exonCount];
                    pop @{$reads{$readName}};
                }else{
                    print STDERR join "\t", ($_, "Filtered_by_spliced\n");
                }
            }
        }
    }else{
        $reads{$readName} = [[$score, $_, $exonCount]];
    }
}

for my $key(keys %reads){
    my @hits = @{$reads{$key}};
    if(@hits == 1){
        print $hits[0]->[1], "\n";
    }else{
        my ($bestCoverage, $bestIdentity) = (split "\t", $hits[0]->[1])[12, 13];
        if($hits[0]->[2] > 1 && $hits[1]->[2] > 1 ||
           $hits[0]->[2] == 1 && $hits[1]->[2] == 1){
            if($hits[1]->[0] / $hits[0]->[0] <= $maxRatio){
                print $hits[0]->[1], "\n";
            }else{
                print STDERR join "\t", ($hits[0]->[1], "High_ratio_H\n");
            }
        }else{
            print $hits[0]->[1], "\n";
        }
        print STDERR join "\t", ($hits[1]->[1], "High_ratio_L\n");
    }
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.bed12+ >filtered.bed12+ 2>discarded.bed12+
    If INPUT isn't specified, input from STDIN
    Discarded reads are put into STDERR with the cause of discarding:
        Low_CV_and_ID: Low Coverage and Identity
        Low_CV: Low Coverage
        Low_ID: Low Identity
        Filtered_by_spliced: single-exon hit with the existence of anthor spliced hit
        High_ratio_H: the best one of the best and second hits with high score ratio
        High_ratio_L: the second one of the best and second hits with high score ratio
Options:
    -c --coverage   DOU Min coverage in percentage[67]
    -i --identity   DOU Min identity in percentage[75]
    -r --ratio      DOU Max ratio of the second hit to the best hit[0.9].
                        Reads with Coverage(second) * Identity(second) / Coverage(best) * Identity(best) > DOU will be discarded
    -h --help           Print this help information
HELP
    exit(-1);
}