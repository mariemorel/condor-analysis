#!/usr/bin/perl -w
use strict;

my $simueem=$ARGV[0];
my $condoreem=$ARGV[1];

my %compeems;

open(IN,$simueem);
<IN>;
while(<IN>){
    chomp;
    my @cols = split(/\t/);

    my $site  = $cols[1];
    my $parent= $cols[2];
    my $child = $cols[3];
    my $eems  = $cols[4];

    $compeems{($site+1)."-".$child}{"simu"} = $eems;
}
close(IN);

open(IN,$condoreem);
<IN>;
while(<IN>){
    chomp;
    my @cols = split(/\t/);
    my $site  = $cols[2];
    my $child = $cols[3];
    my $eems  = $cols[5];

    $compeems{$site."-".$child}{"condor"} = $eems;
}
close(IN);

print "Site\tMutation\tSimuEEMs\tCondorEEMs\n";
for my $k (keys(%compeems)){
    my @cols = split(/-/,$k);
    my $simu=0;
    if(defined $compeems{$k}{"simu"}){
	$simu=$compeems{$k}{"simu"};
    }
    my $condor=0;
    if(defined $compeems{$k}{"condor"}){
	$condor=$compeems{$k}{"condor"};
    }
    print $cols[0]."\t".$cols[1]."\t".$simu."\t".$condor."\n";
}
