#!/usr/bin/perl -w
use strict;
use warnings;

my @runs = @ARGV;

# Vraies mutations
my %truemutations;

$truemutations{"540T"} = 1; # OK
$truemutations{"572Q"} = 1; # OK
$truemutations{"573N"} = 1; # OK
$truemutations{"665N"} = 1; # OK
$truemutations{"731V"} = 1; # OK
$truemutations{"780S"} = 1; # OK
$truemutations{"751F"} = 1; # OK
$truemutations{"623D"} = 1; # OK
$truemutations{"505L"} = 1; # OK
$truemutations{"749T"} = 1; # NO
$truemutations{"761A"} = 1; # OK
$truemutations{"839K"} = 1; # OK

my $nbtrue = keys %truemutations;

my @analyses = ("emergence", "correlation", "condor");

my %tp;
my %fp;
my %fn;
my %tn;
my %recall;
my %precision;
my %f1;
my %type1;

print "Analysis\tReplicate\tTP\tFP\tFN\tTN\tType1\tRecall\tPrecision\tF1Score\n";

sub Type1{
    my ($tp,$fp,$fn,$tn) = @_;

    my $type1= "N/A";
    
    if($tn+$fp > 0){
	$type1 = $fp/($tn+$fp);
    }

    return $type1;
}

sub Precision{
    my ($tp,$fp,$fn,$tn) = @_;

    my $precision = "N/A";

    if($tp+$fp > 0){
	$precision = $tp/($tp+$fp);
    }

    return $precision;
}

sub Recall{
    my ($tp,$fp,$fn,$tn) = @_;

    my $recall = "N/A";

    if($tp+$fn > 0){
	$recall = $tp/($tp+$fn);
    }

    return $recall;
}

sub F1Score{
    my ($precision, $recall) = @_;

    my $f1 = "N/A";

    if($precision ne "N/A" && $recall ne "N/A" && ($precision+$recall) > 0){
	$f1 = 2*($precision*$recall)/($precision+$recall);
    }

    return $f1;
}

my $nbrun=0;
for my $run (@runs) {
    my %seen;
    my $tested=$run."/tested_results.tsv";

    my %temp_tp;
    my %temp_fp;
    my %temp_fn;
    my %temp_tn;

    for my $a (@analyses) {
	$temp_tp{$a} = 0;
	$temp_fp{$a} = 0;
	$temp_fn{$a} = 0;
	$temp_tn{$a} = 0;
    }
    
    open(IN,$tested);
    <IN>;
    while(<IN>){
	chomp;
	my @cols = split(/\t/);
	my $pastml_root = $cols[0];
	my $consensus_root  = $cols[1];
	my $position = $cols[2];
	my $mut = $cols[3];
	my $max_anc = $cols[4];
	my $ref_EEM = $cols[5];
	my $nbseq = $cols[6];
	my $evol_rate = $cols[7];
	my $genetic_distance = $cols[8];
	my $substitution_rate = $cols[9];
	my $findability = $cols[10];
	my $type_substitution = $cols[11];
	my $details = $cols[12];
	my $loss = $cols[13];
	my $loss_details = $cols[14];
	my $max_simu = $cols[15];
	my $variance = $cols[16];
	my $mean = $cols[17];
	my $pvalue_raw = $cols[18];
	my $adjust_pvalue = $cols[19];
	my $adjust_pvalue_fdr = $cols[20];
	my $detected_EEM = $cols[21];
	my $posmut = $cols[22];
	my $log_dep = $cols[23];
	my $log_indep = $cols[24];
	my $BF = $cols[25];
	my $correlation = $cols[26];

	my $mutStr = "".$position.$mut;

	$seen{$mutStr}=1;
	
	if(defined $truemutations{$mutStr}){
	    if($adjust_pvalue <= 0.1){
		$temp_tp{emergence}+=1;
		#print "TP emerg: $mutStr\n";
	    }else{
		$temp_fn{emergence}+=1;
		#print "FN emerg: $mutStr\n";
	    }
	    if($BF >= 2){
		$temp_tp{correlation}+=1;
		#print "TP corr: $mutStr\n";		
	    }else{
		$temp_fn{correlation}+=1;
		#print "FN corr: $mutStr\n";
	    }
	    if($adjust_pvalue <= 0.1 && $BF >= 2){
		$temp_tp{condor}+=1;
		#print "TP cond: $mutStr\n";
	    }else{
		$temp_fn{condor}+=1;
		#print "FN cond: $mutStr\n";
	    }
	}else{
	    if($adjust_pvalue<=0.1){
		$temp_fp{emergence}+=1;
		#print "FP emerg: $mutStr\n";
	    }else{
		$temp_tn{emergence}+=1;
	    }
	    if($BF >= 2){
		$temp_fp{correlation}+=1;
		#print "FP corr: $mutStr\n";
	    }else{
		$temp_tn{correlation}+=1;
	    }
	    if($adjust_pvalue <= 0.1 && $BF >= 2){
		$temp_fp{condor}+=1;
		#print "FP cond: $mutStr\n";
	    }else{
		$temp_tn{condor}+=1;
	    }
	}
    }

    for my $s (keys(%truemutations)){
	if(!defined $seen{$s}){
	    print STDERR "Run $run : Not defined : $s\n";
	}
    }
   

    for my $a (@analyses){
	$tp{$a} += $temp_tp{$a};
	$fp{$a} += $temp_fp{$a};
	$fn{$a} += $temp_fn{$a};
	$tn{$a} += $temp_tn{$a};

	my $temp_precision;
	my $temp_f1;
	my $temp_type1;
	my $temp_recall;

	$temp_recall = Recall($temp_tp{$a}, $temp_fp{$a}, $temp_fn{$a}, $temp_tn{$a});
	$temp_type1 = Type1($temp_tp{$a}, $temp_fp{$a}, $temp_fn{$a}, $temp_tn{$a});
	$temp_precision = Precision($temp_tp{$a}, $temp_fp{$a}, $temp_fn{$a}, $temp_tn{$a});
	$temp_f1 = F1Score($temp_precision, $temp_recall);
	
	print $a."\t".$nbrun."\t".$temp_tp{$a}."\t".$temp_fp{$a}."\t".$temp_fn{$a}."\t";
	print $temp_tn{$a}."\t".$temp_type1."\t".$temp_recall."\t".$temp_precision."\t".$temp_f1."\n";
    }
    $nbrun+=1;
    close(IN);
}

for my $a (@analyses){
 
    my $recall = Recall($tp{$a}/$nbrun, $fp{$a}/$nbrun, $fn{$a}/$nbrun, $tn{$a}/$nbrun);
    my $type1 = Type1($tp{$a}/$nbrun, $fp{$a}/$nbrun, $fn{$a}/$nbrun, $tn{$a}/$nbrun);
    my $precision = Precision($tp{$a}/$nbrun, $fp{$a}/$nbrun, $fn{$a}/$nbrun, $tn{$a}/$nbrun);
    my $f1 = F1Score($precision, $recall);

    print $a."\tAVG\t";
    print $tp{$a}/$nbrun."\t";
    print $fp{$a}/$nbrun."\t";
    print $fn{$a}/$nbrun."\t";
    print $tn{$a}/$nbrun."\t";
    print $type1."\t";
    print $recall."\t";
    print $precision."\t";
    print $f1."\n";
}
