### This is the main script to process go-plus.obo file
### By Haiming Tang, modified 06/08/2018

use strict;
#use warnings;
use 5.10.1;
use Data::Dumper;
use List::Util qw(shuffle);

$Data::Dumper::Sortkeys = 1;

my $goplusfile = "../rawData/go-plus.obo";

my %goType;

############# PARSE FROM GO-PLUS ###########################
my $goType = "../statistics/go_Type.txt" ;


my %termType; 
my $name;
my $term;

if (-e $goType){
    open IN, "< $goType" or die;
    while(<IN>){
	if ($_ =~ /\'(GO:[0-9]+)\' => \'(\w+)\'/){
	    $termType{$1} = $2;
	}
    }
    close IN;
}
else{
    open IN, "< $goplusfile" or die "cannot open $goplusfile\n";
    while(<IN>){
	my $line = $_;
	chomp($line);
    	
	if ($_ =~ /^id: (([A-Za-z]+):[0-9]+)/){
	    $term = $1;
	}
	elsif ($_ =~ /namespace: (.*)/){
	    $name = $1;
	    if ($term =~ /GO:/){
		$termType{$term} = $name;
	    }
	}
    }
    close IN;
    
    open OUT, "> $goType" or die;
    print OUT Dumper(\%termType);
    close OUT;
}


my $cons = "../constraints/goCons.txt" ;
my %goCons;

open IN, "< $cons" or die $!;
while(<IN>){
    if ($_ =~ /\'(GO:[0-9]+)\'/){
	$term = $1;
    }
    elsif ($_ =~ /\'(\w+)\' => \'([\w+:;]+)\'/){
	$goCons{$term}{$1} = $2;
    }
}
close IN;

my %counts;

foreach my $goterm (keys %termType){
    my $type = $termType{$goterm};
    $counts{$type}{'total'}++;
    $counts{'total'}++;
    
    my $gain = $goCons{$goterm}{'Gain'};
    my $loss = $goCons{$goterm}{'Loss'};

    if ($gain.";".$loss =~ /NCBITaxon/){
	$counts{$type}{'withTaxon'}++;
	$counts{'withTaxon'}++;
    }
    if ($loss =~ /NCBITaxon/){
	$counts{$type}{'withLoss'}++;
	$counts{'withLoss'}++;
    }
    my $n =0; my $l = 0;

    while ($gain =~ /(NCBITaxon:[0-9]+)/g){
	$n++;
	if ($1 ne 'NCBITaxon:1'){
	    $l =1;
	}
	last if $n > 1;
    }
    if ($l ==1){
	$counts{$type}{'noneRoot'}++;
	$counts{'noneRoot'}++;
    }
    if ($n >1){
	$counts{$type}{'more1'}++;
	$counts{'more1'}++;
    }

}

print ",Total number of GO terms, GO terms with at least 1 none-root constraint, GO terms with Loss taxonomic ranges, GO terms with more than 1 Gain taxonomic ranges\n";

foreach my $type (keys %counts){
    if ($type =~ /_/){
	print "$type,";
	print $counts{$type}{'total'}.",";
#	print "withTaxon\t".$counts{$type}{'withTaxon'}."\t".($counts{$type}{'withTaxon'}/$counts{$type}{'total'})."\n";

	print $counts{$type}{'noneRoot'}." (".sprintf("%.2f%%",100*($counts{$type}{'noneRoot'}/$counts{$type}{'total'}))."),";
	
	print $counts{$type}{'withLoss'}." (".sprintf("%.2f%%",100*($counts{$type}{'withLoss'}/$counts{$type}{'total'}))."),";

	print $counts{$type}{'more1'}." (".sprintf("%.2f%%",100*($counts{$type}{'more1'}/$counts{$type}{'total'})).")\n";
	
    }
}
print "Total,";
print $counts{'total'}.",";

print $counts{'noneRoot'}." (".sprintf("%.2f%%",100*($counts{'noneRoot'}/$counts{'total'}))."),";
print $counts{'withLoss'}." (".sprintf("%.2f%%",100*($counts{'withLoss'}/$counts{'total'}))."),";
print $counts{'more1'}." (".sprintf("%.2f%%",100*($counts{'more1'}/$counts{'total'})).")\n";
    
	

