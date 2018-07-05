### This is the main script to process go-plus.obo file
### By Haiming Tang, modified 06/08/2018

use strict;
#use warnings;
use 5.10.1;
use Data::Dumper;
use List::Util qw(shuffle);

$Data::Dumper::Sortkeys = 1;

my $goplusfile = "../rawData/go-plus.obo";


############# PARSE FROM GO-PLUS ###########################

my %termName; # $TermName{$id} = $name 
my %alter; # $Alter{$alter_id} = $id
my %strictRelation; # strict relationships: is_a, part_of, occurs_in 
my %otherRelation; # 2 level, first level is from term type to term type
my %relationCount; # inlcuding counts of relationships, from term type to type, counts
my %termType; # counts of each differnt type of terms

my $name;
my @goterm;
my $type;
my %skip;
my %isa; # for is_a of NCBITaxon only

my $maxround = 5;

open IN, "< $goplusfile" or die "cannot open $goplusfile\n";
while(<IN>){
    my $line = $_;
    chomp($line);
    
    if ($line =~ /\[Term\]/){
        if ($name){ 
            foreach my $goterm (@goterm){
                $termName{$goterm} = $name;
            }
        }
        @goterm = ();
    }
    
    elsif ($_ =~ /^id: (([A-Za-z]+):[0-9]+)/){
	$type = $2;
	$termType{$type} ++;
	push(@goterm,$1);
    }
    elsif ($_ =~ /name: (.*)/){
        $name = $1;
    }
    elsif ($_ =~ /alt_id: ([A-Za-z]:[0-9]+)/){
        $alter{$1} = $goterm[0];  # is the direction wrong??? CHECK LATER 
        push(@goterm,$1);
	$termType{$type}++;
    }
    elsif ($_ =~ /is_a:.* (([A-Za-z]+):[0-9]+)/){
        foreach my $goterm (@goterm){
	    my $relation = "is_a:$type:$2";
	    $strictRelation{$goterm}{$1} = $relation;
	    $relationCount{$relation} ++;
	    my $toterm = $1;
	    if ($goterm =~ /NCBITaxon/){
		$isa{$goterm} = $toterm;
	    }
        }
    }
    elsif ($line =~ /relationship: (\w+) (([A-Za-z]+):[0-9]+)/){
      my $relation = "$1:$type:$3";
      $relationCount{$relation}++;
            
      foreach my $goterm (@goterm){
	  if (($1 eq 'part_of') or ($1 eq 'occurs_in')){
	      $strictRelation{$goterm}{$2} = $relation;
	  }
	  else{
	      $otherRelation{$goterm}{$2} = $relation;
	  }
      }
    }
}
close IN;

my $termSta = "../statistics/termStatistics.txt";
unless (-e $termSta){
    open OUT, "> $termSta" or die "cannot open $termSta\n";
    print OUT Dumper(\%termType);
    close OUT;
}

my $nameSta = "../statistics/termName.txt";
unless (-e $nameSta){
    open OUT, "> $nameSta" or die "cannot open $nameSta\n";
    print OUT Dumper(\%termName);
    close OUT;
}

my $typeSta = "../statistics/relationTypeStatistics.txt";
unless (-e $typeSta){
    open OUT, "> $typeSta" or die "cannot open $typeSta\n";
    print OUT Dumper(\%relationCount);
    close OUT;
}

################################

$isa{'NCBITaxon:5476'} = 'NCBITaxon:4892';
$isa{'NCBITaxon:5833'} = 'NCBITaxon:5820';
$isa{'NCBITaxon:3055'} = 'NCBITaxon:33090';

my %goCons;
my %cons; # store the cons from onlyin and manual, not limited to GO terms 
foreach my $term1 (keys %strictRelation){
    next if $term1 =~ /NCBITaxon/;
    foreach my $term2 (keys %{$strictRelation{$term1}}){
	if ($strictRelation{$term1}{$term2} =~ /(\w+):\w+:NCBITaxon/){
	    $cons{$term1} .= $1."|".$term2.";";
	    my $re = $1;
	    if ($term1 =~ /GO:/){
		my $relation = 'Gain';
		if ($re =~ /never/){
		    $relation = 'Loss';
		}
		if ($term2 =~ /NCBITaxon/){
		    $goCons{$term1}{2} .= $relation."|".$term2.";";
		}
	    }
	}
    }
}
foreach my $term1 (keys %otherRelation){
    next if $term1 =~ /NCBITaxon/;
    foreach my $term2 (keys %{$otherRelation{$term1}}){
	if ($otherRelation{$term1}{$term2} =~ /(\w+):\w+:NCBITaxon/){
	    $cons{$term1} .= $1."|".$term2.";";
	}
	my $re = $1;
	if ($term1 =~ /GO:/){
	    my $relation = 'Gain';
	    if ($re =~ /never/){
		$relation = 'Loss';
	    }
	    if ($term2 =~ /NCBITaxon/){
		$goCons{$term1}{2} .= $relation."|".$term2.";";
	    }
	}
    }
}
## get the manual curation list
open MAN, "< ../rawData/manualCurationList_longform.csv" or die;
while(<MAN>){
    chomp;
    my @a = split(/,/);
    my $con = $a[1];
    $con =~ s/\([A-Za-z <>-]+\)//g;
    $cons{$a[0]} = "manual:".$con;
    $goCons{$a[0]}{3} = $con;
}
close MAN;

my $uberonConf = "../constraints/uberon_PO_FAO_CL.constraints.txt";

my %uberonCon; #store the constraints of uberon PO FAO CL terms
my %chebiCon;
my @rounts;
my %othercons;

unless (-e $uberonConf){
    &constructUBERON;

    open UBERON, "> $uberonConf" or die;
    print UBERON Dumper(\%uberonCon);
    close UBERON;
}
else{
    open UBERON, "< $uberonConf" or die;
    while (<UBERON>){
	if ($_ =~ /([\w:]+)\' => \'(.*)\'/){
	    my $go = $1; my $con = $2;
	    next unless $con;
	    unless ($con =~ /;$/){
		$con = $con.";";
	    }
	    $uberonCon{$go} = $con;
	}
    }
    close UBERON;
}


my $goterm;
open IN, "< ../constraints/chebi.go_ancestor.txt" or die;
while(<IN>){
  if ($_ =~ /(GO:[0-9]+)/){
    $goterm = $1;
  }
  my %single = eval($_);
  my $key = (keys %single)[0]; next unless $key;
  if ($key =~ /:/){
    $goCons{$goterm}{4} = ">Chebi|".$single{$key};
  }
}
close IN;

# propagate the uberon constraints to GO terms

foreach my $term1 (keys %strictRelation){
    next unless $term1 =~ /GO:/;
    foreach my $term2 (keys %{$strictRelation{$term1}}){
	if ($uberonCon{$term2}){
	    $goCons{$term1}{1} = $uberonCon{$term2};
	}
    }
}
foreach my $term1 (keys %otherRelation){
    next unless $term1 =~ /GO:/;
    foreach my $term2 (keys %{$otherRelation{$term1}}){
	if ($uberonCon{$term2}){
	    $goCons{$term1}{1} = $uberonCon{$term2};
	}
    }
}

my %combineCons;
my %combine_construct;

# first, a touch up is needed
my $line;
my $GOterm;
my $value;

foreach my $key (keys %termName){
  if ($key =~ /GO:/){
    $GOterm = $key;
    my $name = $termName{$key};
    my $goname = $termName{$key};

    if ($key =~ /GO:0009432/){
      $value = ">Gain|NCBITaxon:2;";
      $goCons{$key}{5} = $key;
      next;
    }

    if ($key =~ /GO:0060361/){
      $value = ">Gain|NCBITaxon:50557;";
      $goCons{$key}{5} = $key;
      next;
    }
    if ($key =~ /GO:0090632/){
      $value = ">Gain|NCBITaxon:40674;>Loss|NCBITaxon:9606";
      $goCons{$key}{5} = $key;
      next;
    }
    
    undef $value;
    
    if ($goname =~ /response to/){
      &process_res($name,$value);
    }
    elsif ($goname =~ /bacteri/){ # be careful of antibacterial
      &process_bac($name,$value);
    }
    elsif ($goname =~ /cellular bud/){
      &process_cel($name,$value);
    }
    elsif ($goname =~ /ukary/){ 
      &process_euk($name,$value);
    }
    elsif ($goname =~ /nuclear/){ # be careful of mononuclear
      &process_nuc($name,$value);
    }
    elsif ($goname =~ /(photosynthe)|(photosystem)/){
      &process_pho($name,$value);
    }
    elsif ($value =~ /NCBI/){
      $goCons{$key}{5} = $key;
    }
    else{
      my $l = $line;
      print "Missing constraint: $l"."\t$name\n";
    }
  }
}

open COM, "> ../constraints/goConsBeforeCombine.txt" or die;
print COM Dumper(\%goCons);
close COM;


&combine_construct;

open COM, "> ../constraints/goCons.txt" or die;
print COM Dumper(\%combineCons);
close COM;

open COM, "> ../constraints/goCons_constructlog.txt" or die;
print COM Dumper(\%combine_construct);
close COM;



sub process_cel{
  my ($name,$value) = @_;
  unless ($value =~ /;/){
    $value .= ";";
  }
  if ($value =~ /4751[^0-9]/){
    return;
  }
  else{
    print "CHECK process_cel !! $name ; $value \n";

    $value = ">Gain|NCBITaxon:4751;";
    $goCons{$GOterm}{5} = $value;
  }
}

sub process_pho{
  my ($name,$value) = @_;

  $value = ">Gain|NCBITaxon:33090;NCBITaxon:1117";
  $goCons{$GOterm}{5} = $value;
  
}

sub process_euk{
  my ($name,$value) = @_;
  unless ($value =~ /;/){
    $value .= ";";
  }
  if ($value =~ /:2759[^0-9]/){
    return;
  }
  elsif ($value =~ /(:1[^0-9])|(:131567[^0-9])/){
    print "check process_euk LUCA!! $name ; $value \n";
    $value = ">Gain|NCBITaxon:2759;";
    $goCons{$GOterm}{5} = $value;
  }
  elsif ($value){
    print "check process_euk OTHER CONS!! $name ; $value \n";
    $value = ">Gain|NCBITaxon:2759;";
    $goCons{$GOterm}{5} = $value;    
  }
  else{
    print "check process_euk NO CONS!! $name ; $value \n";
    $value = ">Gain|NCBITaxon:2759;";
    $goCons{$GOterm}{5} = $value;
  }
}

sub process_nuc{
  my ($name,$value) = @_;
  unless ($value =~ /;/){
    $value .= ";";
  }
  if ($value =~ /:2759[^0-9]/){
    return;
  }
  else{

    if ($value =~ /(:1[^0-9])|(:131567[^0-9])/){
      print "CHECK process_nuc !! LUCA $name ; $value \n";
      $value = ">Gain|NCBITaxon:2759;";
      $goCons{$GOterm}{5} = $value;
    }
    else{
      print "CHECK process_nuc!! OTHER $name ; $value \n";
      #$value = ">Gain|NCBITaxon:2759;";
      # $goCons{$GOterm}{5} = $value;
    }
  }
}

sub process_res{
  my ($name,$value) = @_;
  unless ($value =~ /;/){
    $value .= ";";
  }

  if (($value =~ /NCBITaxon:1[^0-9]/ ) or ($value =~ /NCBITaxon:131567[^0-9]/)){
    
    print "response to!!  $name".  $line;
    return;
  }
  elsif ($value =~ /NCBITaxon/){
    $goCons{$GOterm}{5} = $value;
  }
  else{
    print "repsonse to process_res!! no constraint! $name , $line";
  }
}

sub process_bac{
  my ($name,$value) = @_;
  unless ($value =~ /;/){
    $value .= ";";
  }
  if ($name =~ /antibacter/){
    return;
  }
  elsif ($name =~ /archaeal/){
    return;
  }
  elsif (($name =~ /type/) or ($name =~ /bacterial /) or ($name =~ /bacterium/) or ($name =~ /bacterio/)){
    unless ($value){
      $value = ">Gain|NCBITaxon:2;";
      $goCons{$GOterm}{5} = $value;

    }
    elsif ($value =~ /NCBITaxon:2[^0-9]/){
      return;
    }
    elsif ($value =~ /NCBITaxon:1[^0-9]/){
      $value = ">Gain|NCBITaxon:2;";
      $goCons{$GOterm}{5} = $value;
    }    
    elsif ($value =~ /NCBITaxon:131567[^0-9]/){
      $value = ">Gain|NCBITaxon:2;";
      $goCons{$GOterm}{5} = $value;
    }
    else{
      print "check for $name: $line";
      $goCons{$GOterm}{5} = $value;
    }
  }
  else{
    print "process_bac CHECK!! other condition: ".$name."\t".$line;
  }
}

sub  constructUBERON{
    # the subroutine to get constraints for uberon, po, fao, cl terms
    foreach my $term (keys %termName){
	if ($term =~ /PO:/) {
	    $uberonCon{$term} = 'NCBITaxon:33090'; # plants get Viridiplantae   
	}
	elsif ($term =~ /FAO:/){
	    $uberonCon{$term} = 'NCBITaxon:4751'; # fungi get fungi
	}
	elsif ($term =~ /UBERON:/){
	    next if ($term =~ /0000061/);
	    next if ($term =~ /0000465/);
	    next if ($term =~ /0001062/);
	    next if ($term =~ /0000468/);
	    
	    $uberonCon{$term} = 'NCBITaxon:33208'; # animals get metazoa
	}
	elsif ($term =~ /CL:/){
	    $uberonCon{$term} = '';
	}
	# a few exceptions exist, high level terms
    }

    my %previousCon;

    my $c =0;
    while(1){	
	&uberonCycle;
	my $diff=  &compareCons(\%previousCon,\%uberonCon);
	last if $diff < 1;
	%previousCon = %uberonCon;
	$c++;
	last if $c > $maxround;
    }
}


sub compareCons{
    my $ref1 = shift;
    my $ref2 = shift;
    
    # check if constraints are the same 

    my @set1 = keys %{$ref1};
    my @set2 = keys %{$ref2};

    unless (&identical(\@set1,\@set2)){
	return 1;
    }

    foreach my $term (@set1){
	my $l1 = $ref1->{$term};
	my $l2 = $ref2->{$term};

	my @l1 ; my @l2;
	while($l1 =~ /(NCBITaxon:[0-9]+)/g){
	    push(@l1,$1);
	}
	while($l2 =~ /(NCBITaxon:[0-9]+)/g){
	    push(@l2,$1);
	}
	unless (&identical(\@l1,\@l2)){
	    return 1;
	}
    }
    return 2;
}


sub uberonCycle{
# simple uberon Cycle, where     
    foreach my $uterm (shuffle keys %uberonCon){
	my $hashref = &getRoute($uterm);
	my @taxons;
#	my @set;
	while($uberonCon{$uterm} =~ /(NCBITaxon:[0-9]+)/g){
	    push(@taxons,$1);
	}
	foreach my $level (keys %{$hashref}){
	    foreach my $term (keys %{$hashref->{$level}}){
		if ($cons{$term}){
		    my @tmp = split(/;/,$cons{$term});
		    foreach my $tmp(@tmp){
			my @next = split(/\|/,$tmp);
			next if $next[0] =~ /never/;
			push(@taxons,$next[1]);
		    }
		}
	    }
	}
	my $taxon = &mostTaxon(\@taxons);
	if ($taxon){	 	 
	    $uberonCon{$uterm} = "Gain|".$taxon;
#	    $cons{$uterm} .= "uberonCycle:".$taxon ;
	}
    }
}

sub getRoute{
    my $term = shift;
    my %route;
    my $n = 1;
    my @temp; push(@temp,$term);
    my $curSTerm = \@temp;
    %skip = ();
    while(1){		
	my $terms = &nextlevel($curSTerm);

	my @temp;
	foreach my $iterm (@$terms){
	    $term =~ /(\w+):/; my $type1 = $1;
	    $iterm =~ /(\w+):/; my $type2 = $1;
	    if ($type1 eq $type2){
		unless ($cons{$iterm}){
		    unless ($skip{$iterm}){
			push(@temp,$iterm);
		    }
		}
	    }
	    $route{$n}{$iterm} =1;
	    $skip{$iterm} =1;
	}
	$n++;
	if (@temp){
	    $curSTerm = \@temp;
	}	
	else{
	    last;
	}
    }
    foreach my $key (keys %{$otherRelation{$term}}){
	$route{'other'}{$key} =1;
    }
    
    return \%route;
}

sub nextlevel{
    my $aref = shift;
    my %nextlevel;
    
    foreach my $key (@$aref){
	foreach my $ref (keys %{$strictRelation{$key}}){
	    $nextlevel{$ref} =1;
	}
    }
    my @keys = keys %nextlevel;
    return \@keys;
}


sub combine_construct{

 #   my %previousCon;

    my $c =0;
    while(1){
	print "run combineCycle $c\n";
	&combineCycle($c);
#	my $diff=  &compareCons(\%previousCon,\%combineCons);

#	print "differnece is $diff\n";
#	last if $diff == 2;
#	%previousCon = %combineCons;
	$c++;
	last if $c > $maxround;
    }
}

sub combineCycle{
    my $cycle = shift;
    
    foreach my $goterm (keys %termName){
	next unless $goterm =~ /GO:/;

	my $hashref = &getRoute($goterm);
        
	my %gaintaxons;
	my %losstaxons;
	my %chebitaxons;

	foreach my $n (1..5){
	    my $con = $goCons{$goterm}{$n}; next unless $con;
	    while ($con =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
		my $type = $1;
		my $string = $2;
		while ($string =~ /(NCBITaxon:[0-9]+)/g){
		    if ($n ==4){
			$chebitaxons{$1}{$n}{$goterm} =1;
		    }
		    elsif ($type eq 'Gain'){
			$gaintaxons{$1}{$n}{$goterm} =1;
		    }
		    else{
			$losstaxons{$1}{$n}{$goterm} =1;
		    }
		}
	    }
	}
	
        foreach my $level (sort keys %{$hashref}){
            foreach my $term (keys %{$hashref->{$level}}){
		foreach my $n (1..5){
		    my $con = $goCons{$term}{$n}; next unless $con;
		    while ($con =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
			my $type = $1;
			my $string = $2;
			while ($string =~ /(NCBITaxon:[0-9]+)/g){
			    if ($n ==4){
				$chebitaxons{$1}{$n}{$term} =1;
			    }
			    elsif ($type eq 'Gain'){
				$gaintaxons{$1}{$n}{$term} =1;
			    }
			    else{
				$losstaxons{$1}{$n}{$term} =1;
			    }
			}
		    }
		}
	    }
	}
	my @gain = keys %gaintaxons;
	if (@gain){
	    my $gain = &mostTaxon(\@gain);
	    my $source;
	    while($gain =~ /(NCBITaxon:[0-9]+)/g){
		foreach my $type (keys %{$gaintaxons{$1}}){
		    foreach my $term (keys %{$gaintaxons{$1}{$type}}){
			$source .= "$type:$term;";
		    }
		}
	    }
	    $combine_construct{$goterm}{$cycle}{'gain'} = $source;
	    $goCons{$goterm}{5} .= "Gain|".$gain.";";
	    $combineCons{$goterm}{'Gain'} = $gain;
	}
	else{
	    my @chebi = keys %chebitaxons;
	    if (@chebi){
		my $gain = &mostTaxon(\@chebi);
		my $source;
		while($gain =~ /(NCBITaxon:[0-9]+)/g){
		    foreach my $type (keys %{$chebitaxons{$1}}){
			foreach my $term (keys %{$chebitaxons{$1}{$type}}){
			    $source .= "$type:$term;";
			}
		    }
		}
		$combine_construct{$goterm}{$cycle}{'gain'} = $source;
		$goCons{$goterm}{5} .= "Gain|".$gain.";";
		$combineCons{$goterm}{'Chebi'} = $gain;
	    }
	}
	my @loss = keys %losstaxons;
	if (@loss){
	    # basically, a loss contraint should not be older than the gain constraints
	    my $gaininfo = $combineCons{$goterm}{'Gain'}.";".$combineCons{$goterm}{'Chebi'};
	    
	    my $loss = &mostTaxonLoss($gaininfo,\@loss);
	    my $source;
	    while($loss =~ /(NCBITaxon:[0-9]+)/g){
		foreach my $type (keys %{$losstaxons{$1}}){
		    foreach my $term (keys %{$losstaxons{$1}{$type}}){
			$source .= "$type:$term;";
		    }
		}
	    }
	    $combine_construct{$goterm}{$cycle}{'loss'} = $source;
	    $goCons{$goterm}{5} .= "Loss|".$loss.";";
	    $combineCons{$goterm}{'Loss'} = $loss;
	}
    
    }
}


sub mostTaxon{
    # return the youngest ceancestors from a group
    
  my $aref = shift;
  my @arr = @$aref;
  my %taxons;

  unless ($aref){
      return;
  }
  foreach my $line (@arr){
      while ($line =~ /(NCBITaxon:[0-9]+)/g){
						 
      $taxons{$1} =1;
    }
  }
  my @taxons = keys %taxons;
  my $size = @taxons;
  if ($size ==1){
    return $taxons[0];
  }
  elsif ($size < 1){
      return ;
  }
  else{
    my $rest;
    %skip = ();
    foreach my $i (0..$size-2){
      foreach my $j ($i+1..$size-1){
        if ((exists $skip{$i}) or (exists $skip{$j})){
          next;
        }
        else{
          my $skip_n = &Taxoncompare($taxons[$i],$taxons[$j]);
          if ($skip_n ==1){
            $skip{$j} =1;
#            say "skipping $taxons[$j], mom of $taxons[$i]";
          }
          elsif ($skip_n ==2){
            $skip{$i} =1;
 #           say "skipping $taxons[$i], mom of $taxons[$j]";
          }
        }
      }
    }
    foreach my $i (0..$size-1){
      next if (exists $skip{$i});
      my $taxon = $taxons[$i];
      $rest .= $taxon.";";
    }
    return $rest;
  }
}

sub mostTaxonLoss{
    # return the oldest ceancestors from a group
    
    my $gaininfo = shift;
    my %gainI;

    while($gaininfo =~ /(NCBITaxon:[0-9]+)/g){
	$gainI{$1} =1;
    }
    
    my $aref = shift;
    my @arr = @$aref;
    my %taxons;

    unless ($aref){
	return;
    }
    foreach my $line (@arr){
	while ($line =~ /(NCBITaxon:[0-9]+)/g){	    
	    $taxons{$1} =1;
	}
    }
    
    my @taxons;

    foreach my $taxonL (keys %taxons){
	my $I =1;
	foreach my $taxonG (keys %gainI){	    
	    my $skip_n = &Taxoncompare($taxonG,$taxonL);
	    if ($skip_n == 1 or $skip_n == -1){
		# $taxonL is mom of $taxonG; or $taxonL is equal to $taxonG
		$I = 2; last;		
	    }
	}
	if ($I ==1){
	    push(@taxons,$taxonL);
	}
    }

    my $size = @taxons;
    
    
    if ($size ==1){
	my $taxonL = $taxons[0];
	return $taxons[0];
    }
    elsif ($size < 1){
	return ;
    }
  else{
    my $rest;
    %skip = ();
    foreach my $i (0..$size-2){
      foreach my $j ($i+1..$size-1){
        if ((exists $skip{$i}) or (exists $skip{$j})){
          next;
        }
        else{
          my $skip_n = &Taxoncompare($taxons[$i],$taxons[$j]);
          if ($skip_n ==1){
            $skip{$i} =1;
#            say "skipping $taxons[$i], mom of $taxons[$j]";
          }
          elsif ($skip_n ==2){
            $skip{$j} =1;
 #           say "skipping $taxons[$j], mom of $taxons[$i]";
          }
        }
      }
    }
    foreach my $i (0..$size-1){
      next if (exists $skip{$i});
      my $taxon = $taxons[$i];
      $rest .= $taxon.";";
    }
    return $rest;
  }
}


sub Taxoncompare{
  my $taxon1 = shift;
  my $taxon2 = shift;

  my $mum1;
  my $mum2;

  my $taxon = $taxon1;
  my $current;

  if ($taxon1 eq $taxon2){
    return -1;
  }

  while(1){
    if ($current){
      $taxon = $current;
    }
    if (exists $isa{$taxon}){
      $current = $isa{$taxon};
      $current =~ /(NCBITaxon:[0-9]+)/;
      $current = $1;
      $mum1 .= $current.";";
    }
    else{
      last;
    }
  }

   $taxon = $taxon2;
   $current = '';
  
  while(1){
    if ($current){
      $taxon = $current;
    }
    
    if (exists $isa{$taxon}){
      $current = $isa{$taxon};
      $current =~ /(NCBITaxon:[0-9]+)/;
      $current = $1;
      $mum2 .= $current.";";
    }
    else{
      last;
    }
  }
  
  if ($mum1 =~ /$taxon2[^0-9]/){
    return 1;
  }
  elsif ($mum2 =~ /$taxon1[^0-9]/){
    return 2;
  }
  else{
#    print STDERR "$taxon1 mom: $mum1\n";
 #   print STDERR "$taxon2 mom: $mum2\n\n";
    return 3;
  }
}

sub identical {
    my( $left, $right ) = @_;
    return 0 if scalar @$left != scalar @$right;
    my %hash;
    @hash{ @$left, @$right } = ();
    return scalar keys %hash == scalar @$left;
}
