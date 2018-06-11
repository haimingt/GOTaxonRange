## reconstruct chebi

use HAIMING::Ancestor;
use 5.10.1;
use Data::Dumper;
#use Switch;
## this script is for parsing of go-plus.obo file

$Data::Dumper::Sortkeys = 1;

my $golimit = shift;
my $infile = "go-plus.obo";

my %name; # store name of GO terms as well as others
my %isa; # store isa information of go terms, point go term to mom

my %onlyin; # store constraint info for go terms, specially only in taxonid
my %chebi; # store chebi info for go terms
my %uberon; # store uberon/pato... info for go terms

my $name;
my @goterm;

my %type;
my %type2;

my %regulates;

### basic idea use onlyin chebi uberon as constraint

open IN, "< $infile" or die;
while(<IN>){
    my $line = $_;
    
    if ($line =~ /\[Term\]/){
        if ($name){
            foreach my $goterm (@goterm){
                $name{$goterm} = $name;
            }
        }

        @goterm = ();
    }
    if ($_ =~ /^id: (([A-Za-z]+):[0-9]+)/){

	$type2{$2} =1;
        $goterms{$1} = 1;
	push(@goterm,$1);
    }


    if ($_ =~ /name: (.*)/){
        $name = $1;
    }


    if ($_ =~ /alt_id: ([A-Za-z]:[0-9]+)/){
        $alt{$goterm[0]} = $1;
        push(@goterm,$1);
    }


    if ($_ =~ /is_a:.* ([A-Za-z]+:[0-9]+)/){
        foreach my $goterm (@goterm){
	  $isa{$goterm} .= $1.";";
	  $regulates{$goterm} .= $1.";";
        }
    }
#    if ($_ =~ /intersection_of:.*([A-Za-z]+:[0-9]+)/){
 #       foreach my $goterm (@goterm){
  #          $isa{$goterm} .= $1.";";
   #         #$iso{$goterm} .= $1.";";
    #    }
    #}
    if ($line =~ /relationship:(.*) ([A-Za-z]+:[0-9]+)/){
      my $relation = $1;
      $relation{$relation} = 1;
      my $c = $2;
      foreach my $goterm (@goterm){
	
	if ($goterm ne $c){
	  $isa{$goterm} .= $c.";";
	}

	if ($goterm =~ /GO:/){
	  if ($relation =~ /(regulate)|(part of)/) {	  
	    $regulates{$goterm} .= $c.";";
	    
	  }
	  if ($c =~ /CHEBI/){
	    $chebi{$goterm} .= $c.";";
	  }
	  elsif ($c =~ /NCBITaxon/){
	    
	    $onlyin{$goterm}{$relation} .= $c.";";
	  }
	      
	  else{
	    next if $c =~ /GO:/;
	    $uberon{$goterm} .= $c.";";
	    if ($c =~ /([A-Za-z]+:)/){
	      $type{$1} = 1;
	    }
	  }		
	}
	else{
	  $othercons{$goterm} .= $c.";";
	}
      }	
    }
}
close IN;


$isa{'NCBITaxon:5476'} = 'NCBITaxon:4892';
$isa{'NCBITaxon:5833'} = 'NCBITaxon:5820';
$isa{'NCBITaxon:3055'} = 'NCBITaxon:33090';


=pod

open OTHER, "> otherCons.txt";
open TYPE, "> type.txt";
open NAME, "> name.txt";
open ISA, "> isa.txt";
open CHEBI, "> chebi.txt";
open ONLYIN, "> onlyin.txt";
open UBERON, "> uberon.txt";
open RELA, "> relation.txt";

print OTHER Dumper(\%othercons);
print RELA Dumper(\%relation);
say "Ptype";
print TYPE Dumper(\%type);
say "Ptype2";
print TYPE Dumper(\%type2);

say "Pname";
print NAME Dumper(\%name);

say "Pisa";
print ISA Dumper(\%isa);

say "Pchebi";
print CHEBI Dumper(\%chebi);

say "Ponlyin";
print ONLYIN Dumper(\%onlyin);

say "Puberon";
print UBERON Dumper(\%uberon);
close RELA;
close OTHER;
close TYPE;
close NAME;
close ISA;
close CHEBI;
close ONLYIN;
close UBERON;

=cut

#&chebiConConstruct;
#exit;

my @routs;

my %uberonCon;
my %skip;
my %po_fao_con;

my %pthr_NCBI;

open PN, "< pthr_NCBI.txt" or die;
while(<PN>){
  chomp;
  my @array = split(/\t| +/);
  my $size = @array;
  $pthr_NCBI{lc($array[0])} = "NCBITaxon:".$array[$size-1];

}
close PN;

$pthr_NCBI{'gain'} = 'Gain';
$pthr_NCBI{'loss'} = '>Loss';

print Dumper(\%pthr_NCBI);

my %goCons;

my %chebiCon; # constraint for chebi terms from annotation

&chebiConConstruct;
&constructChebi;

open CHB, "> chebifixed.txt" or die;
print CHB Dumper(\%goCons);
close CHB;
exit;

&constructUBERON;
open UBER, "> uberonConstraint.txt" or die;
print UBER Dumper(\%uberonCon);
close UBER;

open SUM, "< Haiming_summary.csv" or die;
while(<SUM>){
  chomp;
  my @array = split(/,/);
  my $go = $array[1];
  my $info = $array[2];
  $info = lc $info;

  my $out;

  while( $info =~ /([A-Za-z-]+)/g){
    my $target = $pthr_NCBI{$1};

    print $1."\t";


    if (($target =~ /Gain/) or ($target =~/Loss/)){
      $out .= $target."\|";
    }
    else{
      $out .= $target.";";
    }
  }
  print "\n";
  $goCons{$go}{3} = $out;
}
close SUM;

print Dumper(\%goCons);

&constructGO;

open GO, "> goConstraint.txt" or die;
print GO Dumper(\%goCons);
close GO;


my %result;
my %types;

&compareCons;

open GO, "> compareConstraintTypes.txt" or die;
print GO Dumper(\%types);
print GO Dumper(\%result);
close GO;

my @types = sort (keys %types);


open COM, "> compareCons.csv" or die;

print COM "GOid,GO_name,uberon_c,go-plus_taxon,our_construct,chebi,";
foreach my $type (@types){
  print COM $type.",";
}
print COM "\n";

foreach my $key (keys %goCons){
  next unless ($key =~ /GO:/);
  my $name= $name{$key}; 
  $name =~ s/,/;/g;
  my $c1 = $goCons{$key}{1};
  my $c2 = $goCons{$key}{2};
  my $c3 = $goCons{$key}{3};
  my $c4 = $goCons{$key}{4};
  $c1 =~ s/(NCBITaxon:[0-9]+)/$1\($name{$1}\)/g;
  $c2 =~ s/(NCBITaxon:[0-9]+)/$1\($name{$1}\)/g;
  $c3 =~ s/(NCBITaxon:[0-9]+)/$1\($name{$1}\)/g;
  $c4 =~ s/(NCBITaxon:[0-9]+)/$1\($name{$1}\)/g;
  
  print COM "$key,$name,$c1,$c2,$c3,$c4,";
  foreach my $type (@types){
    my @is  = sort {$a <=> $b} (keys %{$result{$key}{$type}});
    my $string;
    print $type."\t".$key."\t";
    print Dumper(\@is);
    
    foreach my $i (@is){
      my $score = $result{$key}{$type}{$i};
      if ($score > 0){
	$string .= $i."_1".";";
	print $string."\t";
      }
    }
    print "\n";
    print COM $string.",";
  }
  print COM "\n";
}
close COM;

## to summarize the constraints from 4 modes into 1
## the hierachy order is 1, 2&3 and 4

my %combine;
&combineConstraints; 

open COM, "> combineConstraints.txt" or die;
print COM Dumper(\%combine);
close COM;

open COM, "> combineConstraints.csv" or die;
foreach my $goterm ( keys %goCons){
  next unless ($goterm =~ /GO:/);
  my $name = $name{$goterm};
  $name =~ s/,/;/g;

  my $combined = $combine{$goterm};
  print COM "$goterm,$name,$combined\n";
}

close COM;


sub combineConstraints{
  foreach my $goterm (keys %goCons){
    next unless ($goterm =~ /GO:/);
    my $goname = $name{$goterm}; 
    
    my $first = $goCons{$goterm}{2};
    my $second = $goCons{$goterm}{3};
    my $uberon_s = $goCons{$goterm}{1};
    my $chebi_s = $goCons{$goterm}{4};

    my %info1;
    my %info2;
    
    while ($first =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
      my $type1 = $1;
      my $string1 = $2;
      while ($string1 =~ /(NCBITaxon:[0-9]+)/g){
        $info1{$type1}{$1} = 1; 
      }
    }
    
    while ($second =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
      my $type2 = $1;
      my $string2 = $2;
      while ($string2 =~ /(NCBITaxon:[0-9]+)/g){
	$info2{$type2}{$1} = 1; 
      }
    }

    my $combined;

    ### first combine info1 and info2
    ## first info1 only in and info2 gain
    
    my %all;
    foreach my $t1 (keys %{$info1{'only_in_taxon'}}){
      $all{$t1} =1;
    }
    foreach my $t1 (keys %{$info2{'Gain'}}){
      $all{$t1} =1;
    }
    my @arr = keys %all;
    
#    if (exists $result{$term}{'only_in-our_gain'}{1} or $result{$term}{'only_in-our_gain'}{2}){      
      $combined = &mostTaxon(\@arr);      
 #   }
  #  else{ # condition that the same or no overlap
    #  $combined = join (';',@arr);
   # }

    $combined = ">Gain|".$combined;
    
    my $loss;
    ### first combine info1 and info2                                                                  
    ## first info1 only in and info2 gain                                                              
    my %all;
    foreach my $t1 (keys %{$info1{'never_in_taxon'}}){
      $all{$t1} =1;
    }
    foreach my $t1 (keys %{$info2{'Loss'}}){
      $all{$t1} =1;
    }
    my @arr = keys %all;
    
#    if (exists $result{$term}{'never_in-our_loss'}{1} or $result{$term}{'never_in-our_loss'}{2}){      
      $loss = &mostTaxonLoss(\@arr);      
#    }
 #   else{ # condition that the same or no overlap
   #   $loss = join (';',@arr);
  #  }

    my @gain;
    my @loss;
    
    while ($combined =~ /(NCBITaxon:[0-9]+)/g){
      push(@gain,$1);
    }
    while ($loss =~ /(NCBITaxon:[0-9]+)/g){
      push(@loss,$1);
    }
    my $keep;
    foreach my $loss_c (@loss){
      foreach my $gain_c (@gain){
	my $i = &Taxoncompare($loss_c,$gain_c);
	if ($i == 1) { # if gain_c is mom of loss_c
	  $keep .= $loss_c.";";
	  print "for $goterm, we keep LOSS $loss_c\n";
	}
	elsif (($i == -1) or ($i == 2)){
	  print STDERR " while processing $goterm, Strange LOSS: $loss_c GAIN: $gain_c \n";
	}
      }
    }
    if ($keep =~ /NCBI/){
      $combined = $combined.">Loss|".$keep;
    }
    
    ########################## 
    ##OK, finished combining only_in never_in info
    ### next step is to combine uberon and only_in gain info
    
    if ($uberon_s =~ /NCBI/){
      my @arr;
      push(@arr,$uberon_s);
      push(@arr,$combined);

      $combined = &mostTaxonSpecial(\@arr);
      
      print "UBERON $uberon_s, COMBINED $combined\n";
    }
    
    if ($chebi_s =~ /NCBI/){
      unless ($combined =~ /NCBI/){
	$combined = $chebi_s;
	print "$goterm only_CHEBI! $chebi_s\n";
      }
    }
    $combine{$goterm} = $combined;
  }
}



sub mostTaxonLoss{
  my $ref = shift;
  my @taxons = @$ref;
  my $size = @taxons;

  if ($size ==1){
    return $taxons[0];
  }

  else{
    my $rest;
    %skip = {};
    foreach my $i (0..$size-2){
      foreach my $j ($i+1..$size-1){
	if ((exists $skip{$i}) or (exists $skip{$j})){
          next;
	}
        else{
          my $skip_n = &Taxoncompare($taxons[$i],$taxons[$j]);
          if ($skip_n ==1){
            $skip{$i} =1;
            say "skipping $taxons[$i]";
          }
          elsif ($skip_n ==2){
            $skip{$j} =1;
            say "skipping $taxons[$j]";
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


sub constructChebi{
  foreach my $term (keys %name){
    next unless $term =~ /GO:/;
    @rounts = ();
    print "\ntracing up $term mode 4\n"; # use chebi terms for tracing up!

    &gotraceupChebi($term); # use several modes for differnt constraints of trace up;
    print Dumper(\@rounts);

    my $theTaxon = &mostTaxon(\@rounts);
    print "traceup to taxon $theTaxon\n";
    if ($theTaxon =~ /NCBI/){
      $goCons{$goterm}{4} = $theTaxon;
    }    
  }
}


sub gotraceupChebi{
  my $orikey = shift;
  my $mode = 4;
  $orikey =~ /([A-Za-z]+:[0-9]+)$/;
  my $key = $1;

  if (&occurTwice($orikey)){
    push(@rounts,$orikey);
  }

  if (exists $goCons{$key}{$mode}){
    $orikey .= "=>".$goCons{$key}{$mode};
    push(@rounts,$orikey);
    return;
  }
  my $relatedterms = $regulates{$key};
  my @array = split(/;/,$relatedterms);

  foreach my $iterm (@array){
    next unless $iterm =~ /GO:/;
    my $thiskey = $orikey."->".$iterm;

    if (exists $goCons{$iterm}{$mode}){
      $thiskey .= ";C=>".$goCons{$iterm}{$mode};
      push(@rounts,$thiskey);
      #return;
    }
    else{
      &gotraceupChebi($thiskey);
    }
  }
}



## this subroutine for comparing "only_in never_in" constraints in go-plus and our "construct"
## only in; never in
## gain; loss

## compare "only in" and "gain" 1
#"only in" and "loss" 2 
# "never in" and "gain" 3
# "never in" and "loss" 4 

## this subroutine for compare taxon constraints between uberon, chebi with "only_in" and "our constructed"
## uberon with "only_in"
## uberon with "our construct"
## uberon with "chebi"
## chebi with "only_in"
## chebi with "our construct"


sub compareCons{
  foreach my $term (keys %goCons){  
    my $first = $goCons{$term}{2};
    my $second = $goCons{$term}{3};
    my $uberon_s = $goCons{$term}{1};
    my $chebi_s = $goCons{$term}{4};
    
    my %info1;
    my %info2;
    
    while ($first =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
      my $type1 = $1;
      my $string1 = $2;
      while ($string1 =~ /(NCBITaxon:[0-9]+)/g){
	$info1{$type1}{$1} = 1;	
      }
    }
    
    while ($second =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
      my $type2 = $1;
      my $string2 = $2;
      while ($string2 =~ /(NCBITaxon:[0-9]+)/g){
	$info2{$type2}{$1} = 1;	
      }
    }
    my %uberon_c;
    my %chebi_c;
    while ($uberon_s =~ /(NCBITaxon:[0-9]+)/g){
      $uberon_c{$1} =1;
    }
    while ($chebi_s =~ /(NCBITaxon:[0-9]+)/g){
      $chebi_c{$1} =1;
    }
    
    foreach my $t1 (keys %{$info1{'only_in_taxon'}}){
      foreach my $t2 (keys %uberon_c){
	my $i = &Taxoncompare($t1,$t2);
	$result{$term}{'only_in-uberon'}{$i}++;
	$types{'only_in-uberon'}++;
      }
    }
    
    foreach my $t1 (keys %{$info1{'only_in_taxon'}}){
      foreach my $t2 (keys %chebi_c){
	my $i = &Taxoncompare($t1,$t2);
	$result{$term}{'only_in-chebi'}{$i}++;
	$types{'only_in-chebi'}++;
      }
    }
    
    foreach my $t1 (keys %{$info1{'never_in_taxon'}}){
      foreach my $t2 (keys %uberon_c){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'never_in-uberon'}{$i}++;
	$types{'never_in-uberon'}++;
      }
    }
    
    foreach my $t1 (keys %{$info1{'never_in_taxon'}}){
      foreach my $t2 (keys %chebi_c){
	my $i = &Taxoncompare($t1,$t2);
	$result{$term}{'never_in-chebi'}{$i}++;	
	$types{'never_in-chebi'}++;
      }
    }
    
    foreach my $t1 (keys %{$info1{'never_in_taxon'}}){
      foreach my $t2 (keys %{$info2{'Loss'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'never_in-our_loss'}{$i}++;	
	$types{'never_in-our_loss'}++;
      }
    }
    foreach my $t1 (keys %{$info1{'never_in_taxon'}}){
      foreach my $t2 (keys %{$info2{'Gain'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'never_in-our_gain'}{$i}++;
	$types{'never_in-our_gain'}++;
      }
    }
    foreach my $t1 (keys %{$info1{'only_in_taxon'}}){
      foreach my $t2 (keys %{$info2{'Loss'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'only_in-our_loss'}{$i}++;
	$types{'only_in-our_loss'}++;
      }
    }
    foreach my $t1 (keys %{$info1{'only_in_taxon'}}){
      foreach my $t2 (keys %{$info2{'Gain'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'only_in-our_gain'}{$i}++;
	$types{'only_in-our_gain'}++;
      }
    }
    foreach my $t1 (keys %uberon_c){
      foreach my $t2 (keys %{$info2{'Gain'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'uberon-our_gain'}{$i}++;
	$types{'uberon-our_gain'}++;
      }
    }
    foreach my $t1 (keys %uberon_c){
      foreach my $t2 (keys %{$info2{'Loss'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'uberon-our_loss'}{$i}++;
	$types{'uberon-our_loss'}++;
      }
    }
    foreach my $t1 (keys %chebi_c){
      foreach my $t2 (keys %{$info2{'Gain'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'chebi-our_gain'}{$i}++;
	$types{'chebi-our_gain'}++;
      }
    }
    foreach my $t1 (keys %chebi_c){
      foreach my $t2 (keys %{$info2{'Loss'}}){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'chebi-our_loss'}{$i}++;
	$types{'chebi-our_loss'}++;
      }
    }
    foreach my $t1 (keys %uberon_c){
      foreach my $t2 (keys %chebi_c){
	my $i = &Taxoncompare($t1,$t2);	
	$result{$term}{'uberon-chebi'}{$i}++;
	$types{'uberon-chebi'}++;
      }
    }   
  }
}




sub constructGO{
  foreach my $goterm (keys %name){
    next unless $goterm =~ /GO:/;
    @rounts = ();
    
    print "\ntracing up $term mode 1\n";
    &gotraceup($goterm,1); # use several mode for differnt constraints of trace up;
    print Dumper(\@rounts);
    
    my $theTaxon = &mostTaxon(\@rounts);
    print "traceup to taxon $theTaxon\n";
    if ($theTaxon =~ /NCBI/){
      $goCons{$goterm}{1} = $theTaxon;
    }

    
    @rounts = ();
    

    print "\ntracing up $term mode 2\n";
    &gotraceup($goterm,2); # use several mode for differnt constraints of trace up;
    print Dumper(\@rounts);
    
    my $theTaxon = &mostTaxonSpecial(\@rounts);
    print "traceup to taxon $theTaxon\n";
    if ($theTaxon =~ /NCBI/){
      $goCons{$goterm}{2} = $theTaxon;
    }    

    @rounts = ();
    
    print "\ntracing up $term mode 3\n";
    &gotraceupOriginal($goterm); # use huaiyu paul defined rules for tracing up;
    print Dumper(\@rounts);

    my $theTaxon = &mostTaxonSpecial(\@rounts);
    print "traceup to taxon $theTaxon\n";
    if ($theTaxon =~ /NCBI/){
      $goCons{$goterm}{3} = $theTaxon;
    }    
  }
}


sub gotraceupOriginal{
  my $orikey =shift;
  my $mode = 3;
  $orikey =~ /([A-Za-z]+:[0-9]+)$/;
  my $key = $1;
  

  if (&occurTwice($orikey)){
    push(@rounts,$orikey);
    return;
  }
  
  if (exists $goCons{$key}{$mode}){
    $orikey .= "=>".$goCons{$key}{$mode};
    push(@rounts,$orikey);
    return;
  }

  my $relatedterms = $isa{$key};
  my @array = split(/;/,$relatedterms);

  foreach my $iterm (@array){
    next unless $iterm =~ /GO:/;
    my $thiskey = $orikey."->".$iterm;

    if (exists $goCons{$iterm}{$mode}){
      $thiskey .= ";G=>".$goCons{$iterm}{$mode};
      push(@rounts,$thiskey);
     # return;
    }
    else{
      &gotraceupOriginal($thiskey);
    }
  }
}


sub  constructUBERON{
  foreach my $term (keys %name){
    if ($term =~ /PO/) {
      $po_fao_con{$term} = 'NCBITaxon:3193';   
    }
    elsif ($term =~ /FAO/){
      $po_fao_con{$term} = 'NCBITaxon:451864';
    }
    next unless $term =~ /UBERON/;

    if (exists $othercons{$term}){
     #  print "$term\n";
      if ($othercons{$term} =~ /(NCBITaxon:[0-9]+)/){
	$uberonCon{$term} = $1;
      }
    }    
  }
    
  foreach my $term (keys %name){
    @rounts = ();
    
    next unless $term =~ /UBERON/;
    next if (exists $uberonCon{$term});

    print "\ntracing up $term\n";
    &traceup($term);
    print Dumper(\@rounts);
    my $theTaxon = &mostTaxon(\@rounts);
    print "traceup to taxon $theTaxon\n";

    if ($theTaxon){
      $uberonCon{$term} = "->".$theTaxon;
    }
  }
}


sub  constructCL{

  foreach my $term (keys %name){
    @rounts = ();
    next unless $term =~ /CL:/;

    if (exists $othercons{$term}){
      
      while ($othercons{$term} =~ /(UBERON:[0-9]+)/g){
	my $taxonc = $uberonCon{$1};
	if ($taxonc){
	  push(@rounts);
	}
      }
      my $theTaxon = &mostTaxon(\@rounts);
      $po_fao_con{$term} = $theTaxon;
    }
  }
  
  foreach my $term (keys %name){
    @rounts = ();
    
    next unless $term =~ /CL:/;
    next if (exists $po_fao_con{$term});

    print "\ntracing up $term\n";
    &traceupCL($term);
    print Dumper(\@rounts);
    my $theTaxon = &mostTaxon(\@rounts);
    print "traceup to taxon $theTaxon\n";

    if ($theTaxon){
      $po_fao_con{$term} = "->".$theTaxon;
    }
  }
}



sub mostTaxonSpecial{
  my $aref = shift;
  my @arr = @$aref;
  my %taxons;

  foreach my $line (@arr){
    while ($line =~ /([A-Za-z_]+)\|([A-Za-z:;0-9]+)/g){
      my $t = $1;
      my $rest = $2;
      while ($rest =~ /(NCBITaxon:[0-9]+)/g){
	$taxons{$t}{$1} =1;
      }
    }
  }
  my $rest;
  
  foreach my $type (keys %taxons){

    $rest .= ">".$type."|";
    my @taxons = keys %{$taxons{$type}};
    my $size = @taxons;
    if ($size ==1){
      $rest .= $taxons[0];
    }
    
    else{
      
      %skip = {};
      foreach my $i (0..$size-2){
	foreach my $j ($i+1..$size-1){
	  if ((exists $skip{$i}) or (exists $skip{$j})){
	    next;
	  }
	  else{
	    my $skip_n = &Taxoncompare($taxons[$i],$taxons[$j]);

	    if (($type =~ /never/) or ($type =~ /Loss/)){
	      if ($skip_n ==1){
		$skip{$i} =1;
		say "skipping $taxons[$i]";
	      }
	      elsif ($skip_n ==2){
		$skip{$j} =1;
		say "skipping $taxons[$j]";
	      }
	    }	    
	    else{
	      if ($skip_n ==1){
		$skip{$j} =1;
		say "skipping $taxons[$j]";
	      }
	      elsif ($skip_n ==2){
		$skip{$i} =1;
		say "skipping $taxons[$i]";
	      }
	    }
	    
	  }
	}
      }
      foreach my $i (0..$size-1){
	next if (exists $skip{$i});
	my $taxon = $taxons[$i];
	$rest .= $taxon.";";
      }
    }    
  }

  return $rest;

}




sub mostTaxon{
  my $aref = shift;
  my @arr = @$aref;
  my %taxons;
  
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
  else{
    my $rest;
    %skip = {};
    foreach my $i (0..$size-2){
      foreach my $j ($i+1..$size-1){
	if ((exists $skip{$i}) or (exists $skip{$j})){
	  next;
	}
	else{
	  my $skip_n = &Taxoncompare($taxons[$i],$taxons[$j]);
	  if ($skip_n ==1){
	    $skip{$j} =1;
	 #   say "skipping $taxons[$j]";
	  }
	  elsif ($skip_n ==2){
	    $skip{$i} =1;
	  #  say "skipping $taxons[$i]";
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

  my $taxon = $taxon2;
  my $current;

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

  if ($mum1 =~ /$taxon2/){
    return 1;
  }
  elsif ($mum2 =~ /$taxon1/){
    return 2;
  }
  else{
#    print STDERR "$taxon1 mom: $mum1\n";
 #   print STDERR "$taxon2 mom: $mum2\n\n";
    return 3;
  }
}

sub gotraceup{

  my $orikey = shift;
  my $mode = shift;
  
  $orikey =~ /([A-Za-z]+:[0-9]+)$/;
  my $key = $1;


  if (&occurTwice($orikey)){
    push(@rounts,$orikey);
    return;
  }

  given($mode){
    when(1){
      ## mode 1 ## consider only uberon pato fao
      if (exists $goCons{$key}{$mode}){
	$orikey .= "=>".$goCons{$key}{$mode};
	push(@rounts,$orikey);
	return;
      }
      elsif (exists $uberon{$key}){
	while ($uberon{$key} =~ /([A-Za-z]+:[0-9]+)/g){
	  my $goto_term = $1;
	  if (exists $po_fao_con{$goto_term}){
	    $po_fao_con{$goto_term} =~ /([A-Za-z]+:[0-9]+)/;
	    $orikey .= "PF=>".$1;
	    push(@rounts,$orikey);
	    return;
	  }
	  elsif (exists $uberonCon{$goto_term}){
	    $uberonCon{$goto_term} =~ /([A-Za-z]+:[0-9]+)/;
	    $orikey .= "U=>".$1;
	    push(@rounts,$orikey);
	    return;
	  }
	}
      }      
    }
    when(2) {
    ## mode 2 ## consider only_in_taxon and also never_in_taxon      
      if (exists $goCons{$key}{$mode}){
        $orikey .= "=>".$goCons{$key}{$mode};
        push(@rounts,$orikey);
        return;
      }
      elsif(exists $onlyin{$key}){
	foreach my $type (keys %{$onlyin{$key}}){
	  $orikey .= "N=>".$type."|".$onlyin{$key}{$type};
	  push(@rounts,$orikey);
	  return;
	}
      }
    }
    
  }
  
  my $relatedterms = $isa{$key};
  my @array = split(/;/,$relatedterms);
#  print "$key,$relatedterms\n";

  foreach my $iterm (@array){
    next unless $iterm =~ /GO:/;
    my $thiskey = $orikey."->".$iterm;

    given ($mode){
      when(1){
	## mode 1 ## consider only uberon pato fao
	if (exists $goCons{$iterm}{$mode}){
	  $thiskey .= ";G=>".$goCons{$iterm}{$mode};
	  push(@rounts,$thiskey);
	  #return;
	}
	elsif (exists $uberon{$iterm}){
	  while ($uberon{$iterm} =~ /([A-Za-z]+:[0-9]+)/g){
	    my $goto_term = $1;
	    
	    if (exists $pato_fao_con{$goto_term}){
	      $thiskey .= ";PF=>".$pato_fao_con{$goto_term};
	      push(@rounts,$thiskey);
 #         #  return;
	    }
	    elsif (exists $uberonCon{$goto_term}){
	      $uberonCon{$goto_term} =~ /([A-Za-z]+:[0-9]+)/;
	      $thiskey .= ";U=>".$1;
	      push(@rounts,$thiskey);
#	    #  return;
	    }
	  }
	}
	else{
	  &gotraceup($thiskey,$mode);
	}
      }
      when(2){
	#2 ## consider only_in_taxon and also never_in_taxon
	my $key = $iterm;
	
	if (exists $goCons{$key}{$mode}){
	  $thiskey .= ";G=>".$goCons{$key}{$mode};
	  push(@rounts,$thiskey);
#	  # return;
	}
	elsif(exists $onlyin{$key}){
	  foreach my $type (keys %{$onlyin{$key}}){
	    $thiskey .= ";N=>".$type."|".$onlyin{$key}{$type};
	    push(@rounts,$thiskey);
#	    # return;
	  }
	}
	else{
	  &gotraceup($thiskey,$mode);

	}
      }      
    }
    
  }
  return;

}


sub traceup{
  my $orikey =shift;
  $orikey =~ /([A-Za-z]+:[0-9]+)$/;
  my $key = $1;

  if (&occurTwice($orikey)){
    push(@rounts,$orikey);
    return;
  }

  if (exists $uberonCon{$key}){
    $orikey .= "=>".$uberonCon{$key};
    push(@rounts,$orikey);
    return;
  }
  my $relatedterms = $isa{$key};
  my @array = split(/;/,$relatedterms);
#  print "$key,$relatedterms\n";

  foreach my $iterm (@array){
    next unless $iterm =~ /:/;
    my $thiskey = $orikey."->".$iterm;
    if (exists $uberonCon{$iterm}){
      $thiskey = $thiskey."=>".$uberonCon{$iterm};
      push(@rounts,$thiskey);
#      return;
    }
    else{
      &traceup($thiskey);
    }
  }
  return;
}

sub traceupCL{
  my $orikey =shift;
  $orikey =~ /([A-Za-z]+:[0-9]+)$/;
  my $key = $1;

  if (&occurTwice($orikey)){
    push(@rounts,$orikey);
    return;
  }

  if (exists $po_fao_con{$key}){
    $orikey .= "=>".$po_fao_con{$key};
    push(@rounts,$orikey);
    return;
  }
  my $relatedterms = $isa{$key};
  my @array = split(/;/,$relatedterms);
#  print "$key,$relatedterms\n";

  foreach my $iterm (@array){
    next unless $iterm =~ /CL:/;
    my $thiskey = $orikey."->".$iterm;
    if (exists $po_fao_con{$iterm}){
      $thiskey = $thiskey."=>".$po_fao_con{$iterm};
      push(@rounts,$thiskey);
#      return;
    }
    else{
      &traceupCL($thiskey);
    }
  }
  return;
}


sub occurTwice{
  my $key = shift;
  my %store;
  
  while ($key =~ /([A-Z]+:[0-9]+)/g){
    $store{$1}++;

    if ($store{$1} > 1){
      return 1;
    }
  }

  return 0;  
}

## concept for chebi constraints
## construct constraint for chebi term
## chebi terms may have a lot of go terms associated with them
## go terms have evidences
## thus chebi terms get the evidences
## then chebi terms get the constraint
## then we go back to go terms, they get the constriants from chebi


sub chebiConConstruct{

  my %chebis; ## %chebis to store GO terms for each chebi term 
  my %go_ancestor;
  my %go;
  my %an;

  my $sumfile = "sum.taxon.extra.csv";

  foreach my $goterm (keys %chebi){
    next unless ($goterm =~ /GO:/);
    my $chebi = $chebi{$goterm};

    while($chebi =~ /(CHEBI:[0-9]+)/g){
      $chebis{$1}{$goterm} = 1;
      $go{$goterm} = 1;
    }    
  }

  open CHEBIOUT, "> chebifixed.new.txt" or die;

  my %evidences;

  open IN, "< $sumfile" or die;
  while(<IN>){
    chomp;
    my @array = split(/,/);
    my $line = $_;
    my $size = @array;
    my $goterm = $array[0];

    next unless (exists $go{$goterm});
    
    foreach my $i (2..$size-1){
	$evidences{$goterm}{$array[$i]} = 1;
    }
  }
  close IN;


  foreach my $chebiterm (keys %chebis){
    my %evis;
    foreach my $goterm (keys %{$chebis{$chebiterm}}){
	foreach my $evi (keys %{$evidences{$goterm}}){
	    $evis{$evi} =1;
	}
    }
    my @eviTotal = keys %evis;
    my $ancestor = pthr10(@eviTotal);
    $ancestor = lc($ancestor);
    my $anTaxon = $pthr_NCBI{$ancestor};
    $an{$ancestor} = $anTaxon;
    if ($anTaxon =~ /NCBITaxon:[0-9]+/){
	$chebiCon{$chebiterm} = $anTaxon;
    }
  }

  open OUT, "> chebi.an.txt" or die;
  print "output of hash an\n";
  print  Dumper(\%an);
  print  OUT Dumper(\%an);
  close OUT;

  open OUT, "> chebi.Con.txt" or die;
  print "output of hash chebiCon\n";
  print Dumper(\%chebiCon);
  print OUT Dumper(\%chebiCon);
  close OUT;

  foreach my $goterm (%chebi){
    if ($golimit){
      next unless ($goterm =~ /$golimit/);
    }
    next if (exists $go_ancestor{$goterm});
    my $chebiline = $chebi{$goterm};
    my %taxons;

    while($chebiline =~ /(CHEBI:[0-9]+)/g){
      my $taxon = $chebiCon{$1};
      if ($taxon){
	$taxons{$taxon} = 1;
      }
    }
    my @taxons = keys %taxons;
    my $taxon;
    if ($taxons[0]){
      $taxon = &mostTaxonLoss(\@taxons);
    }    
    
    foreach my $ta (@taxons){
      print "LIMIT: $goterm $ta $taxon\n";
      $go_ancestor{$goterm}{$ta} = $taxon;
    }
  }
  
  open OUT, "> chebi.go_ancestor.txt" or die;
  #  print Dumper(\%goCons);
  print "output of hash go_ancestor\n";
  print Dumper(\%go_ancestor);
  print OUT Dumper(\%go_ancestor);
  close OUT;
  exit;
  return;
}
