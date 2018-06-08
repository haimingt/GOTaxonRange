## this new version of goplusparse uses the file that lists all children of a go term
## to construct taxon constraints
## the basic idea is that all children of a go term that has a pre-defined taxon constraint should have the same taxon constraint
## then, we can merge these constraints!
## review April 10, 2017

## highest level UBERON term

use List::Util 'shuffle';

# this script for construction of constraint of uberon, plant and fungi terms
my $limit = shift;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

my $infile = "../old/go-plus.obo";

my %name; # store name of GO terms as well as others
my %isa; # store isa information of go terms, point go term to mom
my %relation; # other relationships

my $name;
my @goterm;
my %alt;
my %reverse;
my %newitems;
my @uberonterms;

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
  if ($_ =~ /^id: ((\w+):\w+)/){
    push(@goterm,$1);
    if ($1 =~ /$limit/){
      print "parse goplus, stop at $limit\n";
    }
    if ($1 =~ /UBERON/){
      push(@uberonterms,$1);
    }
  }

  if ($_ =~ /name: (.*)/){
    $name = $1;
  }

  if ($_ =~ /alt_id: (\w+:\w+)/){
    $alt{$goterm[0]} = $1;
    push(@goterm,$1);
  }

  if ($_ =~ /is_a: (\w+:\w+)/){
    foreach my $goterm (@goterm){
      $isa{$goterm} .= $1.";";
      $relation{"is_a"}{$goterm} .= $1.";";
      $reverse{$1} .= $goterm.";";
    }
  }

  if ($line =~ /relationship: (\w+) (\w+:\w+)/){
    my $relation = $1;
    my $c = $2;

    foreach my $goterm (@goterm){
      next unless ($goterm =~ /UBERON|FAO|PO|CL/);
      if ($relation =~ /never_in/){
	$neverin{$goterm} .= $c.";";
      }
      elsif ($relation =~ /taxon/){
	$onlyin{$goterm} .= $c.";";
      }
      else{
	if ($c =~ /UBERON|FAO|PO|CL/){	  
	  $reverse{$c} .= $goterm."*;"; 	  
	}
	$relation{$relation}{$goterm} .= $c.";";
      } 
    }
  }
}
close IN;

my $TRACK;

$isa{'NCBITaxon:5476'} = 'NCBITaxon:4892';
$isa{'NCBITaxon:5833'} = 'NCBITaxon:5820';
$isa{'NCBITaxon:3055'} = 'NCBITaxon:33090';

my %uberonCon;

#print Dumper(%name);
#print Dumper(\%isa);
#print Dumper(\%relation);
#print Dumper(\%onlyin);
#print Dumper(\%neverin);
#print Dumper(\%reverse);

my $count=1;
my @array;
print "$count\tUBERON:0000000;UBERON:0000061\n";

push(@array,"UBERON:0000000");
push(@array,"UBERON:0000061");
#push(@array,"UBERON:0000479");

foreach my $i (1..5){
  my @newarray = @array;
  @array = ();
  foreach my $item (@newarray){
    my $reverse = $reverse{$item};
    while($reverse =~ /(\w+:[0-9]+)/g){
      push(@array,$1);
    }
  }
  $count++;
  print "$count\t";
  print join(";",@array);
  print "\n";
}

#exit;
my @rounts;
&constructUBERON;

sub  constructUBERON{

  &uberonCycle;

  open UBER, "> uberonConstraint1.5.txt" or die;
  print UBER Dumper(\%uberonCon);
  close UBER;

  # manual add

  $uberonCon{'FAO:0000001'} = "NCBITaxon:451864";
  $uberonCon{'PO:0000002'} = 'NCBITaxon:3193';
  $uberonCon{'PO:0000003'} = 'NCBITaxon:3193';

  my %missing;
  foreach my $id (keys %name){
    next unless ($id =~ /UBERON|PO|FAO|CL/);
    my $name = $name{$id};

    unless (exists $uberonCon{$id}){
     $missing{$id} =$name;
    }
  }

  foreach my $key (keys %missing){
    my $name = $name{$key};
    if ($key =~ /FAO:/){
      $uberonCon{$key} = "NCBITaxon:451864;NCBITaxon:5782";
    }
    elsif ($key =~ /PO:/){
      $uberonCon{$key} = "NCBITaxon:3193";
    }
    elsif ($key =~ /UBERON/){
      if ($name =~ /anatomical|cell/){
	next;
      }
      unless (exists $isa{$key}){
	print "no isa? $key, $name\n";
	next;
      }
      else{
	$uberonCon{$key} = "NCBITaxon:6072"; # instead of using 33213, bilateria, use 6072 eumetazoa
      }
    }
    if ($name =~ /Nemato/i){
      $uberonCon{$key} .= "NCBITaxon:6231;";
    }
    if ($name =~ /Protostomia/i){
      $uberonCon{$key} .= "NCBITaxon:33317;";
    }
    if ($name =~ /Diptera/i){
     $uberonCon{$key} .= "NCBITaxon:7147;";
    }
    if ($name =~ /fungal|spore/i){
      $uberonCon{$key} = "NCBITaxon:451864;NCBITaxon:5782";
    }
    if ($name =~ /sperm|oocyte|oogo|egg|blood|visual|animal|gut|lymph|hepato|skeleto|hemato|pancrea|cardio|embryo|neural|neuro/i){
      $uberonCon{$key} = "NCBITaxon:6072";
    }
    if ($name =~ /multicelluar/i){
      $uberonCon{$key} = "NCBITaxon:3193;NCBITaxon:451864;NCBITaxon:5782;NCBITaxon:6072";
    }
  }
  &uberonCycle;
  open UBER, "> uberonConstraint1.6.txt" or die;
  print UBER Dumper(\%uberonCon);
  close UBER;

  my %missing;
  foreach my $id (keys %name){
    next unless ($id =~ /UBERON|PO|FAO|CL/);
    my $name = $name{$id};

    unless (exists $uberonCon{$id}){
      $missing{$id} = $name;
    }
  }
  print Dumper(\%missing);
  
}

sub uberonCycle{

  foreach my $term (keys %onlyin){
    my $current = $term;
    &traceup($current,$onlyin{$term});
  }
  my $count =1;
  while(1){
    print "loop $count\n";
    my @keys = shuffle keys %newitems;
    my $I = @keys;
    last if $I <1 ;
    %newitems = ();
    foreach my $term (@keys){
      my $current = $term;
      &traceup($current,$uberonCon{$term});
    }
    $count++;
  }

  foreach my $term (keys %uberonCon){
    my $current = $term;
    my @array = split(/;/,$uberonCon{$term});
    my $taxon = &mostTaxon(\@array);
    $uberonCon{$current} = $taxon;
  }
}


sub traceup{
  my $current  =shift;
  my $taxon = shift;
  my $target = $reverse{$current};

#  print "CU:$current\tTA:$target\n";

  my @array =split(/;/,$target);
  foreach my $item (@array){
    if ($item =~ /(\w+:\w+)/){
      unless ($uberonCon{$1}){
	$newitems{$1} = 1;
      }
      $uberonCon{$1} .= $taxon.";";
      &traceup($item,$taxon);
    }
  }
}

sub UBERON_go{

  foreach my $thisterm (keys %uberon){
    if ($thisterm =~ $limit){
      print "check point for $goterm\n";
    }
    my $uberonterm = $uberon{$thisterm};
    my %taxons;
    while($uberonterm =~ /([A-Za-z]+:[0-9]+)/g){
      my $taxon = $uberonCon{$1};
      if ($taxon){
	$taxons{$taxon} = 1;
      }
    }
    my @taxons = keys %taxons;
    my $taxon;
    my $size = @taxons;
    if ($size > 1){
      print "UBERON LIMIT: $goterm\n";
      print Dumper(\@taxons);
    }
    if ($taxons[0]){
      $taxon = &mostTaxonLoss(\@taxons);
    }
    if ($taxon){
      $taxon = ">Uberon|".$taxon;
      $goCons{$thisterm}{1} = $taxon;
    }
  }
}

sub ONLYIN_NEVERIN_go{
  foreach my $goterm (keys %onlyin){
    my $con;
    foreach my $type (keys %{$onlyin{$goterm}}){
      if ($type =~ 'only_in'){
	my $value = $onlyin{$goterm}{$type};
	$con .= ">Gain|".$value.";"; 
      }
      elsif ($type =~ 'never_in'){
	my $value = $onlyin{$goterm}{$type};
	$con .= ">Loss|".$value.";"; 
      }
    }
    if ($con){
      $goCons{$goterm}{2} = $con;
    }
  }
}


sub combine_construct_sub{
  foreach my $goterm (keys %allchildren){
    next unless $goterm =~ /GO:/;
    
    my %children;
    my @set1 = keys %{$allchildren{$goterm}};
    my @set2 = keys %{$isachildren{$goterm}};
    foreach my $go (@set1){
      $children{$go}++;
    }
    foreach my $go (@set2){
      $children{$go}++;
    }
    my @children = keys %children;
    foreach my $type (1..5){
      my $con = $goCons{$goterm}{$type};
      if ($con =~ /NCBITaxon:[0-9]+/){
	foreach my $child (@children){
	  if ($limit){
	    next unless( $child =~ /$limit/);
	  }
	  $combine_construct{$child} .= "From $goterm: ".$con.";";
	}
      }
    }
  }

  foreach my $goterm (keys %goCons){

    if ($limit){
      next unless( $goterm =~ /$limit/ );
    }

#    next if (exists $combine_construct{$goterm});
    foreach my $type (1..5){
      my $con = $goCons{$goterm}{$type};
      if ($con =~ /NCBITaxon:[0-9]+/){
	$combine_construct{$goterm} .= "From $goterm: ".$con.";";
      }
    } 
  }

  foreach my $goterm (keys %combine_construct){
    if ($limit){
      next unless( $goterm =~ /$limit/ );
    }
    my $infoline = $combine_construct{$goterm};
    my $newline = &processline($infoline);
    $combineCons{$goterm} = $newline;

    $goCons{$goterm}{5} = $newline;
    print "\n*****\nWorking on $goterm infoline: $infoline\n";
    print "Workong on $goterm newline: $newline\n";
  }
}

sub processline{
  my $infoline = shift;
  my %info;
  while ($infoline =~ /([a-zA-z_]+)\|([NCBITaxon:;0-9]+)/g){
    my $type = $1;
    my $string = $2;
    while ($string =~ /(NCBITaxon:[0-9]+)/g){
      $info{$type}{$1} = 1;
    }
  }
  my $gain;
  my %all;
  foreach my $t1 (keys %{$info{'Gain'}}){
    $all{$t1} =1;
  }
  my @arr = keys %all;
  $gain = &mostTaxon(\@arr);


  my $loss;
  %all = {};
  foreach my $t1 (keys %{$info{'Loss'}}){
    $all{$t1} =1;
  }
  my @arr = keys %all;
  $loss = &mostTaxonLoss(\@arr);
  my @gain;
  my @loss;

  while ($gain =~ /(NCBITaxon:[0-9]+)/g){
    push(@gain,$1);
  }
  while ($loss =~ /(NCBITaxon:[0-9]+)/g){
    push(@loss,$1);
  }
  
  my $uberon;
  my %all;
  foreach my $t1 (keys %{$info{'Uberon'}}){
    $all{$t1} =1;
  }
  my @arr = keys %all;
  $uberon = &mostTaxon(\@arr);
  my @uberon;
  while ($uberon =~ /(NCBITaxon:[0-9]+)/g){
    push(@uberon,$1);
  }

  my %gh;
  my @gain_keep;
  my %visited;
  if (($uberon[0]) and ($gain[0]) ){
    foreach my $gain_c (@gain){
      my $I;
      foreach my $uberon_c (@uberon){
	my $i = &Taxoncompare($gain_c,$uberon_c); 
	if ($i == 1){ # if uberon_c is mom of gain_c
	  $gh{$gain_c} =1 ;
	  $I = 1;
	  $visited{$uberon_c} = 1;
	  last;
	}
	elsif ($i == 2) {# if uberon_c is child of gain_c
	  $gh{$uberon_c} = 1;
	  $I = 2;
	  $visited{$uberon_c} = 1;
	  last;
	}
      }
      if (($I ne 1) and ($I ne 2)){
	$gh{$gain_c} = 1;
      }
    }
    foreach my $uberon_c (@uberon){
      next if (exists $visited{$uberon_c});
      $gh{$uberon_c} = 1;
    }
  }
  @gain_keep = keys %gh;

  unless ($uberon[0]){
    @gain_keep = @gain;
  }
  unless ($gain[0]){
    @gain_keep = @uberon;
  }


  my @loss_keep;
  foreach my $loss_c (@loss){
    foreach my $gain_c (@gain_keep){
      my $i = &Taxoncompare($loss_c,$gain_c);
      if ($i == 1) { # if gain_c is mom of loss_c
	push(@loss_keep,$loss_c);
	print "for $goterm, we keep LOSS $loss_c\n";
      }
      elsif (($i == -1) or ($i == 2)){
	print " while processing $goterm, Strange LOSS: $loss_c GAIN: $gain_c \n";
      }
    }
  }

  my $combined;
  if ($gain_keep[0]){
    $combined = ">Gain|";
    foreach my $gain (@gain_keep){
      $combined .= $gain.";" ;
    }
    if ($loss_keep[0]){
      $combined .= ">Loss|";
      foreach my $loss (@loss_keep){
	$combined .= $loss.";";
      }
    }
  }

  my $chebi;
  my %all;
  foreach my $t1 (keys %{$info{'Chebi'}}){
    $all{$t1} =1;
  }
  my @arr = keys %all;
  $chebi = &mostTaxonLoss(\@arr);
  
  if ($combined){
    return $combined;
  }
  else{
    if ($chebi){
      return ">Chebi|".$chebi;
    }
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
    my $restr;
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
            #say "skipping $taxons[$i]";
          }
          elsif ($skip_n ==2){
            $skip{$j} =1;
            #say "skipping $taxons[$j]";
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
