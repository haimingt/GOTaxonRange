## this new version of goplusparse uses the file that lists all children of a go term
## to construct taxon constraints
## the basic idea is that all children of a go term that has a pre-defined taxon constraint should have the same taxon constraint
## then, we can merge these constraints!
## review April 10, 2017

# only use isa relationship as transitive, and other relationships as one time

use List::Util qw(shuffle);
# this script for construction of constraint of uberon, plant and fungi terms
my $limit = shift;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

my $infile = "../old/go-plus.obo";

my %isa; # store isa information of go terms, point go term to mom
my %relation; # other relationships
my %name;

my $goterm;
open IN, "< $infile" or die;
while(<IN>){
  my $line = $_;
  if ($line =~ /\[Term\]/){
    $goterm = "";
  }
  if ($_ =~ /^id: ((\w+):\w+)/){
    $goterm = $1;
    $name{$goterm} =1;
    if ($1 =~ /$limit/){
      print "parse goplus, stop at $limit\n";
    }
  }
  if ($_ =~ /is_a: (\w+:\w+)/){
    $isa{$goterm} .= $1.";";
    $relation{$goterm}{$1} = "is_a";
  }

  if ($line =~ /relationship: (\w+) (\w+:\w+)/){
    my $relation = $1;
    my $c = $2;
    $relation{$c}{$goterm} = $relation;
    $relation{$goterm}{$c} = $relation;
  } 
}
close IN;
my %ancestors;
my $ancestorlist;

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

#Now begins the part of integration of taxon evidences

my $goCons;

my %pthr_NCBI;
open PN, "< pthr_NCBI.txt" or die;
while(<PN>){
  chomp;
  my @array = split(/\t| +/);
  my $size = @array;
  $pthr_NCBI{lc($array[0])} = "NCBITaxon:".$array[$size-1];
}
close PN;

my $I =1;

open SUM, "< Haiming_summary.csv" or die;
while(<SUM>){
  $I =1;
  chomp;
  my @array = split(/,/);
  my $go = $array[1];
  my $info = $array[2];
  $info = lc $info;
  if ($info =~ /gain:([\w;,-]+)/){
    my $catch = $1;
    my @array = split(/;/,$catch);
    foreach my $item (@array){
      my $id = $pthr_NCBI{$item};
      if ($id){
	$goCons{$go}{"G"}{$id} = "Haiming";
	$I =0;
      }
    }
  }
  if ($info =~ /loss:([\w;,-]+)/){
    my $catch = $1;
    my @array = split(/;/,$catch);
    foreach my $item (@array){
      my $id = $pthr_NCBI{$item};
      $goCons{$go}{"L"}{$id} = "Haiming";
      $I =0;
    }
  }
  if ($I ==1){
    print "check $_\n";
  }
}
close SUM;

open OUT, "> goCons.1.txt" or die;
print OUT Dumper(\%goCons);
#exit;
close OUT;

# now deal with onlyin neverin

foreach my $key (keys %relation){
  next unless ($key =~ /NCBI/);
  foreach my $second (keys %relation){
    my $value = $relation{$key}{$second};
    if ($value =~ /never/){
      $goCons{$second}{"L"}{$key} = "neverin";
    }
    elsif ($value =~ /in_taxon/){
      $goCons{$second}{"G"}{$key} = "onlyin";
    }
  }
}
open OUT, "> goCons.2.txt" or die;
print OUT Dumper(\%goCons);
close OUT;

###############################

foreach my $goterm (keys %name){
  if ($limit){
    next unless ($goterm =~ /^GO:00000/);
  }
  $ancestorlist = "";
  &isaup($goterm);
  $ancestors{$goterm} = $ancestorlist;
}

my %children;

open OUT, "> go.ancestors.txt" or die;
print OUT Dumper(\%ancestors);
close OUT;

foreach my $key (keys %ancestors){
  my $value = $ancestors{$key};
  foreach my $i (split(/;/,$value)){
    $children{$i}{$key} =1;
  }
}

open OUT, "> go.children.txt" or die;
foreach my $item (keys %children){
  my @ks = keys %{$children{$item}};
  my $string = join(";",@ks);
  print OUT "$item\t$string\n";
}
close OUT;

exit;
################################

&goCycle;





############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################


sub goCycle{

  foreach my $term (keys %goCons){
    my $current = $term;
    &traceup($current,$chebiCon{$term});

    print "\n\nfinished: traceup $term\n\n";
  }
  my $count =1;

  while(1){
    print "chebiCycle loop $count\n";
    my @keys = shuffle keys %newitems;
    my $I = @keys;
    last if $I <1 ;
    %newitems = ();
    foreach my $term (@keys){
      my $current = $term;
      %visit = (); 
      &traceup($current,$chebiCon{$term});
      print "\n\nfinished: traceup $term\n\n";
    }
    $count++;
  }

  foreach my $term (keys %chebiCon){
    my $current = $term;
    my @array = split(/;/,$chebiCon{$term});
    my $taxon = &mostTaxon(\@array);
    $chebiCon{$current} = $taxon;
  }
}


sub affiliate{
  my @isaa = split(/;/,$ancestorlist);
  foreach my $goisa (@isaa){
    $ancestorlist .= "ADD:$goisa ->".";";
    foreach my $key (keys %{$relation{$goisa}}){
      $ancestorlist .= $key.";";
    }
  }
}

## find all children
sub isaup{
  my $goterm = shift;
  my $I = &check;
  if ($I == 1){
    return;
  }

  if (exists $ancestors{$goterm}){
    $ancestorlist .= $ancestors{$goterm};
    return;
  }
  else{
    my $isas = $isa{$goterm};
    if ($isas){
      $ancestorlist .= $isas;
      my @isaa = split(/;/,$isas);
      foreach my $goisa (@isaa){
	if ($goisa =~ /:/){
	  &isaup($goisa);
	}
      }
    }
  }
}

sub check{
  my %sum;
  my @isaa = split(/;/,$ancestorlist);
  foreach my $goisa (@isaa){
    if ($goisa =~ /:/){
      $sum{$goisa}++;
    }
  }
  foreach my $value (values %sum){
    if ($value > 1){
      return 1;
    }
  }
}


sub constructCHEBI{
  foreach my $key (keys %chebiother){
    @rounts = ();
    my $value = $chebiother{$key};
    my @array = split(/;/,$value);
    foreach my $item (@array){
      if (exists $uberonCon{$item}){
	push(@rounts,$uberonCon{$item});
      }
    }
    my $taxon = &mostTaxon(\@rounts);
    if ($taxon =~ /NCBI/){
      $chebiCon{$key} = $taxon;
    }
  }
  open UBER, "> chebiCon1.txt" or die;
  print UBER Dumper(\%chebiCon);
  close UBER;
  &chebiCycle;

  open UBER, "> chebiCon2.txt" or die;
  print UBER Dumper(\%chebiCon);
  close UBER;

  foreach my $key (keys %chebigo){ 
    my @gos = keys %{$chebigo{$key}};
    @rounts = ();
    foreach my $go (@gos){
      if (exists $YES{$go}){
	push(@rounts,$YES{$go});
      }
    }

    my $taxon;
    if (@rounts){
      $taxon = &taxonCA(\@rounts);
    }
    if ($taxon =~ /NCBI/){
      if (exists $chebiCon{$key}){
	@rounts =();
	push(@rounts,$taxon);
	push(@rounts,$chebiCon{$key});
	$taxon = &taxonCA(\@rounts);
	$chebiCon{$key} = $taxon.";";
      }
      else{
	$chebiCon{$key} .= $taxon.";";
      }
    }
  }
  open UBER, "> chebiCon3.txt" or die;
  print UBER Dumper(\%chebiCon);
  close UBER;
#  &chebiCycle;

#  open UBER, "> chebiCon4.txt" or die;
#  print UBER Dumper(\%chebiCon);
#  close UBER;

}

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
    next unless ($id =~ /UBERON|PO|FAO/);
    my $name = $name{$id};

    unless (exists $uberonCon{$id}){
     $missing{$id} =$name;
    }
  }

  foreach my $key (keys %missing){
    my $name = $name{$key};
    if ($key =~ /FAO:/){
      $uberonCon{$key} = "NCBITaxon:451864";
    }
    elsif ($key =~ /PO:/){
      $uberonCon{$key} = "NCBITaxon:3193";
    }
    elsif ($key =~ /UBERON/){
      if ($name =~ /anatomical|cell/){
	next;
      }
      else{
	$uberonCon{$key} = "NCBITaxon:33213";
      }
    }

  }
  &uberonCycle;
  open UBER, "> uberonConstraint1.6.txt" or die;
  print UBER Dumper(\%uberonCon);
  close UBER;

  my %missing;
  foreach my $id (keys %name){
    next unless ($id =~ /UBERON|PO|FAO/);
    my $name = $name{$id};

    unless (exists $uberonCon{$id}){
      $missing{$id} = $name;
    }
  }
  print Dumper(\%missing);
  
}

sub chebiCycle{

  foreach my $term (keys %chebiCon){
    my $current = $term;
    %visit = ();
    &traceup($current,$chebiCon{$term});

    print "\n\nfinished: traceup $term\n\n";
  }
  my $count =1;

  while(1){
    print "chebiCycle loop $count\n";
    my @keys = shuffle keys %newitems;
    my $I = @keys;
    last if $I <1 ;
    %newitems = ();
    foreach my $term (@keys){
      my $current = $term;
      %visit = (); 
      &traceup($current,$chebiCon{$term});
      print "\n\nfinished: traceup $term\n\n";
    }
    $count++;
  }

  foreach my $term (keys %chebiCon){
    my $current = $term;
    my @array = split(/;/,$chebiCon{$term});
    my $taxon = &mostTaxon(\@array);
    $chebiCon{$current} = $taxon;
  }
}


sub traceup{
  my $current  =shift;
  $visit{$current}++;
  my $taxon = shift;
  my $target = $reverse{$current};
  my $num = $visit{$current};
  print "CU:$current\t$num\tTA:$target\n";

  return if ($current =~ /\*/);
  return if ($num > 1);

  my @array =split(/;/,$target);
  foreach my $item (@array){
    if ($item =~ /(\w+:\w+)/){
      unless ($chebiCon{$1}){
	$newitems{$1} = 1;
      }
      $chebiCon{$1} .= $taxon.";";
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
sub taxonCA{
  my $aref = shift;
  my @arr = @$aref;
  my %taxons;
  
  foreach my $line (@arr){
    while ($line =~ /(NCBITaxon:[0-9]+)/g){
      my $id = $1;
      if (exists $taxonalt{$id}){
	$id = $taxonalt{$id};
      }
      $taxons{$id} =1;
    }
  }
  my @taxons = keys %taxons;
  my $size = @taxons;
  if ($size ==1){
    return $taxons[0];
  }
  else{
    my @anlist;
    my $current = $taxons[0];
    push(@anlist,$current);

    while(1){
      if ($current){
	$taxon = $current;
      }
      if (exists $isa{$taxon}){
	$current = $isa{$taxon};
	$current =~ /(NCBITaxon:[0-9]+)/;
	$current = $1;
	push(@anlist,$current);
      }
      else{
	last;
      }
    }
    my %store;
    foreach my $i (1..$size-1){
      my $current = $taxons[$i];
      
      while(1){
	if ($current){
	  $taxon = $current;
	}
    
	if (exists $isa{$taxon}){
	  $current = $isa{$taxon};
	  $current =~ /(NCBITaxon:[0-9]+)/;
	  $current = $1;
	  $store{$taxons[$i]}{$current} =1;
	}
	else{
	  last;
	}
      }
    }

    foreach my $mom (@anlist){
      my $I = 1;
      foreach my $i (1..$size-1){
	unless (exists $store{$taxons[$i]}{$mom}){
	  $I =0;
	  last;
	}
      }
      if ($I == 1){
	return $mom;
      }
    }
  }
  my $string = join(";",@taxons);
 # &taxonCA(\@taxons);
  return "no mom found for $string";

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
