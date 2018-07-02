use Data::Dumper;
my %type;
my $name;
my $namespace;
my %alt;

my $infile = "./old/go-plus.obo";
open IN, "< $infile" or die;
while(<IN>){
  my $line = $_;
  if ($_ =~ /alt: (GO:[0-9]+)/){
    $alt{$1} = 1;
  }

  if ($_ =~ /\[Term\]/){
    if ($name =~ /GO/){
      next if $namespace =~ /ncbi/;
      $type{$name} = $namespace;
    }
  }
  if ($_ =~ /^id: (GO:[0-9]+)/){
    $name = $1;
  }

  if ($_ =~ /namespace: (.*)/){
    $namespace = $1;
  }

}
close IN;

my %count1; 
my %count2;

my %isa; # isa for taxon id

open IN, "< ./REDO/isa.txt" or die;
while(<IN>){
  my %single = eval($_);
  my $key = (keys %single)[0];
  my $value = $single{$key};
  $value =~ s/;//g;
  $isa{$key} = $value;
}
close IN;


my $yesfile = shift;
my $nofile = shift;
my $confile = shift;
my $countfile = shift;
open OUT, "> $countfile" or die;

my %YES; #  
my %NO;  
my %relationstore;

my %alt;

open IN, "< missingNCBI.dat" or die;
while(<IN>){
  chomp;
  my @array = split(/\t/);
  $alt{$array[0]} = $array[1];
}
close IN;

open IN, "< $yesfile" or die;
while(<IN>){
  chomp;
  my @array = split(/\t|;/);
  my $size =@array;

  foreach my $i (1..$size-1){
    my $ncbi = $array[$i];
    if ($ncbi =~ /NCBI/){
      if ($alt{$ncbi}){
	$ncbi = $alt{$ncbi};
      }
      else {
	unless(exists $isa{$ncbi}){
	  print "check $ncbi\n";
	}
      }

      $YES{$array[0]}{$ncbi} =1;
    }
  }
}
close IN;

open IN, "< $nofile" or die;
while(<IN>){
  chomp;
  my @array = split(/\t|;/);
  my $size =@array;

  foreach my $i (1..$size-1){
    my $ncbi = $array[$i];
    if ($ncbi =~ /NCBI/){
      if ($alt{$ncbi}){
	$ncbi = $alt{$ncbi};
      }
      else {
	unless(exists $isa{$ncbi}){
	  print "check $ncbi\n";
	}
      }

      $NO{$array[0]}{$ncbi} =1;
    }
  }
}
close IN;


open IN, "< $confile" or die;
while(<IN>){
  chomp;
  my $line = $_;
 
  my ($key,$value);
  if ($line =~ /\'(GO:[0-9]+)\' => \'(.*)\'/){
    $key = $1;
    $value = $2;
  }
  my $type = $type{$key};

  my %store;
  while($value =~ />(\w+)\|([NCBITaxon:0-9;]+)/g){
    my $type = $1;
    unless ($type =~ /Loss/){
      $type = "Gain";
    }
    my $value = $2;
    my @varray = split(/;/,$value);
    foreach my $v (@varray){
      if ($v =~ /NCBI/){
	$store{$type}{$v} =1;
      }
    }
  }
  
  my @yesset = keys %{$YES{$key}};
  my @noset = keys %{$NO{$key}};

  ## conflict test 1: is any evidence in yesset a subset or equal to one of the gain ncbitaxon?;
  foreach my $yesid (@yesset){
    my $I =0;
    foreach my $conid (keys %{$store{"Gain"}}){
      if (&isatest($yesid,$conid)){
	$I =1;
	last;
      }
    }
    unless ($I == 1){
      my $toprint = "CONFLICT type 1: $line\tYES: $yesid\n";
      $count1{1}++;
      $count2{1}{$yesid}++;
      $count3{1}{$type}{$key} =1;
      print $toprint;
    }

    if (&isatest($conid,$yesid)){
      my $toprint = "CONFLICT type 1_: $line\tYES: $yesid\n";
      print $toprint;
      $count1{"1_"}++;
      $count2{"1_"}{$yesid}++;
      $count3{"1_"}{$type}{$key} =1;
    }
  }

 ## conflict test 2: is any evidence in yesset a subset or equal to one of the loss ncbitaxon?;
  foreach my $yesid (@yesset){
    next unless $line =~ /Loss/;
    foreach my $conid (keys %{$store{"Loss"}}){
      if (&isatest($yesid,$conid)){

	my $toprint =  "CONFLICT type 2a: $line\tYES: $yesid,child of $conid\n";
	print $toprint;
	$count1{2}++;
	$count2{2}{$yesid}++;
	$count3{2}{$type}{$key} =1;
      }
      elsif (&isatest($conid,$yesid)){

	my $toprint = "CONFLICT type 2b: $line\tYES :$yesid,mom of $conid\n";
	print $toprint;
	$count1{"2_"}++;
	$count2{"2_"}{$yesid}++;
	$count3{"2_"}{$type}{$key} =1;

      }
    }
  }

  ## conflict test 3: is any evidence in noset a subset or equal to one of the gain ncbitaxon?;
  foreach my $noid (@noset){
    my $I =0;
    foreach my $conid (keys %{$store{"Gain"}}){
      if (&isatest($noid,$conid)){
	my $name =~ $name{$noid};
	my $toprint = "CONFLICT type 3a: $line\tNO: $noid,child of $conid\n";

	print $toprint;
	$count1{3}++;
        $count2{3}{$noid}++;
        $count3{3}{$type}{$key} =1;

      }
      elsif (&isatest($conid,$noid)){
	my $name = $name{$noid};
	my $toprint = "CONFLICT type 3b: $line\tNO: $noid,mom of $conid\n";
	print $toprint;
	$count1{"3_"}++;
        $count2{"3_"}{$noid}++;
        $count3{"3_"}{$type}{$key} =1;

      }
    }
  }
  
}
close IN;

# test if the first is a child of the second
sub isatest{
  my $first = shift;
  my $second = shift;
  
  if (exists $relationstore{$first}{$second}){
    return 1;
  }
  if ($first eq $second){
    return 1;
  }
  while(1){
    my $mom = $isa{$first};
    if ($mom eq $second){
      $relationship{$first}{$second} =1;
      return 1;
    }
    elsif ($mom =~ /NCBI/){
      $first = $mom;
    }
    else{
      return 0;
    }
  }
}


print OUT Dumper(\%count1);
#print OUT Dumper(\%count2);

foreach my $key (keys %count2){
  my $count ;
  foreach my $go (sort { $count2{$key}{$b} <=> $count2{$key}{$a} } keys %{$count2{$key}}){
    my $v = $count2{$key}{$go};
    print OUT "$key,$go,$v\n";
    last if $count > 5;
    $count ++;
  }
}

foreach my $key (keys %count3){
  foreach my $type (keys %{$count3{$key}}){
   my @gos = keys %{$count3{$key}{$type}};
   my $n = @gos;
   print OUT "$key,$type,$n\n";
  }
}
print OUT Dumper(\%count3);
close OUT;
