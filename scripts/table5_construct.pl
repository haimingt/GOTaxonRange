my %isa; # isa for taxon id

open IN, "< ../rawData/GOA/isa.txt" or die;
while(<IN>){
  my %single = eval($_);
  my $key = (keys %single)[0];
  my $value = $single{$key};
  $value =~ s/;//g;
  if ($key =~ /NCBITaxon/){
      $isa{$key} = $value;
  }
}
close IN;


my $yesfile = "../rawData/GOA/go_yes.txt";
my $nofile = "../rawData/GOA/go_no.txt";
my $confile = "../constraints/goCons.txt";

open STDOUT, "> ../statistics/conflicts.goa.txt";

my %YES; #  
my %NO;  
my %relationstore;

my %alt;

open IN, "< ../rawData/GOA/missingNCBI.dat" or die;
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

my $key;
my %store;
my %keepinfo;
my $info;

open IN, "< $confile" or die;
while(<IN>){
  my $line = $_;

  if ($line =~ /\}/){
      $keepinfo{$key} = $info;
      undef $info;
  }
  
  if ($line =~ /\'(GO:[0-9]+)\' => \{/){
      $key = $1;
      $info = $key." constraints:";
  }
  
  if($line =~ /\'([A-Za-z]+)\' => \'([\w;:]+)\'/){
      my $type = $1;
      my $value = $2;
      my @varray = split(/;/,$value);
      foreach my $v (@varray){
	  if ($v =~ /NCBI/){
	      $store{$key}{$type}{$v} =1;
	      $info .= "$type:$v;";
	  }
      }
  }
}
close IN;

foreach my $key (keys %store){
  my @yesset = keys %{$YES{$key}};
  my @noset = keys %{$NO{$key}};
  my $info = $keepinfo{$key};

  foreach my $yesid (@yesset){
    foreach my $conid (keys %{$store{$key}{"Loss"}}){
      if (&isatest($yesid,$conid)){

	my $toprint =  "type_1a: $info\tYES: $yesid,child of $conid\n";
	
	print $toprint;
      }
      elsif (&isatest($conid,$yesid)){

	my $toprint = "type_1b: $info\tYES :$yesid,mom of $conid\n";
	print $toprint;
      }
    }
  }

  ## conflict test 3: is any evidence in noset a subset or equal to one of the gain ncbitaxon?;
  foreach my $noid (@noset){
    my $I =0;
    foreach my $conid (keys %{$store{$key}{"Gain"}}){
      if (&isatest($noid,$conid)){
	my $name =~ $name{$noid};
	my $toprint = "type_2a: $info\tNO: $noid,child of $conid\n";

	print $toprint;
      }
      elsif (&isatest($conid,$noid)){
	my $name = $name{$noid};
	my $toprint = "type_2b: $info\tNO: $noid,mom of $conid\n";

	print $toprint;
      }
    }
  }
  
}


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
