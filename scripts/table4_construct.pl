# get the sources of GO terms
# in the seed list

use Data::Dumper;

my $goType = "../statistics/go_Type.txt" ;

my %termType;
open IN, "< $goType" or die;
while(<IN>){
    if ($_ =~ /\'(GO:[0-9]+)\' => \'(\w+)\'/){
        $termType{$1} = $2;
    }
}
close IN;

my %goterms;
my %statistics;

my $manualfile = "../rawData/manualCurationList_longform.csv";

open IN, "< $manualfile" or die;
while(<IN>){
    chomp;
    my @a = split(/,/);
    $statistics{$termType{$a[0]}}{'ManualCuration'}{$a[0]} = 1;
    $goterms{$termType{$a[0]}}{$a[0]} =1;
}
close IN;

my %chebi;
my $chebifile = "../constraints/chebi.Con.txt";

open IN, "< $chebifile" or die;
while(<IN>){
    chomp;
    if ($_ =~ /(CHEBI:[0-9]+)/){
	$chebi{$1} =1;
    }
}
close IN;

my %uberon;
my $uberonfile = "../constraints/uberon_PO_FAO_CL.constraints.txt";

open IN, "< $uberonfile" or die;
while(<IN>){
    chomp;
    if ($_ =~ /(\w+:[0-9]+)/){
	$uberon{$1} =1;
    }
}
close IN;

my $goplusfile = "../rawData/go-plus.obo";

my $name;
my $term;

open IN, "< $goplusfile" or die "cannot open $goplusfile\n";
while(<IN>){
    my $line = $_;
    chomp($line);

    if ($_ =~ /\[Term\]/){
	$term = '';
    }
    elsif ($_ =~ /^id: (GO:[0-9]+)/){
	$term = $1;
    }
    elsif ($line =~ /relationship: (\w+) ([A-Za-z]+:[0-9]+)/){
      my $relation = $1;
      my $m = $2;

#      if ($m =~ /UBERON:/ or $m =~ /PO:/ or $m =~ /FAO/){
#	  if ($termType{$term} =~ /molecular/){
#	      print OUT "$term\t$line\n";
#	  }
#      }
      
      next unless $term =~ /GO:/;
      if (exists $chebi{$m}){
	  $statistics{$termType{$term}}{'Chebi/PR'}{$term} = 1;
	  $goterms{$termType{$term}}{$term} =1;
      }
      if (exists $uberon{$m}){
	  $statistics{$termType{$term}}{'UBERON/PO/FAO/CL'}{$term} = 1;
	  $goterms{$termType{$term}}{$term} =1;
      }
      if ($m =~ /NCBITaxon/){
	  $statistics{$termType{$term}}{'Only/Never_in taxon'}{$term} = 1;
	  $goterms{$termType{$term}}{$term} =1;
      }      
    }
}
close IN;

my @type = ("molecular_function","cellular_component","biological_process");
my @source = ("ManualCuration","Only/Never_in taxon","UBERON/PO/FAO/CL","Chebi/PR");

print ",";
foreach my $source (@source){
    print "$source,";
}
print "Total from all sources\n";

my $total;

foreach my $type (@type){
    my @keys = keys %{$goterms{$type}};
    my $n = @keys;
    $total += $n;
}

my %col;

foreach my $type (@type){
    my $t;
    print "$type,";
    foreach my $source (@source){
	my @keys = keys %{$statistics{$type}{$source}};
	my $n = @keys;
	$t += $n;
	print "$n,";
	$col{$source} += $n;
    }
    my $per =  sprintf("%.2f%%",100*$t/$total);

    print "$t \($per\)\n";
}

print "Total,";

foreach my $source(@source){

    my $n = $col{$source};
    my $per =  sprintf("%.2f%%",100*$n/$total);

    print "$n \($per\),";
}

print "$total\n";

