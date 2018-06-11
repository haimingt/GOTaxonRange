# this script is for creating table 3,
# the top used NCBITaxon

my $goType = "../statistics/go_Type.txt" ;

open IN, "< $goType" or die;
while(<IN>){
    if ($_ =~ /\'(GO:[0-9]+)\' => \'(\w+)\'/){
	$termType{$1} = $2;
    }
}
close IN;


my %count;

my $cons = "../constraints/goCons.txt" ;
my %goCons;
my $total;

open IN, "< $cons" or die $!;
while(<IN>){
    if ($_ =~ /\'(GO:[0-9]+)\'/){
        $term = $1;
    }
    elsif ($_ =~ /\'(\w+)\' => \'([\w+:;]+)\'/){
	my $type = $1; my $taxon = $2;
	my $ns = $termType{$term};
	while($taxon =~ /(NCBITaxon:[0-9]+)/g){
	    $count{$1}{$type}{$ns}++;
	    $count{$1}{'total'}++;
	    $total++;
	}
    }
}
close IN;

my %name;

my $namefile = "../statistics/termName.txt";
open IN, "< $namefile" or die;
while(<IN>){
    if ($_ =~ /\'(NCBITaxon:[0-9]+)\' => \'(.*)\'/){
	$name{$1} = $2;
    }
}
close IN;

my @sortTerm = sort {$count{$b}{'total'} <=> $count{$a}{'total'}} keys %count;

print "Taxon id and name,,biological_process,molecular_function,cellular_component,Total,Gain/Loss Total\n";

my ($V1,$V2,$V3,$V4,$V5,$V6);

foreach my $i (0..9){
    my $term = $sortTerm[$i];
    my $v1 = $count{$term}{'Gain'}{'biological_process'};
    my $v2 = $count{$term}{'Gain'}{'molecular_function'};
    my $v3 = $count{$term}{'Gain'}{'cellular_component'};
    my $v4 = $count{$term}{'Loss'}{'biological_process'};
    my $v5 = $count{$term}{'Loss'}{'molecular_function'};
    my $v6 = $count{$term}{'Loss'}{'cellular_component'};

    my $t = $v1+$v2+$v3+$v4+$v5+$v6;
    my $per = sprintf("%.2f%%",100*$t/$total); my $name = $name{$term};
    print "$term,Gain,$v1,$v2,$v3,".($v1+$v2+$v3).",".$t."\n";
    print "($name),Loss,$v4,$v5,$v6,".($v4+$v5+$v6).","."$per\n";
    $V1 += $v1;
    $V2 += $v2;
    $V3 += $v3;
    $V4 += $v4;
    $V5 += $v5;
    $V6 += $v6;
    
}

my ($v1,$v2,$v3,$v4,$v5,$v6,$t,$per);
foreach my $i (10..(scalar @sortTerm -1)){
    my $term = $sortTerm[$i];
    $v1 += $count{$term}{'Gain'}{'biological_process'};
    $v2 += $count{$term}{'Gain'}{'molecular_function'};
    $v3 += $count{$term}{'Gain'}{'cellular_component'};
    $v4 += $count{$term}{'Loss'}{'biological_process'};
    $v5 += $count{$term}{'Loss'}{'molecular_function'};
    $v6 += $count{$term}{'Loss'}{'cellular_component'};
}

$t =  $v1+$v2+$v3+$v4+$v5+$v6;
$per =  sprintf("%.2f%%",100*$t/$total);

print "other taxa,Gain,$v1,$v2,$v3,".($v1+$v2+$v3).",".$t."\n";
print ",Loss,$v4,$v5,$v6,".($v4+$v5+$v6).","."$per\n";

$V1 += $v1;
$V2 += $v2;
$V3 += $v3;
$V4 += $v4;
$V5 += $v5;
$V6 += $v6;

print "Total,Gain,$V1,$V2,$V3,".($V1+$V2+$V3).",".$total."\n";
print ",Loss,$V4,$V5,$V6,".($V4+$V5+$V6).","."\n";

