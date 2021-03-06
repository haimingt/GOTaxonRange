# get the sources of GO terms

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

my $sourceFile = "../constraints/goCons_constructlog.txt" ;
my %source;

my $term ;
my $round;

open IN, "< $sourceFile" or die;
while(<IN>){
    if ($_ =~ /\'(GO:[0-9]+)\' => /){
	$term =$1;
    }
    elsif ($_ =~ /'([0-9])' => /){
	$round = $1;
    }
    elsif ($_ =~ /\'(\w+)\' => \'(.*)\'/){
	my $type = $1;
	my $info = $2;
	next unless $round <1;
	while($info =~ /([0-9]):(GO:[0-9]+)/g){
	    $source{$type}{$term}{$round}{$2}{$1} = 1;
	}
    }
}
close IN;

my %goSource;

foreach my $type (sort keys %source){
    foreach my $term (sort keys %{$source{$type}}){
	foreach my $key (keys %{$source{$type}{$term}{0}}){
	    foreach my $n (keys %{$source{$type}{$term}{0}{$key}}){
		if ($n ne 5){
		    $goSource{$type}{$term}{$n} =1;
		}
	    }
	}
	
	my @tmp ; push (@tmp,$term);
	my $current = \@tmp;
	while(1){
	    my ($termref,$sourceref) = &nextLevel($type,$current);
	    foreach my $source (@$sourceref){
		$goSource{$type}{$term}{$source} =1;
	    }

	    if ($termref->[0] =~ /GO:/){
		$current = $termref;
	    }
	    else{
		last;
	    }
	}
    }
}

my %types;

foreach my $type(keys %goSource){
    foreach my $term (keys %{$goSource{$type}}){
	my @source = keys %{$goSource{$type}{$term}};

	my $string = join (';',sort @source);
	$types{$type}{$string}++;
    }
}
print "1: From UBERON,PO,FAO,CL\n";
print "2: From only_in never_in\n";
print "3: From manual curation\n";
print "4: From Chebi\n";

print Dumper(\%types);

sub nextLevel{
    my $type = shift;
    my $aref = shift;

    my %tmpsource; my %tmpTerms;
    foreach my $term (@$aref){
	foreach my $key (keys %{$source{$type}{$term}{0}}){
	    my $indicator = 0;
	    foreach my $n (keys %{$source{$type}{$term}{0}{$key}}){
		if ($n ne 5){
		    $tmpsource{$n} =1;
		    $indicator =1;
		}
	    }
	    if ($indicator ne 1){		
		if ($term ne $key){
		    $tmpTerms{$key} =1;
		}		
	    }
	}
    }
    
    my @tmpterm = keys %tmpTerms;
    my @tmpsource = keys %tmpsource;
    return (\@tmpterm,\@tmpsource);
}
