my $limit = shift;

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
            $source{$term}{$2} = $type.":".$1.":".$2;
        }
    }
}
close IN;


my $n=0;
my $current = $limit;
while(1){
    my $next = &process($current);

    if ($next){
	$current = $next;
    }
    else{
	last;
    }
    $n++;
    last if $n > 20;
}

sub process{
    my $input = shift;
    my %tmp;
    while($input =~ /(GO:[0-9]+)/g){
	$tmp{$1} = 1;
    }
    my %target;

    foreach my $key (keys %tmp){
	foreach my $next (keys %{$source{$key}}){
	    my $val = $source{$key}{$next};
	    my @t = split(/:/,$val);


	    print "$key -> $next -> $val\n";		
	    if ($t[1] ne 5){
	    }
	    else{
		$target{$next} = 1;
	    }
	}
    }
    my $returnstring;
    
    foreach my $key (keys %target){
	$returnstring .= $key.";";
    }
    return $returnstring;
}
