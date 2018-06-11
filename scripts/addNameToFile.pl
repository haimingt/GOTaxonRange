## the script is to add name to files with ids
my $inFile = shift;

my $nameFile = "../statistics/termName.txt" ;

my %termName;
open IN, "< $nameFile" or die;
while(<IN>){
    if ($_ =~ /\'(.*)\' => \'(.*)\'/){
	$termName{$1} = $2;
    }
}
close IN;

open IN, "< $inFile" or die;
while(<IN>){
    s/(\w+:[0-9]+)/$1\($termName{$1}\)/g;
    print;
}
close IN;
