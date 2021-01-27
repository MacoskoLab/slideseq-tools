#!/usr/bin/perl -w

$numArgs = @ARGV;
if($numArgs != 2)
{   
	print "the number of parameters is not correct!";
	exit(1);
}

$infile	= "$ARGV[0]";
$outfile= "$ARGV[1]";

open(my $fh1, '<', $infile) or die "Could not open file '$infile' $!";
open(my $fh2, '>', $outfile) or die "Could not write file '$outfile' $!";
while (my $row = <$fh1>) {
	chomp $row;
	if (length($row)==0 || index($row, "#") == 0) {
		next;
	}
	
	my @values = split('\t', $row);
	$newline = '';
	foreach my $val (@values) {
		$newline.=$val.',';
	}
	$newline=substr($newline, 0, length($newline)-1);
	print $fh2 $newline."\n";
}
close($fh1);
close($fh2);
