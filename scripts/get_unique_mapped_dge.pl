#!/usr/bin/perl -w

#use strict ;
#use warnings ;

$numArgs = @ARGV;
if($numArgs != 4)
{   
	print "the number of parameters is not correct!";
	exit(1);
}

$infile			= "$ARGV[0]";  # original dge
$outfile1		= "$ARGV[1]";  # processed dge
$outfile2 		= "$ARGV[2]";  # gene names
$outfile3 		= "$ARGV[3]";  # bead

$lt=localtime();
print $lt."\n";

open(my $fh1, '<', $infile) or die "Could not open file '$infile' $!";
open(my $fh2, '>', $outfile1) or die "Could not write file '$outfile1' $!";
open(my $fh3, '>', $outfile2) or die "Could not write file '$outfile2' $!";
open(my $fh4, '>', $outfile3) or die "Could not write file '$outfile3' $!";

$i = 0;
while (my $row = <$fh1>) {
	chomp $row;
	if (length($row)==0 || index($row, "#") == 0) {
		next;
	}
	
	$i = $i + 1;
	my @values = split('\t', $row);
	if ($i == 1) {
		for (my $j = 1; $j < @values; $j++) {
			print $fh4 $values[$j]."\n";
		}
	}
	else
	{
		print $fh3 $values[0]."\n";
		$newline = $values[1];
		for (my $j = 2; $j < @values; $j++) {
			$newline.="\t".$values[$j];
		}
		print $fh2 $newline."\n";
	}
}
close($fh1);
close($fh2);
close($fh3);
close($fh4);

$lt=localtime();
print $lt."\n";

