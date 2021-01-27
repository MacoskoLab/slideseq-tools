############################################################################
# Function: convert a .txt file to .csv file                               #
#                                                                          #
# Parameters:                                                              #
# 1) input .txt file                                                       #
# 2) output .csv file                                                      #
#                                                                          #
# Author: Jilong Li                                                        #
# Creation date: 7/25/2019                                                 #
############################################################################

#!/usr/bin/perl -w

#use strict ;
#use warnings ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

#perl(perlfile,arg1,...,argN)

$numArgs = @ARGV;
if($numArgs != 2)
{   
	print "the number of parameters is not correct!";
	exit(1);
}

$infile1	= "$ARGV[0]";
$outfile 	= "$ARGV[1]";

$lt=localtime();
print $lt."\n";

$infile=substr($infile1, 0, length($infile1)-3);
gunzip $infile1 => $infile or die "gunzip failed: $infile1\n";

$lt=localtime();
print $lt."\n";

#$infile
#Puck_181026_13.digital_expression.txt
#Puck_181026_13.digital_expression_summary.txt
$summaryfile=substr($infile, 0, length($infile)-4)."_summary.txt";

open(my $fh3, '<', $summaryfile) or die "Could not open file '$summaryfile' $!";
my %barcodes;
my $i = 0;
while (my $line = <$fh3>) {
	chomp $line;
	if (length($line)==0 || index($line, "#") == 0) {
		next;
	}
	
	$i = $i + 1;
	if ($i == 1) {
		next;
	}
	
	my @values = split('\t', $line);
	if ($values[2] >= 10) {
		$barcodes{$values[0]} = 0;
	}
}
close($fh3);

open(my $fh1, '<', $infile) or die "Could not open file '$infile' $!";
open(my $fh2, '>', $outfile) or die "Could not write file '$outfile' $!";

my @idx = ();
$i = 0;
while (my $row = <$fh1>) {
	chomp $row;
	if (length($row)==0 || index($row, "#") == 0) {
		next;
	}
	
	$i = $i + 1;
	my @values = split('\t', $row);
	if ($i == 1) {
		for (my $j = 0; $j < @values; $j++) {
			if (exists $barcodes{$values[$j]}) {
				push(@idx, $j);
			}
		}
	}
	
	$newline = $values[0];
	foreach my $j (@idx) {
		$newline.=','.$values[$j];
	}
	print $fh2 $newline."\n";
}
close($fh1);
close($fh2);

$lt=localtime();
print $lt."\n";

