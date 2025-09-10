####
#Refine_CytoBand_Liftover.pl
#Written by Nicolas Altemose, 9-9-2025
#
#Uses information from cytoband liftOver output files, CenSat annotations, and chromosome 
#sizes to lift cytoband coordinates from CHM13v2.0 to other karyotypically 'normal' human
#genome assemblies. Optimized for CHM13v2.0 to HG002v1.1 conversion. May not work well on
#edge cases in other assemblies. Uses many heuristic steps to provide rough cytoband
#coordinates.
#
#Outputs a final cytoband bedfile ready for browser upload, as well as a file identifying 
#which bands required special heuristic steps (*ISSUES.txt) and a file providing columns 
#comparing bands to chm13v2.0 for QC purposes (*QC.txt).
#
#Expects new assembly's chromosome nomenclature to be of the form chrN_*
#
#####
use warnings;
use strict;

my $usage = "perl Refine_Cytoband_Liftover_v1.0.pl /path/to/chm13v2.0_cytobands_allchrs.bed /path/to/NewAssembly_CytoBands_LiftedOver.bed /path/to/NewAssembly_CenSat_Track.bed  /path/to/NewAssembly_faidxfile.fai /path/to/outputfile.bed";

#read in file names from command line
my $chm13file = "";
if(defined $ARGV[0]){
	$chm13file = $ARGV[0];
	chomp($chm13file);
}else{
	die "no chm13v2.0 cytobands file specified!\n\nusage: $usage\n\n";
}
my $liftedfile = "";
if(defined $ARGV[1]){
	$liftedfile = $ARGV[1];
	chomp($liftedfile);
}else{
	die "no cytobands file lifted over to the new reference was specified!\n\nusage: $usage\n\n";
}
my $censatfile = "";
if(defined $ARGV[2]){
	$censatfile = $ARGV[2];
	chomp($censatfile);
}else{
	die "no censat file for the new reference was specified!\n\nusage: $usage\n\n";
}
my $chrnamesfile = "";
if(defined $ARGV[3]){
	$chrnamesfile = $ARGV[3];
	chomp($chrnamesfile);
}else{
	die "no reference faidx file was specified! this is needed to get total chromosome sizes \n\nusage: $usage\n\n";
}
my $outfile="";
if(defined $ARGV[4]){
	$outfile = $ARGV[4];
	chomp($outfile);
}else{
	die "no output file was specified!\n\nusage: $usage\n\n";
}


#Read in new assembly chromosome name and size information and store in hashes
my %chrnameconvertNtO; #use hashes to convert between chm13 and new assembly names
my %chrnameconvertOtN;
my %chrlengths;
my @newchrs;
my $prevchr1 = "NONE";

open(IN2, $chrnamesfile) or die "ERROR: Could not open faidx file at $chrnamesfile. Try putting file in this directory.\n";
while(my $line = <IN2>){
	chomp($line);
	$line=~/^(\S+)\t(\S+)\t/ or die "ERROR in faidx file. Failed to parse: $line\n";
	my $chr = $1;
	$chrlengths{$chr}=$2;
	if($chr ne $prevchr1){
		push(@newchrs,$chr);
	}
	if($chr=~/^(\S+)\_/){ #assumes chromosome names are in the format chr1_MATERNAL
		$chrnameconvertOtN{$1}=$chr;
		$chrnameconvertNtO{$chr}=$1;
	}
	$prevchr1=$chr;
}
close IN2;


#Read in CHM13 Cytoband information and store in hashes
my @bands;
my %bandcolor;
my %bandcoordschm13;
my $prevchr0 = 'NONE';
my $prevband = 'NONE';
my %firstbands;
my %lastbands;
my %bandlengths;
my @chrs;
my %passflags;

open(IN1, $chm13file) or die "ERROR: Could not open chm13 cytobands file at $chm13file. Try putting file in this directory.\n";
while(my $line = <IN1>){
	chomp($line);
	$line=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or die "ERROR in cytobands file. Expecting 6 tab separated columns. Failed to parse: $line\n";
	my $chr = $1;
	if(exists $chrnameconvertOtN{$chr}){
		my $band = $chr.' '.$4;
		push(@bands,$band);
		if($chr ne $prevchr0){
			$firstbands{$chr}=$band;
			$lastbands{$prevchr0}=$prevband;
			push(@chrs,$chr);
		}
		$bandcolor{$band}=$5;
		$bandcoordschm13{$band}="$1\t$2\t$3";
		$bandlengths{$band} = $3-$2;
		$prevchr0=$chr;
		$prevband=$band;
		$passflags{$band}=0;
	}
}
$lastbands{$prevchr0}=$prevband;
close IN1;


#Read in new reference's CenSat information and define positions for centromeres, 
#rDNAs, chromosome ends, and large heterochromatic bands on 1, 9, Y
my %censtarts;
my %cenends;
my %rDNAstarts;
my %rDNAends;
my %otherstarts;
my %otherends;
open(IN3, $censatfile) or die "ERROR: Could not open censat file at $censatfile. Try putting file in this directory.\n";
while(my $line = <IN3>){
	chomp($line);
	unless($line=~/^track/){
		$line=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or die "ERROR in censat file. Expecting 9 tab separated columns. Failed to parse: $line\n";
		my $chr = $1;
		my $start = $2;
		my $end = $3;
		my $len = $3-$2;
		my $name = $4;
		my $oldchr = 'NA';
		if( defined $chrnameconvertNtO{$chr}){
			$oldchr=$chrnameconvertNtO{$chr};
		}
		if($name =~ /^active_hor/){
			unless(defined $censtarts{$chr}){
				$censtarts{$chr}=$start;
			}
			$cenends{$chr}=$end;
		}elsif(($name =~ /^rDNA/ || $name =~ /^GAP/) && $len>10000){
			unless(defined $rDNAstarts{$chr}){
				$rDNAstarts{$chr}=$start;
			}
			$rDNAends{$chr}=$end;
		}elsif($oldchr eq 'chr1' && $name eq 'HSat2' && $len >500000){
			$otherends{'chr1 q12'}=$end;
			$otherstarts{'chr1 q21.1'}=$end;
		}elsif($oldchr eq 'chr9' && $name eq 'HSat3' && $len >500000){
			$otherends{'chr9 q12'}=$end;
			$otherstarts{'chr9 q13'}=$end;
		}elsif($oldchr eq 'chrY' && $name eq 'HSat3' && $start >15000000 && $len >100000){
			$otherends{'chrY q11.23'}=$end;
			$otherstarts{'chrY q12'}=$end;
		}
	}
}
close IN3;
#find the midpoints of each centromeric region to serve as the p-q junction
my %cenmids;
foreach my $chr2(@newchrs){
	my $cenS = $censtarts{$chr2};
	my $cenE = $cenends{$chr2};
	my $cenM = $cenS+int(($cenE-$cenS)/2);
	$cenmids{$chr2}=$cenM;
}


#Read in cytoband coordinates lifted over from chm13 to new reference using liftOver
my %liftedbandchrs;
my %liftedbandstarts;
my %liftedbandends;

open(IN4, $liftedfile);
my $prevchr = '';
while(my $line = <IN4>){
	chomp($line);
	$line=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or die "ERROR in lifted over cytobands file. Expecting 6 tab separated columns. Failed to parse: $line\n";
	my $chr = $chrnameconvertNtO{$1};
	my $band = $chr.' '.$4;
	$liftedbandchrs{$band}=$1;
	$liftedbandstarts{$band}=$2;
	$liftedbandends{$band}=$3;
}
close IN4;


#Pass 1: Adjust coordinates for centromeres, rDNAs, beginning, and end of each chr
my %liftedbandchrsOUT = %liftedbandchrs;
my %liftedbandstartsOUT= %liftedbandstarts;
my %liftedbandendsOUT = %liftedbandends;
my $ct = 0;
foreach my $band (@bands){
	$band =~ /^(\S+)\s(\S+)$/;
	my $chr = $1;
	my $bandonly = $2;
	
	if(defined $chrnameconvertOtN{$chr}){
	
		my $newchr = $chrnameconvertOtN{$chr};
		my $col = $bandcolor{$band};
		my $bandlen = $bandlengths{$band};
		#print "$band\t$newchr\t$col\t$bandlen\n";
		if($col eq 'acen'){
			#print "here!\n";
			my $censtart = $censtarts{$newchr};
			my $cenmid = $cenmids{$newchr};
			my $cenend = $cenends{$newchr};
			my $nextband = $bands[$ct+1];
			if($bandonly=~/^p/){
				$liftedbandstartsOUT{$band}= "$censtart";
				$liftedbandendsOUT{$band}= "$cenmid";
				$liftedbandendsOUT{$prevband}= "$censtart";
			}elsif($bandonly=~/^q/){
				$liftedbandstartsOUT{$band}= "$cenmid";
				$liftedbandendsOUT{$band}= "$cenend";
				$liftedbandstartsOUT{$nextband}= "$cenend";
			}
		}elsif($col eq 'stalk'){
			#print "here!\n";
			my $rDNAstart = $rDNAstarts{$newchr};
			my $rDNAend = $rDNAends{$newchr};
			my $nextband = $bands[$ct+1];
			$liftedbandstartsOUT{$band}= "$rDNAstart";
			$liftedbandendsOUT{$band}= "$rDNAend";
			$liftedbandendsOUT{$prevband}= "$rDNAstart";
			$liftedbandstartsOUT{$nextband}= "$rDNAend";
		}elsif($band eq $firstbands{$chr}){
			$liftedbandstartsOUT{$band}= 0;
		}elsif($band eq $lastbands{$chr}){
			$liftedbandendsOUT{$band}= $chrlengths{$newchr};
		}elsif(defined $liftedbandchrs{$band}){
			#print "here\n";
		}elsif(defined $liftedbandstartsOUT{$band}){
			#print "here\n";
		}elsif(defined $liftedbandendsOUT{$band}){
			#print "here\n";
		}
		else{
			$liftedbandstartsOUT{$band}= "NA";
			$liftedbandendsOUT{$band}= "NA";
		}
		#override lifted over coords for censat defined band boundaries (e.g. 1q12)
		if(defined $otherstarts{$band}){
			$liftedbandstartsOUT{$band}=$otherstarts{$band};
		}
		if(defined $otherends{$band}){
			$liftedbandendsOUT{$band}=$otherends{$band};
		}
		$ct++;
	}
}


#Pass 2: Fill in missing boundary coordinates next to positioned bands
#e.g. now that cen bands are placed, the bands next to the cens can get 
#valid start/end positions
my $prevend=0;
$prevchr="NONE";
my $ct2=0;
my $total = scalar(@bands);
my %bandinfo;
foreach my $band (@bands){
	$band =~ /^(\S+)\s(\S+)$/;
	my $chr = $1;
	my $bandonly = $2;
	#print "$band $chr $bandonly\n";
	if(defined $chrnameconvertOtN{$chr}){
	my $newchr = $chrnameconvertOtN{$chr};

	my $col = $bandcolor{$band};
	my $bandlen = $bandlengths{$band};
	my $bandstart = "NA";
	if(defined $liftedbandstartsOUT{$band}){
		$bandstart = $liftedbandstartsOUT{$band};
	}
	my $bandend = "NA";
	if(defined $liftedbandendsOUT{$band}){
		$bandend = $liftedbandendsOUT{$band};
	}
	if($prevchr eq $newchr && $prevend ne 'NA' && $bandstart eq 'NA'){
		$bandstart = $prevend;
	}
	if($band ne $lastbands{$chr}){
		my $nextband = $bands[$ct2+1];
		if(defined $liftedbandstartsOUT{$nextband} && $bandend eq 'NA'){
			$bandend = $liftedbandstartsOUT{$nextband};
		}
	}
	my $newbandlen = 'NA';
	my $bandlendiff = 'NA';
	my $overlap = 'NA';
	if($bandstart ne 'NA' && $bandend ne 'NA'){
		$newbandlen = $bandend-$bandstart;
		$bandlendiff = $newbandlen - $bandlen;
		if($bandstart ne "0" && "prevend" ne "NA"){
			$overlap = $bandstart - $prevend;
		}
	}
	$liftedbandstartsOUT{$band}=$bandstart;
	$liftedbandendsOUT{$band}=$bandend;
	$bandinfo{$band}= "$newchr\t$bandstart\t$bandend\t$bandonly\t$col\t$bandlen\n";
	$prevend=$bandend;
	$prevchr=$newchr;
	$ct2++;
	}
	
}
#Fix the 9p13.1-9p12-9p11.2 region to have fixed band lengths
if(defined $liftedbandstartsOUT{"chr9 p13.1"} && $liftedbandstartsOUT{"chr9 p13.1"} ne 'NA'){
	$liftedbandendsOUT{"chr9 p13.1"}=$liftedbandstartsOUT{"chr9 p13.1"}+1100000;
	$liftedbandstartsOUT{"chr9 p12"}=$liftedbandstartsOUT{"chr9 p13.1"}+1100000;
	$liftedbandendsOUT{"chr9 p12"}=$liftedbandstartsOUT{"chr9 p13.1"}+2100000;
	$liftedbandstartsOUT{"chr9 p11.2"}=$liftedbandstartsOUT{"chr9 p13.1"}+2100000;
	$passflags{"chr9 p13.1"}=5;
	$passflags{"chr9 p12"}=5;
	$passflags{"chr9 p11.2"}=5;
}else{
	print "ERROR resolving 9p13.1-9p12-9p11.2\n";
}


#Pass 3: For adjacent bands separated by a boundary that failed to lift over,
#provide an approximate boundary location that preserves the relative size
#of each band in CHM13

$prevend=0;
my %printinfo;

for(my $i = 0; $i<scalar(@bands);$i++){
	my $band = $bands[$i];
	$band =~ /^(\S+)\s(\S+)$/;
	my $chr = $1;
	my $bandonly = $2;
	my $newchr = $chrnameconvertOtN{$chr};
	my $col = $bandcolor{$band};
	my $bandlen = $bandlengths{$band};
	my $bandstart = $liftedbandstartsOUT{$band};
	my $bandend = $liftedbandendsOUT{$band};
	my $passflag = $passflags{$band};

	if($i>0 && $i<scalar(@bands)-1){
		my $prevband = $bands[$i-1];
		my $nextband = $bands[$i+1];
		$prevband =~ /^(\S+)\s(\S+)$/;
		my $prevchr0 = $1;
		$nextband =~ /^(\S+)\s(\S+)$/;
		my $nextchr0 = $1;

		my $nextbandstart = $liftedbandstartsOUT{$nextband};
		my $prevbandstart = $liftedbandstartsOUT{$prevband};
		my $nextbandend = $liftedbandendsOUT{$nextband};
		my $prevbandend = $liftedbandendsOUT{$prevband};
		my $nextbandlen = $bandlengths{$nextband};
		my $prevbandlen = $bandlengths{$prevband};
	
		if($chr eq $nextchr0 && $bandend eq 'NA' && $nextbandstart eq 'NA' && $bandstart ne 'NA' && $nextbandend ne 'NA'){
			my $diff = $nextbandend - $bandstart;
			my $newlength = int($bandlen*$diff/($bandlen+$nextbandlen));
			$liftedbandendsOUT{$band}=$bandstart+$newlength;
			$liftedbandstartsOUT{$nextband}=$bandstart+$newlength;
			$passflags{$band}=3;
			$passflags{$nextband}=3;
		}elsif($i<scalar(@bands)-2 && $bandstart eq 'NA' && $bandend eq 'NA' && $chr eq $nextchr0){
			my $j=0;
			my $anchorflag = 0;
			my $regionend=0;
			my $preliftoverTotalLength=$prevbandlen+$bandlen;
			until($anchorflag==1){
				$j++;
				my $band0 = $bands[$i+$j];
				$preliftoverTotalLength+= $bandlengths{$band0};
				if($liftedbandendsOUT{$band0} ne 'NA'){
					$regionend = $liftedbandendsOUT{$band0};
					$anchorflag = 1;
				}
			}
			my $regionlength = $regionend-$prevbandstart;
			$liftedbandendsOUT{$prevband} = $prevbandstart + int($prevbandlen*$regionlength/$preliftoverTotalLength);
			$liftedbandstartsOUT{$band} = $liftedbandendsOUT{$prevband};
			$liftedbandendsOUT{$band} = $liftedbandendsOUT{$prevband} + int($bandlen*$regionlength/$preliftoverTotalLength);
			$liftedbandstartsOUT{$nextband} = $liftedbandendsOUT{$prevband} + int($bandlen*$regionlength/$preliftoverTotalLength);
			if($j>=2){
				for(my $k=2;$k<=$j;$k++){
					my $band0 = $bands[$i+$k];
					my $bandlen0 = $bandlengths{$band0};
					my $prevband0 = $bands[$i+$k-1];
					my $prevband0len = $bandlengths{$prevband0};
					$liftedbandendsOUT{$prevband0} = $liftedbandstartsOUT{$prevband0} + int($prevband0len*$regionlength/$preliftoverTotalLength);
					$liftedbandstartsOUT{$band0} = $liftedbandendsOUT{$prevband0};
					if($k<$j){
						$liftedbandendsOUT{$band0} = $liftedbandstartsOUT{$band0}+ int($bandlen0*$regionlength/$preliftoverTotalLength);
					}
					$passflags{$band0}=4;
				}
			}
			$passflags{$prevband}=4;
			$passflags{$band}=4;
			$passflags{$nextband}=4;
		}		
	}
}


#Pass 4: Adjust adjacent bands to be flush against each other
my %banddiffs;
for(my $i = 0; $i<scalar(@bands);$i++){
	my $band = $bands[$i];
	$band =~ /^(\S+)\s(\S+)$/;
	my $chr = $1;
	my $bandonly = $2;
	my $newchr = $chrnameconvertOtN{$chr};
	my $col = $bandcolor{$band};
	my $bandlen = $bandlengths{$band};
	my $bandstart = $liftedbandstartsOUT{$band};
	my $bandend = $liftedbandendsOUT{$band};
	my $passflag = $passflags{$band};
	$chr=~/^chr(\S+)$/;
	my $chrnum = $1;
	my $newbandname = $chrnum.$bandonly;
	my $newbandlen = $bandend-$bandstart;
	my $gap=0;
	my $prevband = $bands[$i-1];
	my $prevbandend = $liftedbandendsOUT{$prevband};
	my $diff = $newbandlen-$bandlen;
	$banddiffs{$band}=$diff;
	
	unless($band eq $firstbands{$chr}){
		$gap = $bandstart-$prevbandend;
	
		if($gap<10){
			$liftedbandstartsOUT{$band}=$liftedbandendsOUT{$prevband};
		}else{
			my $prevdiff = $banddiffs{$prevband};
			if($diff < $prevdiff){
				$liftedbandstartsOUT{$band}=$liftedbandendsOUT{$prevband};
			}else{
				$liftedbandendsOUT{$prevband}=$liftedbandstartsOUT{$band};
			}
			$passflags{$band}=6;
		}
	}
}


#Now print the final coordinates for each band
my $issuesfile = $outfile."ISSUES.txt";
my $qcfile = $outfile."QC.txt";
open(OUT,'>'.$outfile);
open(ISSUES,'>'.$issuesfile);
open(QC,'>'.$qcfile);
for(my $i = 0; $i<scalar(@bands);$i++){
	my $band = $bands[$i];
	$band =~ /^(\S+)\s(\S+)$/;
	my $chr = $1;
	my $bandonly = $2;
	my $newchr = $chrnameconvertOtN{$chr};
	my $col = $bandcolor{$band};
	my $bandlen = $bandlengths{$band};
	my $bandstart = $liftedbandstartsOUT{$band};
	my $bandend = $liftedbandendsOUT{$band};
	my $passflag = $passflags{$band};
	$chr=~/^chr(\S+)$/;
	my $chrnum = $1;
	my $newbandname = $chrnum.$bandonly;
	my $newbandlen = $bandend-$bandstart;
	my $gap=0;
	unless($band eq $firstbands{$chr}){
		my $prevband = $bands[$i-1];
		my $prevbandend = $liftedbandendsOUT{$prevband};
		$gap = $bandstart-$prevbandend;
	}
	my $ratio = sprintf("%3d", 100 * abs($newbandlen/$bandlen-1));
	my $diff = abs($newbandlen-$bandlen);
	print OUT "$newchr\t$bandstart\t$bandend\t$bandonly\t$col\t$newbandname\n";
	print QC "$newchr\t$bandstart\t$bandend\t$bandonly\t$col\t$bandlen\t$newbandlen\t$diff\t$ratio\t$gap\n";
	if($passflag == 3){
		print ISSUES "$newchr\t$bandstart\t$bandend\t$bandonly\tSingleBoundarySetToPreserveRelativeLengths\n";
	}
	if($passflag == 4){
		print ISSUES "$newchr\t$bandstart\t$bandend\t$bandonly\tMultipleBoundariesSetToPreserveRelativeLengths\n";
	}
	if($passflag == 5){
		print ISSUES "$newchr\t$bandstart\t$bandend\t$bandonly\tMultipleBoundariesSetToPreserveFixedLengths\n";
	}
	if($passflag == 6){
		print ISSUES "$newchr\t$bandstart\t$bandend\t$bandonly\tBoundaryMovedToBeFlushWithNeighbor\n";
	}
	if($bandstart eq 'NA' || $bandend eq 'NA'){
		print ISSUES "$newchr\t$bandstart\t$bandend\t$bandonly\tBoundaryNotConvertedSuccessfully_ManualAdjustmentNeeded\n";
		print "$newbandname\tBoundaryNotConvertedSuccessfully_ManualAdjustmentNeeded\n";
	}
}
close OUT;
close ISSUES;
close QC;

print "Printed to $outfile\n";
print "Issues raised in $issuesfile\n";
print "QC metrics printed to $qcfile, with columns CHM13_BandLength Lifted_BandLength AbsDiffInLength AbsPctDiffInLength GapFromPreviousBand\n";



exit;