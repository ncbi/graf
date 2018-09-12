#!/usr/local/bin/perl

my $disclaim = << "EOF";
    ==================================================================================== 
                               PUBLIC DOMAIN NOTICE
                National Center for Biotechnology Information

        This software/database is a "United States Government Work" under the
        terms of the United States Copyright Act.  It was written as part of
        the author's official duties as a United States Government employee and
        thus cannot be copyrighted.  This software/database is freely available
        to the public for use. The National Library of Medicine and the U.S.
        Government have not placed any restriction on its use or reproduction.
        Although all reasonable efforts have been taken to ensure the accuracy
        and reliability of the software and data, the NLM and the U.S.
        Government do not and cannot warrant the performance or results that
        may be obtained by using this software or data. The NLM and the U.S.
        Government disclaim all warranties, express or implied, including
        warranties of performance, merchantability or fitness for any particular
        purpose.

        Please cite the author in any work or product based on this material.

        Author: Yumi (Jimmy) Jin (jinyu\@ncbi.nlm.nih.gov)
        File Description: script to plot graphs to show subject populations.
        Date: 04/12/2018
    ====================================================================================
EOF

use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use strict;
use Time::HiRes qw(gettimeofday);
use Data::Dumper;
use GD;
use GD::Text;
use GD::Graph;
use GD::Graph::colour;
use GD::Graph::lines;

#------------------------------------------- Take and set parameters ----------------------------------------#
my $pi = 3.1415926535;
my $epsilon = 0.0000000001;
my $maxRounds = 100;
my @startTime = gettimeofday;
my @finishTime = gettimeofday;

my ($inGrafFile, $outputFile, $params) = GetParameters(@ARGV);

unless ($inGrafFile && $outputFile) {
    print "\n$disclaim\n";
    print "\nUsage: PlotPopulations.pl <input file> <output file> [Options]\n\n";
    print "Note:\n";
    print "      Output file should be either a .png file or a .txt file.\n";
    print "      If the output file is a .png file, the script will plot the results to a graph and save the graph to the file.\n";
    print "      If the output file is a .txt file, the script will save the calculated subject ancestry components to the file.\n\n";
    print "Options:";
    print "\n    Set window size in pixels\n";
    print "        -gw      graph width\n";
    print "\n    Set graph axis limits\n";
    print "        -xmin    min x value\n";
    print "        -xmax    max x value\n";
    print "        -ymin    min y value\n";
    print "        -ymax    max y value\n";
    print "\n    Set a rectangle area to retrieve subjects for graph of DG1 vs. DG2\n";
    print "        -xcmin   min x value\n";
    print "        -xcmax   max x value\n";
    print "        -ycmin   min y value\n";
    print "        -ycmax   max y value\n";
    print "        -isByd   0 or 1 \n";
    print "                 0:  retrieve subjects whose values are within the above rectangle (default value)\n";
    print "                 1:  retrieve subjects whose values are beyond the above rectangle\n";
    print "\n    Set minimum and maximum numbers of genotyped fingerprint SNPs for samples to be processed\n";
    print "        -minsnp  minimum number of SNPs with genotypes\n";
    print "        -maxsnp  maximum number of SNPs with genotypes\n";
    print "\n    Set population cutoff lines\n";
    print "        -ecut    proportion: cutoff European proportion dividing Europeans from other population. Default 87%.\n";
    print "        -fcut    proportion: cutoff African proportion dividing Africans from other population. Default 95%.\n";
    print "                             Set it to -1 to combine African and African American populations\n";
    print "        -acut    proportion: cutoff East Asian proportion dividing East Asians from other populations. Default 95%.\n";
    print "                             Set it to -1 to combine East Asian and Other Asian populations\n";
    print "        -ohcut   proportion: cutoff African proportion dividing Hispanics from Other population. Default 13%.\n";
    print "        -fhcut   proportion: cutoff African proportion dividing Hispanics from African Americans. Default 40%.\n";
    print "\n    Select some self-reported populations (by IDs) to be highlighted on the graph\n";
    print "        -pops    comma separated population IDs, e.g., -pops 1,3,4 -> highlight populations #1, #3 and #4 \n";
    print "\n    Select self-reported populations (by IDs) to show areas including 95% dbGaP subjects with genotypes of at least 4000 fingerprint SNPs\n";
    print "        -areas   comma separated dbGaP self-population IDs, e.g., -areas 1,3\n";
    print "                     -> show areas that include 95% dbGaP subjects with self-reported populations #1 and #3\n";
    print "                        1: European/White/Caucasian       2: African (Ghana/Yoruba)\n";
    print "                        3: East Asian (Chinese/Japanese)  4: African American/Black\n";
    print "                        5: Puerto Rican/Dominican         6: Mexican/Latino\n";
    print "                        7: Asian/Pacific Islander         8: Asian Indian/Pakistani\n";
    print "\n    Select which score to show on the y-axis\n";
    print "        -gd4     1 or 0.  1: show GD4 on y-axis;  0: show GD2\n";
    print "\n    Set population cutoff lines\n";
    print "        -cutoff  1 or 0.  1: show cutoff lines;  0: hide cutoff lines\n";
    print "\n    Rotate the plot on x-axis by a certain angle\n";
    print "        -rotx    angle in degrees \n";
    print "\n    Set the size (diameter) of each dot that represents each subject\n";
    print "        -dot     pixels\n";
    print "\n    The input file with self-reported subject race information\n";
    print "        -spf     a file with two columns: subject and self-reported population\n";
    print "\n";
    exit 0;
}

my $showSbjs = 0;
if ($outputFile =~ /\.png$/) {
    $showSbjs = 0;
}
elsif ($outputFile =~ /\.txt$/) {
    $showSbjs = 1;
}
else {
    print "\nERROR: Invalid output file $outputFile\n";
    print "       It should be either a .png file (for saving results in a graph) or a .txt file (for saving results in a table).\n\n";
    exit;
}

my $gWidth   = $params->{"gw"};
my $xMin     = $params->{"xmin"};
my $xMax     = $params->{"xmax"};
my $yMin     = $params->{"ymin"};
my $yMax     = $params->{"ymax"};
my $xCutMin  = $params->{"xcmin"};
my $xCutMax  = $params->{"xcmax"};
my $yCutMin  = $params->{"ycmin"};
my $yCutMax  = $params->{"ycmax"};
my $minSnps  = $params->{"minsnp"};
my $maxSnps  = $params->{"maxsnp"};

my $showGd4  = $params->{"gd4"};
my $rotx     = $params->{"rotx"};
my $dotSize  = $params->{"dot"};
my $selRaces = $params->{"pops"};
my $isBeyond = $params->{"isByd"};
my $raceFile = $params->{"spf"};
my $showCutoff = $params->{"cutoff"};
my $showAreas = $params->{"areas"};

$xMin = 1.2 unless ($xMin);
$xMax = 2.0 unless ($xMax);

if ($showSbjs) {
    $showGd4 = 0;
}

if ($showGd4) {
    $yMin = -0.3 unless ($yMin);
    $yMax =  0.5 unless ($yMax);
}
else {
    $yMin = 1.0 unless ($yMin);
    $yMax = 1.8 unless ($yMax);
}
$yMin = 0.8 if (!$showGd4 && $rotx && $rotx > 120 && $rotx < 240);

# Round to the tenth
$xMin = int($xMin * 10 + 0.5) / 10;
$xMax = int($xMax * 10 + 0.5) / 10;
$yMax = int($yMax * 10 + 0.5) / 10;
if ($yMin < 0) {
    $yMin = int($yMin * 10 - 0.5) / 10;
}
else {
    $yMin = int($yMin * 10 + 0.5) / 10;
}

if ($xMin && $xMax && $xMin >= $xMax) {
    print "\nERROR: x-min ($xMin) is not smaller than x-max ($xMax).\n\n";
    exit;
}
if ($yMin && $yMax && $yMin >= $yMax) {
    print "\nERROR: y-min ($yMin) is not smaller than y-max ($yMax).\n\n";
    exit;
}


my $xRange = $xMax - $xMin;
my $yRange = $yMax - $yMin;
my $gHeight  = $gWidth * $yRange / $xRange;

$xCutMin = 0 unless ($xCutMin);
$xCutMax = 0 unless ($xCutMax);
$yCutMin = 0 unless ($yCutMin);
$yCutMax = 0 unless ($yCutMax);

$minSnps = 0 unless ($minSnps);
$maxSnps = 10000 unless ($maxSnps);

if ($xCutMin && $xCutMax && $xCutMin > $xCutMax) {
    $xCutMin = 0;
    $xCutMax = 0;
}

if ($yCutMin && $yCutMax && $yCutMin > $yCutMax) {
    $yCutMin = 0;
    $yCutMax = 0;
}

if ($dotSize && $dotSize =~ /^\d+$/) {
    $dotSize = 1 if ($dotSize < 1);
    $dotSize = 10 if ($dotSize > 10);
}
else {
    $dotSize = 2;
}

# When the user doesn't select races to show, nor which regions to get subjects, retrive all subjects
my $isSelAll = 1;
$isSelAll = 0 if ($selRaces || $xCutMin || $xCutMax || $yCutMin || $yCutMax);

my ($sec, $min, $hour, $day, $month, $year, $wday, $yday, $isdst) = localtime;
$year += 1900;
$month += 1;

if ($showSbjs) {
    if (-e $outputFile) {
	print "\nOutput file $outputFile already exists. Overwrite it? (y/n)\n";
	my $ans = <STDIN>;
	if ($ans !~ /^Y/i) {
	    my $base = $outputFile;
	    $base = $1 if ($outputFile =~ /(.+)\.txt$/);
	    $outputFile = sprintf("%s_%04d%02d%02d_%02d%02d.txt", $base, $year, $month, $day, $hour, $min);
	}
    }
}
else {
    if (-e $outputFile) {
	print "\nOutput file $outputFile already exists. Overwrite it? (y/n)\n";
	my $ans = <STDIN>;
	if ($ans !~ /^Y/i) {
	    my $base = $outputFile;
	    $base = $1 if ($outputFile =~ /(.+)\.png$/);
	    $outputFile = sprintf("%s_%04d%02d%02d_%02d%02d.png", $base, $year, $month, $day, $hour, $min);
	}
    }
}

my $xMinorTicks = $params->{"minorTicks"};
$xMinorTicks = 5 unless ($xMinorTicks);

# Set population cutoff values.
# Variables are named white, black and asian for European, East Asian and Africans
# so that they have the same lengths
my $whiteCut = $params->{"ecut"};
my $blackCut = $params->{"fcut"};
my $asianCut = $params->{"acut"};
my $oHispCut = $params->{"ohcut"}; # Hispanics and 'Other'
my $bHispCut = $params->{"fhcut"}; # Hispanics and African Americans

# The default cutoff values (%)
$whiteCut = 87 unless ($whiteCut);
$blackCut = 95 unless ($blackCut);
$asianCut = 95 unless ($asianCut);
$oHispCut = 13 unless ($oHispCut);
$bHispCut = 40 unless ($bHispCut);

my $meanSasx = 1.69;
my $meanSasy = 0.0617558659013189;
my $sdSasx = 0.0146156873670083;
my $sdSasy = 0.00500353360060671;

my $numSasSds = 4;
my $sasCutBasey = $meanSasy - $numSasSds * $sdSasy;
my $sasCutaVal = 5;

my $asianCutBasex = 1.73;
my $meanAsiany = 0;
my $asianCutaVal = 30;

# Do some (not all) sanity checks
my $cutErr = "";
$cutErr .= "\tecut < 0\n"  if ($whiteCut < 0);
$cutErr .= "\tohcut < 0\n" if ($oHispCut < 0);
$cutErr .= "\tfhcut < 0\n" if ($bHispCut < 0);

$cutErr .= "\tecut > 100\n"  if ($whiteCut > 100);
$cutErr .= "\tacut > 100\n"  if ($asianCut > 100);
$cutErr .= "\tfcut > 100\n"  if ($blackCut > 100);
$cutErr .= "\tohcut > 100\n" if ($oHispCut > 100);
$cutErr .= "\tfhcut > 100\n" if ($bHispCut > 100);

$cutErr .= "\tecut + ohcut < 100\n" if ($whiteCut + $oHispCut < 100);
$cutErr .= "\tfcut + ohcut < 100\n" if ($blackCut + $oHispCut < 100 && $blackCut > 0);
$cutErr .= "\tacut + ohcut < 100\n" if ($asianCut + $oHispCut < 100 && $asianCut > 0);
$cutErr .= "\tfhcut > fcut\n"       if ($bHispCut > $blackCut && $blackCut > 0);
$cutErr .= "\tfhcut + ecut < 100\n" if ($bHispCut < 100 - $whiteCut);

if ($cutErr) {
    print "\nERROR in population cutoff line settings:\n$cutErr\n";
    exit;
}

# The cutoff line will not be drawn if they are set to be negative
$blackCut = 200 if ($blackCut < 0);
$asianCut = 200 if ($asianCut < 0);

# When calculating populations, use cutoff values in (0, 1)
$whiteCut /= 100.0;
$blackCut /= 100.0;
$asianCut /= 100.0;
$bHispCut /= 100.0;
$oHispCut /= 100.0;

my $sasMexCut  = 0.04;  # For separating South Asians from Hispanics using a different score
my $cutLineLen = 0.07;  # Length of the cutoff lines outside the triangle

# Decide to show expected areas for which populations
my @showPredAreas = ();
@showPredAreas = split /\,/, $showAreas;
my %predAreaHash = ();
for my $area (@showPredAreas) {
    $predAreaHash{$area} = 1;
}
my $whitePs = $predAreaHash{1} ? 1 : 0;
my $ghanaPs = $predAreaHash{2} ? 1 : 0;
my $chinaPs = $predAreaHash{3} ? 1 : 0;
my $blackPs = $predAreaHash{4} ? 1 : 0;
my $hisp1Ps = $predAreaHash{5} ? 1 : 0;
my $mexicPs = $predAreaHash{6} ? 1 : 0;
my $asianPs = $predAreaHash{7} ? 1 : 0;
my $indiaPs = $predAreaHash{8} ? 1 : 0;

if ($showGd4) {
    $rotx = 0;
    $whitePs = 0;
    $ghanaPs = 0;
    $chinaPs = 0;
    $hisp1Ps = 0;
    $mexicPs = 0;
    $asianPs = 0;
    $indiaPs = 0;
}

#-------------------------------- Read subject populations from file generated by GRAF-pop  -------------------------------#

my @selRaceVals = split /\,/, $selRaces; # The race IDs selected to display on the graph
my ($minFpSnps, $maxFpSnps, $totFpSnps) = (10000, 0, 0);

my $unkRace = "NOT REPORTED";
my $sbjNo = 1;
my %sbjPopScores = ();
my %allSbjs = ();
open FILE, $inGrafFile or die "Couldn't open $inGrafFile!\n";
my $header = <FILE>;
if ($header !~ /^DS No\.\tSample\t\#SNPs\tGD1\tGD2\tGD3\tGD4/) {
    print "\nERROR: Invalid input file.  Expected following columns:\n" .
          "\tDS No.\n\tSample\n\t\#SNPs\n\tGD1\n\tGD2\n\tGD3\n\tGD4\n\n";
    exit;
}
while(<FILE>) {
    chomp;
    my ($dsNo, $sbj, $numSnps, $xVal, $yVal, $zVal, $gd4, $fPct, $ePct, $aPct, $black, $white, $asian, $mexican, $indPak) = split /\t/, $_;
    if ($sbj && $numSnps && $black && $white && $asian) {
	my %info = (subject => $sbj, race => $unkRace, snps => $numSnps,
		    x => $xVal, y => $yVal, z => $zVal, gd4 => $gd4, fPct => $fPct, ePct => $ePct, aPct => $aPct,
		    white => $white, black => $black, asian => $asian, indPak => $indPak, mexican => $mexican);
	if ($numSnps >= $minSnps && $numSnps <= $maxSnps) {
	    $sbjPopScores{$sbjNo} = \%info;
	    $allSbjs{$sbj} = 1;
	}

	$minFpSnps = $numSnps if ($numSnps < $minFpSnps);
	$maxFpSnps = $numSnps if ($numSnps > $maxFpSnps);
	$totFpSnps += $numSnps;

	$sbjNo++;
    }
}
close FILE;
my $totSbjs = $sbjNo - 1;
my $numSbjs = keys %sbjPopScores;
my $meanFpSnps = $totFpSnps * 1.0 / $totSbjs;

if ($totSbjs < 1) {
    print "\nERROR: No subject found in $inGrafFile!\n";
    exit;
}

print "Found $totSbjs subjects with population scores in file $inGrafFile\n";
printf("Mean FP SNPs: %d (%d - %d)\n", $meanFpSnps, $minFpSnps, $maxFpSnps);
if ($minSnps > 0 || $maxSnps < 10000) {
    if ($numSbjs > 0) {
	print "  $numSbjs subjects have genotyped FP SNPs between $minSnps and $maxSnps.\n";
    }
    else {
	print "\nNo subjects with genotype FP SNPs between $minSnps and $maxSnps found in the input file.\n\n";
	exit;
    }
}

# If file exists, read self-repoted ancestries from the input file
my $hasRace = 0;
my %sbjRaces = (); # Ancestry of each subject
my %allRaces = ();
if ($raceFile) {
    unless (-e $raceFile) {
	print "\nERROR: didn't find subject race file $raceFile!\n";
	exit;
    }

    open RACEFILE, $raceFile or die "Couldn't open $raceFile!\n";
    while (<RACEFILE>) {
	chomp;
	my ($sbj, $race) = split /\t/, $_;
	$race =~ s/\s*$//;
	$race = $1 if ($race =~ /^\s*\"(.+)\"\s*$/);
	if ($sbj && $allSbjs{$sbj}) {
	    $race = $unkRace if (!$race || $race !~ /\S/);
	    $hasRace = 1 if ($race);
	    $sbjRaces{$sbj} = $race;
	    $allRaces{$race} = 1;
	}
    }
    close RACEFILE;

    my $numSbjs = keys %sbjRaces;
    my $numRaces = keys %allRaces;
    if ($numRaces > 0) {
	print "\nRead $numRaces ancestries from $numSbjs subjects in $raceFile\n";
    }
    else {
	print "\nNOTE: No ancestries found in $raceFile.\n";
    }
}

# Get subject counts for different self-reported ancestries
my %raceSbjCnts = ();
my $minRaceSbjs = 0;

foreach my $sbjNo (keys %sbjPopScores) {
    my %info = %{$sbjPopScores{$sbjNo}};
    my $sbj = $sbjPopScores{$sbjNo}->{subject};
    my $race = $unkRace;
    $race = $sbjRaces{$sbj} if ($sbjRaces{$sbj});

    if ($raceSbjCnts{$race}) {
	$raceSbjCnts{$race}++;
    }
    else {
	$raceSbjCnts{$race} = 1;
    }
}

# Sort self-reported ancestries by subject counts
my @sortRaces = ();
foreach my $race (sort {$raceSbjCnts{$b} <=> $raceSbjCnts{$a}} keys %raceSbjCnts) {
    my $cnt = $raceSbjCnts{$race};
    push @sortRaces, $race if ($cnt > $minRaceSbjs);
}

if ($selRaces !~ /\S/) {
    for my $raceNo (0 .. $#sortRaces) {
	my $raceId = $raceNo + 1;
	my $race = $sortRaces[$raceNo];
	push @selRaceVals, $raceId;
	$selRaces .= "$raceId,";
	last if ($raceNo > 8);
    }
    $selRaces =~ s/\,\s*$//;
}

# Create a hash to find ancestry name from the selected ancestry ID.
my %selRaceNames = ();
my %selRaceIds = ();
for my $val (@selRaceVals) {
    my $id = $val - 1;
    my $race = $sortRaces[$id];
    $selRaceIds{$id} = 1;
    $selRaceNames{$race} = $id;
}
my $numSelRaces = keys %selRaceIds;


#------------------------------------ Set parameters for the graph -----------------------------------#

my $leftEdge = 20;        # The left edge of the graph

# Width and height of the plot (the box excluding header and footer)
my $graphWidth  = 600;
my $graphHeight = 600;
$graphWidth  = $gWidth  if ($gWidth);
$graphHeight = $gHeight if ($gHeight);
$graphHeight = $graphWidth * $yRange * 1.0 / $xRange;

my $rightEdge = 20;       # The right edge
my $gyTop = 10;           # The top edge
my $bottomEdge = 50;      # The bottom edge
my $graphTitleHt = 20;    # Height of the graph title

my $majorTickLen = 6;  # major tick length (in pixels)
my $minorTickLen = 3;  # minor tick length

my $mbFontWidth  = gdMediumBoldFont->width;
my $mbFontHeight = gdMediumBoldFont->height;
my $mlFontWidth  = gdLargeFont->width;

my $axisLabelGap = 20;  # Distance between tick labels and axis labels
my $tickGap = 5;        # Gap between the tick and tick label

# Left edge of the plot (the rectangular box)
my $gxLeft = $leftEdge + $mbFontHeight + $axisLabelGap + $mbFontWidth*2 + $tickGap + $majorTickLen;
my $gxRight  = $gxLeft + $graphWidth;
my $gyBottom = $gyTop + $graphHeight;

# Width and height of the whole graph
my $imageWidth = $gxLeft + $graphWidth + $rightEdge;
my $imageHeight = $gyTop + $graphHeight + $graphTitleHt + $bottomEdge;

my $topGap = 100;


#------------------------------------------- Set colors  ----------------------------------------------#
my $cgi = new CGI;
my $im = new GD::Image($imageWidth, $imageHeight);

my $white   = $im->colorAllocate(255, 255, 255);
my $red     = $im->colorAllocate(255,   0,   0);
my $maroon  = $im->colorAllocate(255,   0,   0);
my $green   = $im->colorAllocate(255, 128,   0); # not green anymore
my $blue    = $im->colorAllocate(  0,   0, 255);
my $navy    = $im->colorAllocate(  0,   0, 128);
my $black   = $im->colorAllocate(  0,   0,   0);
my $gray    = $im->colorAllocate(128, 128, 128);
my $yellow  = $im->colorAllocate(255, 255,   0);
my $olive   = $im->colorAllocate(128, 128,   0);
my $purple  = $im->colorAllocate(128,   0, 128);
my $magenta = $im->colorAllocate(255,   0, 255);
my $orange  = $im->colorAllocate(255, 165,   0);
my $cyan    = $im->colorAllocate(  0, 255, 255);
my $teal    = $im->colorAllocate(  0, 128, 128);
my $gold    = $im->colorAllocate(204, 153,  80);
my $expAreaColor = $im->colorAllocate(  0,150,150);

# Colors for self-reported ancestries
my $numShowRaces = 10;
my @raceColors = ();
for my $raceNo (1 .. $numShowRaces) {
    push @raceColors, "";
}
$raceColors[0] = $yellow;
$raceColors[1] = $blue;
$raceColors[2] = $red;
$raceColors[3] = $olive;
$raceColors[4] = $purple;
$raceColors[5] = $cyan;
$raceColors[6] = $green;
$raceColors[7] = $teal;
$raceColors[8] = $navy;
$raceColors[9] = $magenta;
$raceColors[10] = $gold;


# Create a hash to find color for a self-reported ancestry
my %raceColorIds = ();
my $selRaceNo = 1;
foreach my $raceId (0 .. $#sortRaces) {
    my $race = $sortRaces[$raceId];

    my $colorNo = 0;
    if (defined $selRaceIds{$raceId}) {
	$colorNo = $selRaceNo;
	$selRaceNo++;
    }
    $colorNo = $colorNo % $numShowRaces;
    $raceColorIds{$race} = $colorNo;
}

my $raceNo = 0;
my %raceIdNos = ();
foreach my $raceId (sort {$a <=> $b} keys %selRaceIds) {
    $raceIdNos{$raceId} = $raceNo;
    $raceNo++;
}



#------------------------------------------- Plot the graph ----------------------------------------------#

# Save scores, colors, self-reported ancestries, etc. to global arrays
my @allSbjPvalues = ();
my @allSbjColors = ();
my @allSbjColorIds = ();
my @allSbjIds = ();
my @allSbjRaces = ();
my @allSbjSnps = ();
my @popMeanPvalues = ();
my @sbjIndPaks = ();
my @sbjMexicans = ();
GetSubjectPositions();
my $numShowSbjs = @allSbjPvalues;
print "Read scores for $numShowSbjs subjects\n";

#my @aCtr = (1.8911, 1.1641, 0.8219, 1);
#my @wCtr = (1.6290, 1.4420, 0.8219, 1);
#my @bCtr = (1.2963, 1.1641, 0.8219, 1);

my @aCtr = (1.8972, 1.1673, 0.8227, 1);
my @wCtr = (1.6286, 1.4513, 0.8227, 1);
my @bCtr = (1.2993, 1.1673, 0.8227, 1);

# The x, y, z coordinates for all subjects to be plotted
my @allPlotPts = ();
for my $i (0 .. $#allSbjPvalues) {
    push @allPlotPts, $allSbjPvalues[$i];
}
# The x, y, z coordinates of the three vertices of the triangle
push @allPlotPts, \@aCtr;
push @allPlotPts, \@wCtr;
push @allPlotPts, \@bCtr;

# Move and rotate the values to project them on the EFA plan
#TransformAllSubjects();

# Locations of the 3 vertices of the triangle after transformation.
# a: Asian; w: White (European); b: Black (African)
my ($ax0, $ay0, $az0, $junk1) = @{$allPlotPts[$numShowSbjs+0]};
my ($wx0, $wy0, $wz0, $junk2) = @{$allPlotPts[$numShowSbjs+1]};
my ($bx0, $by0, $bz0, $junk3) = @{$allPlotPts[$numShowSbjs+2]};

# The angles (in degree) of the two sides of the triangle
my $leftEdgeAng = atan2(($wx0-$bx0)*1.0, ($wy0-$by0)) * 180 / $pi;  # The Black/White side
my $rightEdgeAng = atan2(($wx0-$ax0)*1.0, ($wy0-$ay0)) * 180 / $pi; # The Asian/White side

# The three vertices before rotation on the sides
my @aPt0 = ($ax0, $ay0, $az0, $junk1);
my @wPt0 = ($wx0, $wy0, $wz0, $junk2);
my @bPt0 = ($bx0, $by0, $bz0, $junk3);

# The three vertices after rotation on the left side
my @aPt1 = RotateOnLeftEdge(\@aPt0);
my @wPt1 = RotateOnLeftEdge(\@wPt0);
my @bPt1 = RotateOnLeftEdge(\@bPt0);

my ($ax1, $ay1, $az1, $j11) = @aPt1;
my ($wx1, $wy1, $wz1, $j12) = @wPt1;
my ($bx1, $by1, $bz1, $j13) = @bPt1;

# The three vertices after rotation on the right side
my @aPt2 = RotateOnRightEdge(\@aPt0);
my @wPt2 = RotateOnRightEdge(\@wPt0);
my @bPt2 = RotateOnRightEdge(\@bPt0);

my ($ax2, $ay2, $az2, $j21) = @aPt2;
my ($wx2, $wy2, $wz2, $j22) = @wPt2;
my ($bx2, $by2, $bz2, $j23) = @bPt2;

my $cutoffLineColor = $gold;

my %genoPops = (
		1 => "European",
		2 => "African",
		3 => "East Asian",
		4 => "African American",
		5 => "Hispanic1",
		6 => "Hispanic2",
		7 => "Other Asian or Pacific Islander",
		8 => "South Asian",
		9 => "Other",
	       );

$genoPops{"4"} = "African or African American" if ($blackCut > 1); #Geno IDs 3 and 4 are merged into one when blackCut is greater than 1
$genoPops{"7"} = "Asian or Pacific Islander"    if ($asianCut > 1); #Geno IDs 3 and 7 are merged into one when asianCut is greater than 1

if ($rotx) {
    RotateSubjectsOnx($rotx);
}

# Plot the scores or output the subjects
if ($showSbjs) {
    ShowSubjects($outputFile); # Save subjects to a file
}
else {
    PlotPopulations($outputFile);      # Plot score to a png file
}

my @finishTime = gettimeofday;
ShowTimeDifference(\@startTime, \@finishTime);



#------------------------------------------- Subroutines ----------------------------------------------#

#
# Plot graph to show positions of each subject
#
sub PlotPopulations
{
    my $outPngFile = shift;

#    print $cgi->header( {-type=>'image/png' });

    open(IMG, ">$outPngFile") or die $!;
    binmode IMG;

    my $xLbl = "GD1";
    my $yLbl = "GD2";
    if ($rotx) {
	$yLbl = "GD2, GD3: rotated by $rotx degree";
    }
    if ($showGd4) {
	$yLbl = "GD4";
    }

    PlotAxes($gyTop, $xLbl, $yLbl);

    # Legends to show what color represents which self-reported ancestry
    my $lgdx = $gxLeft + 10;
    my $lgdy = $gyTop + 10;
    my $lgdSide = 10;
    my $lgdGap = $lgdSide + 10;

    my $snpLgdGap = 105; # Room to show "FP SNPs/Sbj"
    $snpLgdGap = 85 if ($graphWidth < 401);

    my $hasOther = 0;
    my $numSelSbjs = 0;
    my $numOthSbjs = 0;
    my $raceNo = 0;
    foreach my $raceId (0 .. $#sortRaces) {
	my $raceNum = $raceId + 1;
	my $race = $sortRaces[$raceId];
	my $cnt = $raceSbjCnts{$race};
	my $dispRace = TruncateDisplayRace($race, $cnt, $lgdSide, $snpLgdGap);
	my $colorNo = $raceColorIds{$race};
	my $color = $raceColors[$colorNo];
	$color = $black unless ($hasRace);
	if ($colorNo > 0) {
	    $im->filledRectangle($lgdx, $lgdy, $lgdx + $lgdSide, $lgdy + $lgdSide, $color);
	    $im->string(gdMediumBoldFont, $lgdx + $lgdSide + 10, $lgdy, "$dispRace ($cnt)", $black);
	    $lgdy += $lgdGap;
	    $numSelSbjs += $cnt;
	}
	else {
	    $hasOther = 1;
	    $numOthSbjs += $cnt;
	}

	$raceNo++;
    }

    if ($hasOther) {
	my $color = $raceColors[0];
	$im->filledRectangle($lgdx, $lgdy, $lgdx + $lgdSide, $lgdy + $lgdSide, $color);
	$im->string(gdMediumBoldFont, $lgdx + $lgdSide + 10, $lgdy, "unselected ($numOthSbjs)", $black);
	$lgdy += $lgdGap;
    }

    $im->string(gdMediumBoldFont, $lgdx, $lgdy, "Total $numSbjs subjects", $black);

    # Show min and max FP SNPs with genotypes per subject
    my $snpLgdx = $gxLeft + $graphWidth - $snpLgdGap;
    my $snpLgdy = $gyTop + 10;
    my $snpLgdGap = $mbFontHeight * 1.5;

    $im->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, "FP SNPs/Sbj:", $black);

    $snpLgdx += $mbFontWidth;
    $snpLgdy += $snpLgdGap;
    my $minSnpText = sprintf(" Min %5d", $minFpSnps);
    $im->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $minSnpText, $black);

    $snpLgdy += $snpLgdGap;
    my $maxSnpText = sprintf(" Max %5d", $maxFpSnps);
    $im->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $maxSnpText, $black);

    $snpLgdy += $snpLgdGap;
    my $meanSnpText = sprintf("Mean %5.0f", $meanFpSnps);
    $im->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $meanSnpText, $black);

    PlotAllSubjects();

    print IMG $im->png;
    print "Graph saved to $outPngFile\n\n";

    exit;
}

#
# Check if there is enough room to show the self-reported ancestry.
# If not, truncate the ancestry.
#
sub TruncateDisplayRace
{
    my ($race, $cnt, $lgdSide, $snpLgdGap) = @_;

    my $dispRace = $race;

    my $raceLen = length($race);
    my $cntLen = length($cnt);
    my $lblLen = $raceLen + $cntLen + 4;

    my $extraLen = 10 + $lgdSide + $lblLen * $mbFontWidth + $snpLgdGap - $graphWidth;
    if ($extraLen > 0) {
	my $cutChars = int($extraLen * 1.0 / $mbFontWidth) + 1;
	$dispRace = substr($race, $raceLen);
	$dispRace = substr($race, 0, $raceLen - $cutChars - 3);
	$dispRace .= "...";
    }

    return $dispRace;
}

#
# Given a position (x, y) of a subject, calculate the ancestry proportions for each main population
#
sub GetPopulationComponent
{
    my ($xVal, $yVal) = @_;

    my @chkPt0 = ($xVal, $yVal, 1, 1);
    my ($xVal1, $yVal1, $junk1, $junk2) = RotateOnLeftEdge(\@chkPt0);
    my ($xVal2, $yVal2, $junk3, $junk4) = RotateOnRightEdge(\@chkPt0);

    # Height from each vertex to the opposite side
    my $Hw = $wy0 - $by0;
    my $Hb = $ax0 - $bx2;
    my $Ha = $ax1 - $bx0;

    # Distance from this subject to the above side
    my $hw = $yVal  - $by0;
    my $hb = $ax0   - $xVal2;
    my $ha = $xVal1 - $bx0;

    # Set value to 0 if subject is not on the same side as the vertex is
    $hw = 0 if ($hw < 0);
    $hb = 0 if ($hb < 0);
    $ha = 0 if ($ha < 0);

    my $pSum = $hw/$Hw + $ha/$Ha + $hb/$Hb;
    my $wPct = $hw/($Hw * $pSum);
    my $bPct = $hb/($Hb * $pSum);
    my $aPct = $ha/($Ha * $pSum);

    return ($bPct, $wPct, $aPct);
}

#
# Plot cutoff lines to separate SAS, ASN from other populations, on GD4 vs. GD1 graph
#
sub PlotSasCutoffLines
{
    my $halfWid = 0.1;
    my $numPts = 1000;

    my ($x1, $y1) = (0, 0);
    for my $i (0 .. $numPts) {
	my $dx = ($i*2 - $numPts)/$numPts * $halfWid;
	my $xVal = $meanSasx + $dx;
	my $yVal = $sasCutBasey + $sasCutaVal * $dx * $dx;

	my $x2 = $gxLeft   + ($xVal - $xMin) * $graphWidth  / ($xMax - $xMin);
	my $y2 = $gyBottom - ($yVal - $yMin) * $graphHeight / ($yMax - $yMin);
	if ($i > 0) {
	    $im->line($x1, $y1, $x2, $y2, $cutoffLineColor);
	}
	$x1 = $x2;
	$y1 = $y2;
    }

    my $halfHt = 0.085;
    ($x1, $y1) = (0, 0);
    for my $i (0 .. $numPts) {
	my $dy = ($i*2 - $numPts)/$numPts * $halfHt;
	my $yVal = $meanAsiany + $dy;
	my $xVal = $asianCutBasex + $asianCutaVal * $dy * $dy;

	my $y2 = $gyBottom - ($yVal - $yMin) * $graphHeight / ($yMax - $yMin);
	my $x2 = $gxLeft   + ($xVal - $xMin) * $graphWidth  / ($xMax - $xMin);
	if ($i > 0) {
	    $im->line($x1, $y1, $x2, $y2, $cutoffLineColor);
	}
	$x1 = $x2;
	$y1 = $y2;
    }
}


#
# Plot cutoff lines to divide different populations
#
sub PlotCutoffLines
{
    # East Asian line inside the triangle
    my $axVal1 = $ax0 - ($ax0 - $bx0) * (1 - $asianCut);
    my $ayVal1 = $ay0;
    my $axVal2 = $ax0 - ($ax0 - $wx0) * (1 - $asianCut);
    my $ayVal2 = $ay0 + ($wy0 - $ay0) * (1 - $asianCut);
    PlotOneCutoffLine($axVal1, $ayVal1, $axVal2, $ayVal2) if ($asianCut < 1);
#    $im->string(gdMediumBoldFont, 40, 260, "ax0 $ax0 ay0 $ay0 wx0 $wx0 wy0 $wy0 cut $asianCut", $black);

    # African line inside the triangle
    my $bxVal1 = $bx0 + ($ax0 - $bx0) * (1 - $blackCut);
    my $byVal1 = $by0;
    my $bxVal2 = $bx0 + ($wx0 - $bx0) * (1 - $blackCut);
    my $byVal2 = $by0 + ($wy0 - $by0) * (1 - $blackCut);
    PlotOneCutoffLine($bxVal1, $byVal1, $bxVal2, $byVal2) if ($blackCut < 1);

    # European line inside the triangle
    my $wxVal1 = $wx0 - ($wx0 - $bx0) * (1 - $whiteCut);
    my $wyVal1 = $wy0 - ($wy0 - $by0) * (1 - $whiteCut);
    my $wxVal2 = $wx0 + ($ax0 - $wx0) * (1 - $whiteCut);
    my $wyVal2 = $wyVal1;
    PlotOneCutoffLine($wxVal1, $wyVal1, $wxVal2, $wyVal2);

    # Cutoff lines outside the triangle, extended from the above three lines
    PlotOuterCutoff($axVal1, $ayVal1, 1);
    PlotOuterCutoff($bxVal1, $byVal1, 1);
    PlotOuterCutoff($axVal2, $ayVal2, 2);
    PlotOuterCutoff($wxVal1, $wyVal1, 3);
    PlotOuterCutoff($wxVal2, $wyVal2, 2);
    PlotOuterCutoff($bxVal2, $byVal2, 3);

    # Cutoff lines outside the triangle, extended from the other lines inside the triangle
    my $bxVal3 = $ax0 - ($ax0 - $bx0) * (1 - $oHispCut);
    my $bxVal4 = $ax0 - ($ax0 - $bx0) * $oHispCut;
    my $byVal3 = $by0;
    my $byVal4 = $by0;
    my $wxVal3 = $wx0 - ($wx0 - $bx0) * $bHispCut;
    my $wyVal3 = $wy0 - ($wy0 - $by0) * $bHispCut;

    PlotOuterCutoff($bxVal3, $byVal3, 1);
    PlotOuterCutoff($bxVal4, $byVal4, 1);
#    PlotOuterCutoff($axVal3, $ayVal3, 2);
    PlotOuterCutoff($wxVal3, $wyVal3, 3);

    # Cutoff lines outside the triangle, separating Hispanics from other populations
    my $obxVal1 = $bxVal3;
    my $obyVal1 = $by0;
    my $obxVal2 = $wx0 + ($ax0 - $wx0) * $oHispCut;
    my $obyVal2 = $wy0 - ($wy0 - $ay0) * $oHispCut;

    my $oaxVal1 = $bxVal4;
    my $oayVal1 = $by0;
    my $oaxVal2 = $wx0 - ($wx0 - $bx0) * $oHispCut;
    my $oayVal2 = $wy0 - ($wy0 - $by0) * $oHispCut;

    my $obhxVal = $obxVal2 - ($wx0 - $bx0) * $oHispCut;
    my $obhyVal = $obyVal2 - ($wy0 - $by0) * $oHispCut;
    PlotOneCutoffLine($oaxVal1, $oayVal1, $obhxVal, $obhyVal);
    PlotOneCutoffLine($obxVal1, $obyVal1, $obhxVal, $obhyVal);

    # Cutoff line separating the two Hispanic sub-populations
    my $abmx = ($ax0 + $bx0) / 2;
    my $whxVal = $wx0 - ($wx0 - $abmx) * (1 - $whiteCut);
    my $whyVal = $wy0 - ($wy0 - $by0)  * (1 - $whiteCut);
    PlotOneCutoffLine($whxVal, $whyVal, $obhxVal, $obhyVal);

    # Cutoff lines separating Hispanics from African Americans
    my $bhxVal1 = $wxVal3;
    my $bhyVal1 = $wyVal3;
    my $bhxVal2 = $wxVal3 + ($ax0 - $wx0) * $oHispCut;
    my $bhyVal2 = $wyVal3 - ($wy0 - $ay0) * $oHispCut;
    PlotOneCutoffLine($bhxVal1, $bhyVal1, $bhxVal2, $bhyVal2);
}

#
# Plot cutoff line outside the triangle, starting from a certain point
#
sub PlotOuterCutoff
{
    my ($x1, $y1, $popId) = @_;

    my ($dx, $dy) = (0, 0);

    if ($popId == 1) {
	$dx = $x1 - $wx0;
	$dy = $y1 - $wy0;
    }
    elsif ($popId == 2) {
	$dx = $x1 - $bx0;
	$dy = $y1 - $by0;
    }
    elsif ($popId == 3) {
	$dx = $x1 - $ax0;
	$dy = $y1 - $ay0;
    }

    my $dl = sqrt($dx**2 + $dy**2);
    my $x2 = $x1 + $dx * $cutLineLen / $dl;
    my $y2 = $y1 + $dy * $cutLineLen / $dl;

    PlotOneCutoffLine($x1, $y1, $x2, $y2);
}

#
# Plot one cutoff line given two points
#
sub PlotOneCutoffLine
{
    my ($x1, $y1, $x2, $y2) = @_;

    my $xPos1 = $gxLeft   + ($x1 - $xMin) * $graphWidth  / $xRange;
    my $yPos1 = $gyBottom - ($y1 - $yMin) * $graphHeight / $yRange;
    my $xPos2 = $gxLeft   + ($x2 - $xMin) * $graphWidth  / $xRange;
    my $yPos2 = $gyBottom - ($y2 - $yMin) * $graphHeight / $yRange;

    $im->line($xPos1, $yPos1, $xPos2, $yPos2, $cutoffLineColor);
}

#
# Plot subject dots on the graph
#
sub ShowSubjects
{
    my $outFile = shift;

    my $numShowSbjs = 0;
    my @sbjOutputLines = ();
    my %subPopSbjs = ();
    my $totSbjs = @allSbjPvalues;

    for my $sbjNo (0 .. $#allSbjPvalues) {
	my ($gd1, $gd2, $gd3, $junk) = @{$allPlotPts[$sbjNo]};

	# The two scores separating South Asians from Hispanics
	my $indPak  = $sbjIndPaks[$sbjNo];
	my $mexican = $sbjMexicans[$sbjNo];
	my $gd4 = $mexican - $indPak;

	my $sbjId   = $allSbjIds[$sbjNo];
	my $numSnps = $allSbjSnps[$sbjNo];
	my $race    = $allSbjRaces[$sbjNo];
	my $raceId  = $selRaceNames{$race};

	# Check if the subject is within the selected area or not
	my $isWithin = 1;
	$isWithin = 0 if ($xCutMin && $gd1 < $xCutMin);
	$isWithin = 0 if ($yCutMin && $gd2 < $yCutMin);
	$isWithin = 0 if ($xCutMax && $gd1 > $xCutMax);
	$isWithin = 0 if ($yCutMax && $gd2 > $yCutMax);

	# Get ancestry proportions and calculate genotype populations based on the proportions
	my ($bPct, $wPct, $aPct) = GetPopulationComponent($gd1, $gd2);
	my $genoPopId = 0;
	if ($numSnps > 1000) {
	    $genoPopId = GetGenoPopId($gd1, $gd4, $wPct, $bPct, $aPct);
	}

	# Decide whether to include the subject or not based on the selected area
	my $showData = 1;
	$showData = 0 if (($isBeyond && $isWithin) || (!$isBeyond && !$isWithin));
	if ($isSelAll || (defined $raceId && $showData) ) {
	    my $race  = $sbjRaces{$sbjId};
	    my $showGd1  = sprintf("%4.2f", $gd1);
	    my $showGd2  = sprintf("%4.2f", $gd2);
	    my $showGd3  = sprintf("%4.2f", $gd3);
	    my $showGd4  = sprintf("%4.2f", $gd4);
	    my $showbPct = sprintf("%5.1f", $bPct * 100);
	    my $showwPct = sprintf("%5.1f", $wPct * 100);
	    my $showaPct = sprintf("%5.1f", $aPct * 100);
	    my $genoPop  = $genoPops{$genoPopId};
	    push @sbjOutputLines,
	        "$sbjId\t$race\t$showGd1\t$showGd2\t$showGd3\t$showGd4\t$showbPct\t$showwPct\t$showaPct\t$genoPopId\t$genoPop";
	    if ($subPopSbjs{$genoPopId}) {
		$subPopSbjs{$genoPopId}++;
	    }
	    else {
		$subPopSbjs{$genoPopId} = 1;
	    }
	    $numShowSbjs++;
	}
    }

    # Save subjects to the output file
    if ($numShowSbjs > 0) {
	open FILE, ">$outFile" or die "\nERROR: Couldn't open $outFile for writing!\n";
	print FILE "Subject\tSelf-reported ancestry\tGD1\tGD2\tGD3\tGD4\tP_f (%)\tP_e (%)\tP_a (%)\tPopID\tComputed population\n";
	for my $line (@sbjOutputLines) {
	    print FILE "$line\n";
	}
	close FILE;
	print "\nTotal $totSbjs subjects. $numShowSbjs subjects were selected and saved to $outFile\n\n";

	print "\tPopID\t#Subjs\tPopulation\n";
	foreach my $genoPopId (sort {$subPopSbjs{$b} <=> $subPopSbjs{$a}} keys %subPopSbjs) {
	    my $genoPop = $genoPops{$genoPopId};
	    my $numSbjs = $subPopSbjs{$genoPopId};
	    printf("\t%d\t%6d\t%s\n", $genoPopId, $numSbjs, $genoPop);
	}
	print "\n";
    }
    else {
	print "\nNo subjects are found in the specified area.\n\n";
    }
}

#
# Rotate subject around the African vertex so that the EF line is parallel to the y-axis
#
sub RotateOnLeftEdge
{
    my $p1Ref = shift;
    my @p1 = @$p1Ref;

    my @p2 = MovePoint3D(\@p1, -$bx0, -$by0, -$bz0);
    my @p3 = RotatePoint3D(\@p2, $leftEdgeAng, "z");
    my @p4 = MovePoint3D(\@p3, $bx0, $by0, $bz0);

    return @p4;
}

#
# Rotate subject around the Asian vertex so that the EA line is parallel to the y-axis
#
sub RotateOnRightEdge
{
    my $p1Ref = shift;
    my @p1 = @$p1Ref;

    my @p2 = MovePoint3D(\@p1, -$ax0, -$ay0, -$az0);
    my @p3 = RotatePoint3D(\@p2, $rightEdgeAng, "z");
    my @p4 = MovePoint3D(\@p3, $ax0, $ay0, $az0);

    return @p4;
}

#
# Move and rotate all subjects to project them on the EFA plane, and the FA line is parallel to the x-axis
#
sub TransformAllSubjects
{
    my ($bx0, $by0, $bz0, $junk1) = @{$allPlotPts[$numShowSbjs+2]};
    my ($ax0, $ay0, $az0, $junk2) = @{$allPlotPts[$numShowSbjs+0]};
    my ($wx0, $wy0, $wz0, $junk2) = @{$allPlotPts[$numShowSbjs+1]};

    MoveShape3D(-$bx0, -$by0, -$bz0);

    ($ax0, $ay0, $az0, $junk2) = @{$allPlotPts[$numShowSbjs+0]};
    my $theta1 = atan2($ay0*1.0, $ax0) * -180 / $pi;
    RotateShape3D($theta1, "z");

    ($ax0, $ay0, $az0, $junk2) = @{$allPlotPts[$numShowSbjs+0]};
    my $theta2 = atan2($az0*1.0, $ax0) * 180 / $pi;

    RotateShape3D($theta2, "y");
    ($ax0, $ay0, $az0, $junk2) = @{$allPlotPts[$numShowSbjs+0]};

    ($wx0, $wy0, $wz0, $junk2) = @{$allPlotPts[$numShowSbjs+1]};
    my $theta3 = atan2($wy0*1.0, ($wz0+0)) * 180 / $pi;
    $theta3 -= 90;
    $theta3 += $rotx if ($rotx);

    RotateShape3D($theta3, "x");
    MoveShape3D($bx0, $by0, $bz0);
}

#
# Move and rotate all subjects to project them on the EFA plane, and the FA line is parallel to the x-axis
#
sub RotateSubjectsOnx
{
    my $theta = shift;

    my ($bx0, $by0, $bz0, $junk1) = @{$allPlotPts[$numShowSbjs+2]};
    my ($ax0, $ay0, $az0, $junk2) = @{$allPlotPts[$numShowSbjs+0]};
    my ($wx0, $wy0, $wz0, $junk2) = @{$allPlotPts[$numShowSbjs+1]};

    MoveShape3D(-$bx0, -$by0, -$bz0);
    RotateShape3D(-$theta, "x");
    MoveShape3D($bx0, $by0, $bz0);
}

#
# Plot a filled circle
#
sub PlotOneDot
{
    my ($xVal, $yVal, $size, $color) = @_;

    my $x = $gxLeft   + ($xVal - $xMin) * $graphWidth  * 1.0 / $xRange;
    my $y = $gyBottom - ($yVal - $yMin) * $graphHeight * 1.0 / $yRange;
    $x = int($x + 0.5);
    $y = int($y + 0.5);

    # Acoid using fill() since:
    # 1. It is hard to find spots to fill when circles overlap
    # 2. It may overspill, especially when the dots are close to the edges
    $im->setPixel($x, $y, $color);

    if ($size > 1) {
	$im->arc($x, $y, $size, $size, 0, 360, $color);
	if ($size > 4) {
	    $im->setPixel($x+1, $y, $color);
	    $im->setPixel($x-1, $y, $color);
	    $im->setPixel($x, $y+1, $color);
	    $im->setPixel($x, $y-1, $color);
	}
	if ($size > 5) {
	    for my $iSize (5 .. $size-1) {
		$im->arc($x, $y, $iSize, $iSize, 0, 360, $color);
	    }
	}
    }
}

#
# Plot the ellipse showing the expected area
#
sub PlotEllipse
{
    my ($xc, $yc, $a, $b, $rXc, $rYc, $rotAng, $color) = @_;

    my $yPos = $gyTop;
    my $yBottom = $gyBottom;

    my $numPts = 360;

    my $dotSize = 8;
    my $xcp = $gxLeft  + ($xc - $xMin) * $graphWidth  / ($xMax - $xMin);
    my $ycp = $yBottom - ($yc - $yMin) * $graphHeight / ($yMax - $yMin);
#    $im->arc($xcp, $ycp, $dotSize, $dotSize, 0, 360, $color);

    my $x1 = $gxLeft  + ($xc + $a - $xMin) * $graphWidth  / ($xMax - $xMin);
    my $y1 = $yBottom - ($yc - $yMin) * $graphHeight / ($yMax - $yMin);

    for my $i (0 .. $numPts) {
	my $beta = $i * 2 * $pi / $numPts;
	my $x = $a * cos($beta);
	my $y = $b * sin($beta);

	my $x2 = $gxLeft  + ($xc + $x - $xMin) * $graphWidth  / ($xMax - $xMin);
	my $y2 = $yBottom - ($yc + $y - $yMin) * $graphHeight / ($yMax - $yMin);

	if ($i > 0) {
#	    $im->line($x1, $y1, $x2, $y2, $color);
	}

	$x1 = $x2;
	$y1 = $y2;
    }

    my $dx = -1 * $rXc;
    my $dy = -1 * $rYc;

    my ($xcr, $ycr) = TransformPoint2D($xc, $yc, $dx, $dy, $rotAng);

    $xcp = $gxLeft  + ($xcr - $xMin) * $graphWidth  / ($xMax - $xMin);
    $ycp = $yBottom - ($ycr - $yMin) * $graphHeight / ($yMax - $yMin);
#    $im->arc($xcp, $ycp, $dotSize, $dotSize, 0, 360, $black);

    my ($xv1, $yv1) = TransformPoint2D($xc+$a, $yc, $dx, $dy, $rotAng);
    $x1 = $gxLeft   + ($xv1 - $xMin) * $graphWidth  / ($xMax - $xMin);
    $y1 = $yBottom - ($yv1 - $yMin) * $graphHeight / ($yMax - $yMin);

    for my $i (0 .. $numPts) {
	my $beta = $i * 2 * $pi / $numPts;
	my $x = $a * cos($beta);
	my $y = $b * sin($beta);

	my ($xv2, $yv2) = TransformPoint2D($xc+$x, $yc+$y, $dx, $dy, $rotAng);
	my $x2 = $gxLeft   + ($xv2 - $xMin) * $graphWidth  / ($xMax - $xMin);
	my $y2 = $yBottom - ($yv2 - $yMin) * $graphHeight / ($yMax - $yMin);

	if ($i > 0) {
	    $im->line($x1, $y1, $x2, $y2, $color);
	}

	$x1 = $x2;
	$y1 = $y2;
    }
}

#
# Plot population labels for the three vertices of the triangle
#
sub PlotVertexLabel
{
    my ($vLbl, $vType, $vLbly, $extraGap) = @_;

    my $lblLen = length($vLbl);
    my ($xVal, $yVal, $zVal, $junk) = @{$allPlotPts[$numShowSbjs+$vType-1]};
    my $xPos = $gxLeft  + ($xVal - $xMin) * $graphWidth  / ($xMax - $xMin) - $mbFontWidth*$lblLen/2.0;
    my $yPos = $gyBottom - ($vLbly - $yMin) * $graphHeight / ($yMax - $yMin) + $extraGap;
    $yPos -=  $mbFontHeight / 2; # The position of the middle of the label
    $im->string(gdMediumBoldFont, $xPos, $yPos, $vLbl, $black);
}

#
# Plot all subjects, and the triangle
#
sub PlotAllSubjects
{
    # Find the top and bottom positions to show population labels without blocking data points
    my $eMaxy = 0;
    my $fMiny = 10;
    my $aMiny = 10;

    # Z-values determine which subjects are in the front
    my %sbjZvals = ();
    foreach my $sbjNo (0 .. $numShowSbjs-1) {
	my ($xVal, $yVal, $zVal, $junk) = @{$allPlotPts[$sbjNo]};
	$sbjZvals{$sbjNo} = $zVal;

	$eMaxy = $yVal if ($xVal > 1.56 && $xVal < 1.66 && $yVal > $eMaxy);
	$fMiny = $yVal if ($xVal > 1.26 && $xVal < 1.36 && $yVal < $fMiny);
	$aMiny = $yVal if ($xVal > 1.84 && $xVal < 1.96 && $yVal < $aMiny);
    }

    my @sortSbjNos = ();
    foreach my $sbjNo (sort {$sbjZvals{$b} <=> $sbjZvals{$a}} keys %sbjZvals) {
	push @sortSbjNos, $sbjNo;
    }

    # Plot the subjects of non-selected ancestries to the back
    my $halfDot = $dotSize/2;
    foreach my $sbjNo (0 .. @sortSbjNos) {
	unless ($allSbjColorIds[$sbjNo]) {
	    my ($xVal, $yVal, $zVal, $junk) = @{$allPlotPts[$sbjNo]};

	    if ($xVal > $xMin && $xVal < $xMax && $yVal > $yMin && $yVal < $yMax) {
		my $color = $allSbjColors[$sbjNo];
		$color = $black unless ($hasRace);
		PlotOneDot($xVal, $yVal, $dotSize, $color);
	    }
	}
    }

    # Plot the subjects of selected ancestries to the front
    foreach my $sbjNo (0 .. @sortSbjNos) {
	if ($allSbjColorIds[$sbjNo]) {
	    my ($xVal, $yVal, $zVal, $junk) = @{$allPlotPts[$sbjNo]};

	    if ($xVal > $xMin && $xVal < $xMax && $yVal > $yMin && $yVal < $yMax) {
		my $color = $allSbjColors[$sbjNo];
		$color = $black unless ($hasRace);
		PlotOneDot($xVal, $yVal, $dotSize, $color);
	    }
	}
    }

    # Plot cutoff lines
    if ($showGd4) {
	PlotSasCutoffLines() if ($showCutoff);
    }
    else {
	my @vet0 = @{$allPlotPts[$numShowSbjs]};
	my @vet1 = @{$allPlotPts[$numShowSbjs+1]};
	my @vet2 = @{$allPlotPts[$numShowSbjs+2]};

	# Plot the triangle if all vertices are in the drawing area
	if ( $vet0[0] > $xMin && $vet0[0] < $xMax && $vet0[1] > $yMin && $vet0[1] < $yMax &&
	     $vet1[0] > $xMin && $vet1[0] < $xMax && $vet1[1] > $yMin && $vet1[1] < $yMax &&
	     $vet2[0] > $xMin && $vet2[0] < $xMax && $vet2[1] > $yMin && $vet2[1] < $yMax ) {
	    my ($x1, $y1) = (0, 0);
	    for my $i (0 .. 3) {
		my $popNo = $i % 3;
		my ($xVal, $yVal, $zVal, $junk) = @{$allPlotPts[$numShowSbjs+$popNo]};
		my $x2 = $gxLeft   + ($xVal - $xMin) * $graphWidth  / $xRange;
		my $y2 = $gyBottom - ($yVal - $yMin) * $graphHeight / $yRange;
		if ($i > 0) {
		    $im->line($x1, $y1, $x2, $y2, $black);
		}
		$x1 = $x2;
		$y1 = $y2;
	    }

	    my ($eLbl, $fLbl, $aLbl) = ("European", "African", "East Asian");
	    my $vetGap = 0.02;
	    my ($eLbly, $fLbly, $aLbly) = (1.52, 1.12, 1.12);
	    $eLbly = $eMaxy + $vetGap if ($eLbly < $eMaxy + $vetGap);
	    $fLbly = $fMiny - $vetGap if ($fLbly > $fMiny - $vetGap);
	    $aLbly = $aMiny - $vetGap if ($aLbly > $aMiny - $vetGap);

	    if (!$rotx && !$showCutoff) {
		PlotVertexLabel($eLbl, 2, $eLbly);
		PlotVertexLabel($fLbl, 3, $fLbly);
		PlotVertexLabel($aLbl, 1, $aLbly);
	    }

	    # Plot the cutoff lines
	    if ($showCutoff && !$rotx) {
		PlotCutoffLines();
	    }
	}

	my %whitePa = (xm => 0.9891, ym => 1.4703, ea99 => 0.0343, eb99 => 0.0391, ea95 => 0.0247, eb95 => 0.0282, ang => -1.6577);
	my %blackPa = (xm => 1.3016, ym => 1.2504, ea99 => 0.0324, eb99 => 0.1954, ea95 => 0.0210, eb95 => 0.1270, ang => -0.9152);
	my %ghanaPa = (xm => 1.2997, ym => 1.1671, ea99 => 0.0174, eb99 => 0.0263, ea95 => 0.0122, eb95 => 0.0185, ang => -1.6387);
	my %hisp1Pa = (xm => 1.2935, ym => 1.4873, ea99 => 0.0427, eb99 => 0.2384, ea95 => 0.0305, eb95 => 0.1705, ang => -0.9952);
	my %asianPa = (xm => 1.0451, ym => 1.6954, ea99 => 0.0278, eb99 => 0.0997, ea95 => 0.0165, eb95 => 0.0591, ang => -2.0103);
	my %chinaPa = (xm => 1.1253, ym => 1.7390, ea99 => 0.0163, eb99 => 0.0240, ea95 => 0.0132, eb95 => 0.0194, ang => -1.8665);
	my %indiaPa = (xm => 0.9053, ym => 1.3183, ea99 => 0.0262, eb99 => 0.0662, ea95 => 0.0196, eb95 => 0.0496, ang => -2.2960);
	my %mexicPa = (xm => 0.8988, ym => 1.3423, ea99 => 0.0402, eb99 => 0.1625, ea95 => 0.0291, eb95 => 0.1175, ang => -2.3087);
	my $bx0 = $bCtr[0];
	my $by0 = $bCtr[1];

	if (!$rotx && $whitePs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($whitePa{xm}, $whitePa{ym}, $whitePa{ea95}, $whitePa{eb95}, $whitePa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $blackPs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($blackPa{xm}, $blackPa{ym}, $blackPa{ea95}, $blackPa{eb95}, $blackPa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $asianPs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($asianPa{xm}, $asianPa{ym}, $asianPa{ea95}, $asianPa{eb95}, $asianPa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $ghanaPs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($ghanaPa{xm}, $ghanaPa{ym}, $ghanaPa{ea95}, $ghanaPa{eb95}, $ghanaPa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $chinaPs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($chinaPa{xm}, $chinaPa{ym}, $chinaPa{ea95}, $chinaPa{eb95}, $chinaPa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $indiaPs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($indiaPa{xm}, $indiaPa{ym}, $indiaPa{ea95}, $indiaPa{eb95}, $indiaPa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $mexicPs) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($mexicPa{xm}, $mexicPa{ym}, $mexicPa{ea95}, $mexicPa{eb95}, $mexicPa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
	if (!$rotx && $hisp1Ps) {
	    my ($meanx, $meany, $ea9500, $eb9500, $backAng) = ($hisp1Pa{xm}, $hisp1Pa{ym}, $hisp1Pa{ea95}, $hisp1Pa{eb95}, $hisp1Pa{ang});
	    PlotEllipse ($meanx, $meany, $ea9500, $eb9500, $bx0, $by0, $backAng, $expAreaColor);
	}
    }
}

#
# Move and rotate a 2D point (given x, y coordinates
#
sub TransformPoint2D
{
    my ($x1, $y1, $dx, $dy, $alpha) = @_;

    my $xt = $x1 + $dx;
    my $yt = $y1 + $dy;

    my $x = $xt * cos($alpha) - $yt * sin($alpha);
    my $y = $yt * cos($alpha) + $xt * sin($alpha);

    $x = $x - $dx;
    $y = $y - $dy;

    return ($x, $y);
}

#
# Rotate a 3D point (array reference) around an axis
#
sub RotatePoint3D
{
    my ($pRef, $deg, $axis) = @_; # deg = angle in degree; axis = x, y or z

    my $rad = $deg * $pi / 180;
    my $sind = sin($rad);
    my $cosd = cos($rad);

    my $t1 = [];
    if ($axis =~ /x/) {
	$t1 = [
	      [1, 0,      0,     0],
	      [0, $cosd, -$sind, 0],
	      [0, $sind,  $cosd, 0],
	      [0, 0,      0,     1]  ];
    }
    elsif ($axis =~ /y/) {
	$t1 = [
	       [ $cosd, 0,  $sind, 0],
	       [ 0,     1,  0,     0],
	       [-$sind, 0,  $cosd, 0],
	       [ 0,     0,  0,     1]  ];
    }
    elsif ($axis =~ /z/) {
	$t1 = [
	       [$cosd, -$sind, 0, 0],
	       [$sind,  $cosd, 0, 0],
	       [0,      0,     1, 0],
	       [0,      0,     0, 1]  ];
    }

    my @p2 = TransformPoint3D($pRef, $t1);

    return @p2;
}

#
# Move a 3D point
#
sub MovePoint3D
{
    my ($pRef, $dx, $dy, $dz) = @_;

    my $t1 = [ [1, 0, 0, $dx],
	       [0, 1, 0, $dy],
	       [0, 0, 1, $dz],
	       [0, 0, 0, 1  ]  ];

    my @p2 = TransformPoint3D($pRef, $t1);

    return @p2;
}

#
# Transform a 3D point given a transformation matrix
#
sub TransformPoint3D
{
    my ($pRef, $tRef) = @_;

    my @p = @$pRef;
    my @t = @$tRef;

    my @q = (0, 0, 0, 0);

    for my $i (0 .. 3) {
	for my $j (0 .. 3) {
	    $q[$i] += $t[$i]->[$j] * $p[$j];
	}
    }

    return @q;
}

#
# Rotate all points in the global array around an axis
#
sub RotateShape3D
{
    my ($deg, $axis, $test) = @_;

    my $rad = $deg * $pi / 180;
    my $sind = sin($rad);
    my $cosd = cos($rad);

    my $t1 = [];
    if ($axis =~ /x/) {
	$t1 = [
	      [1, 0,      0,     0],
	      [0, $cosd, -$sind, 0],
	      [0, $sind,  $cosd, 0],
	      [0, 0,      0,     1]  ];
    }
    elsif ($axis =~ /y/) {
	$t1 = [
	       [ $cosd, 0,  $sind, 0],
	       [ 0,     1,  0,     0],
	       [-$sind, 0,  $cosd, 0],
	       [ 0,     0,  0,     1]  ];
    }
    elsif ($axis =~ /z/) {
	$t1 = [
	       [$cosd, -$sind, 0, 0],
	       [$sind,  $cosd, 0, 0],
	       [0,      0,     1, 0],
	       [0,      0,     0, 1]  ];
    }

    TransformShape3D($t1, $test);
}

#
# Scale all points in the global array
#
sub ScaleShape3D
{
    my ($sx, $sy, $sz) = @_;
    $sy = $sx unless ($sy);
    $sz = $sx unless ($sz);

    my $t1 = [ [$sx, 0,   0,   0],
	       [0,   $sy, 0,   0],
	       [0,   0,   $sz, 0],
	       [0,   0,   0,   1]  ];

    TransformShape3D($t1);
}

#
# Move all points in the global array
#
sub MoveShape3D
{
    my ($dx, $dy, $dz) = @_;

    my $t1 = [ [1, 0, 0, $dx],
	       [0, 1, 0, $dy],
	       [0, 0, 1, $dz],
	       [0, 0, 0, 1  ]  ];

    TransformShape3D($t1);
}


#
# Move and rotate all points in the global array to project them on the EFA plane
#
sub TransformShape3D
{
    my ($tRef) = @_;

     for my $i (0 .. $#allPlotPts) {
	my $msgy = 0;
	my @p = @{$allPlotPts[$i]};
	if ($i > $#allPlotPts-2) {
	    my @q = TransformPoint3D(\@p, $tRef);
	    $allPlotPts[$i] = \@q;
	}
	else {
	    my @q = TransformPoint3D(\@p, $tRef);
	    $allPlotPts[$i] = \@q;
	}
    }
}

#
# Plot the x and y axes, and the tick marks
#
sub PlotAxes
{
    my ($yTop, $xLabelStr, $yLabelStr) = @_;

    my $yBottom = $yTop + $graphHeight;

    # Plot y-axis
    $im->line($gxLeft, $yTop, $gxLeft, $yBottom, $black);
    $im->line($gxRight, $yTop, $gxRight, $yBottom, $black);

    # Plot ticks
    my $yStep = 0.1;
    my $yTicks = $yRange / $yStep;
    my $yMinorTicks = 5;

    my $xMajor = $gxLeft - $majorTickLen;
    my $dy = $graphHeight/$yTicks;
    my $minorDy = $graphHeight/($yTicks*$yMinorTicks);

    my $tickLblLen = 4;
    my $numDigits = $tickLblLen;;
    my $format = "%4.1f";
    my $xString = $gxLeft - $mbFontWidth*$tickLblLen - $majorTickLen - $tickGap;

    for my $i (0 .. $yTicks+1) {
	my $val = $yMin + $yStep * $i;
	my $valStr = sprintf($format, $val);
	my $y = int($yBottom - $dy*$i + 0.5);
	$y = $yBottom if ($y > $yBottom);
	if ($y >= $yTop) {
	    $im->line($gxLeft, $y, $xMajor, $y, $black);
	    $im->string(gdMediumBoldFont, $xString, $y-$mbFontHeight/2, $valStr, $black);
	}

	# Plot minor ticks
	if ($i < $yTicks) {
	    my $yMinor = $y;
	    my $xMinor = $gxLeft - $minorTickLen;
	    for my $j (1 .. $yMinorTicks-1) {
		$yMinor -= $minorDy;
		$im->line($gxLeft, $yMinor, $xMinor, $yMinor, $black) if ($yMinor > $yTop);
	    }
	}
    }

    # Draw y-axis label
    my $lblLen = length($yLabelStr);
    my $yLabel = $yBottom - ($graphHeight - $mbFontWidth*$lblLen)/2;
    $im->stringUp(gdMediumBoldFont, $xString-20, $yLabel, $yLabelStr, $black);


    # Plot x-axis
    $im->line($gxLeft, $yBottom, $gxRight, $yBottom, $black);
    $im->line($gxLeft, $yTop, $gxRight, $yTop, $black);

    # Plot ticks
    my $xStep = 0.1;
    my $xTicks = int($xRange / $xStep + 0.05);
    my $xMinorTicks = 5;

    my $yMajor = $yBottom + $majorTickLen;
    my $dx = $graphWidth/$xTicks;
    my $minorDx = $graphWidth/($xTicks*$xMinorTicks);

    my $yString = $yMajor + 5;

    for my $i (0 .. $xTicks) {
	my $val = $xMin + $xStep * $i;
	my $valStr = sprintf($format, $val);
	my $x = int($gxLeft + $dx*$i + 0.5);
	$im->line($x, $yBottom, $x, $yMajor, $black);
	$im->string(gdMediumBoldFont, $x-$mbFontWidth*2, $yString, $valStr, $black);

	# Plot minor ticks
	if ($i < $xTicks) {
	    my $xMinor = $x;
	    my $yMinor = $yBottom + $minorTickLen;
	    for my $j (1 .. $xMinorTicks-1) {
		$xMinor += $minorDx;
		$im->line($xMinor, $yBottom, $xMinor, $yMinor, $black);
	    }
	}
    }

    # Draw x-axis label
    my $lblLen = length($xLabelStr);
    my $xLabel = $gxLeft + ($graphWidth - $mbFontWidth*$lblLen)/2;
    $im->string(gdMediumBoldFont, $xLabel, $yMajor + 25, $xLabelStr, $black);
}

#
# Save scores, colors, ancestries, etc, to arrays
#
sub GetSubjectPositions
{
    foreach my $sbjNo (keys %sbjPopScores) {
	my %info = %{$sbjPopScores{$sbjNo}};
	my $sbj  = $sbjPopScores{$sbjNo}->{subject};
	my $snps = $sbjPopScores{$sbjNo}->{snps};
	my $race = $unkRace;
	$race = $sbjRaces{$sbj} if ($sbjRaces{$sbj});

	my $indPak  = $info{indPak}  ? $info{indPak}  : 100;
	my $mexican = $info{mexican} ? $info{mexican} : 0;

	my $gd1 = $info{x};
	my $gd2 = $info{y};
	my $gd3 = $info{z};
	my $gd4 = $mexican - $indPak;

	my @pvals = ($gd1, $gd2, $gd3, 1);
	if ($showGd4) {
	    @pvals = ($gd1, $gd4, $gd3, 1);
	}
	push @allSbjPvalues, \@pvals;

	my $raceId = $selRaceNames{$race};
	my $raceNo = $raceIdNos{$raceId};

	my $colorNo = $raceColorIds{$race};
	my $color = $raceColors[$colorNo];

	push @allSbjColors, $color;
	push @allSbjColorIds, $colorNo;
	push @allSbjIds, $sbj;
	push @allSbjRaces, $race;
	push @allSbjSnps, $snps;
	push @sbjIndPaks, $indPak;
	push @sbjMexicans, $mexican;
    }
}

#
# Show time passed
#
sub ShowTimeDifference
{
    my ($t1, $t2) = @_;
	
    my $t1sec = $$t1[0];
    my $t1us  = $$t1[1];
    my $t2sec = $$t2[0];
    my $t2us  = $$t2[1];
	
    my $ds = $t2sec - $t1sec;
    my $dus = $t2us - $t1us;
	
    if ($dus < 0) {
	$dus += 1000000;
	$ds -= 1;
    }
		
    if ($ds < 0 || ($ds == 0 && $dus < 0)) {
	print "Error: T1 > T2\n";
	return;
    }
	
    print "Time used: ";

    my $dds = $ds % 60;
    $dds += $dus/1000000;

    if ($ds == 0) {
	print "$dus micro seconds\n";
	return;
    }

    my $dh = int($ds/3600);
    $ds = $ds % 3600;
    my $dm = int($ds/60);

    print "$dh hour" if ($dh > 0);
    print "s" if ($dh > 1);
    print " $dm minute" if ($dh > 0 || $dm > 0);
    print "s" if ($dm > 1);
    print " $dds seconds\n";
}

#
# Read parameters from command line
#
sub GetParameters
{
    my @args = @_;

    my ($inFile, $outFile) = ("", "");
    my %params = ();

    if (@args > 1) {
	$inFile = $args[0];
	$outFile = $args[1];

	for my $i (2 .. $#args) {
	    my $arg = $args[$i];
	    if ($arg =~ /^-([A-Za-z]\w*)$/) {
		my $param = $1;
		if ($i == $#args) {
		    $params{$param} = 1;
		}
		else {
		    my $val = $args[$i+1];
		    if ($val =~ /^-[A-Za-z]\w*$/) {
			$params{$param} = 1;
		    }
		    else {
			$params{$param} = $val;
			$i++;
		    }
		}
	    }
	}
    }

    if ($inFile =~ /^-/ || $outFile =~ /^-/) {
	($inFile, $outFile) = ("", "");
    }

    return ($inFile, $outFile, \%params);
}

#
# Assign population based on the calculated ancestry protortions and GD4
#
sub GetGenoPopId
{
    my ($gd1, $gd4, $ePct, $fPct, $aPct) = @_;

    my $genoPopId = 0;

    my $isSas = 0;
    my $y = $sasCutBasey + $sasCutaVal * ($gd1 - $meanSasx)**2;
    $isSas = 1 if ($gd4 > $y);

    my $isAsn = 0;
    my $x = $asianCutBasex + $asianCutaVal * ($gd4 - $meanAsiany)**2;
    $isAsn = 1 if ($gd1 > $x);

    if ($ePct > $whiteCut) {
	$genoPopId = 1;
    }
    elsif ($fPct > $blackCut) {
	$genoPopId = 2;
    }
    elsif ($aPct > $asianCut) {
	$genoPopId = 3;
    }
    elsif ($isSas) {
	$genoPopId = 8;
    }
    elsif ($fPct < $oHispCut) {
	if ($isAsn) {
	    $genoPopId = 7;
	}
	else {
	    if ($aPct < $fPct) {
		$genoPopId = 5;
	    }
	    else {
		if ($gd4 < 0) {
		    $genoPopId = 6; # Hispanic2 has gd4 < 0
		}
		else {
		    # Set pop ID to other. These are not Hispanics,
		    # probably European/Asian mistures. Set to Other.
		    $genoPopId = 9;
		}
	    }
	}
    }
    elsif ($aPct < $oHispCut) {
	if ($fPct > $bHispCut) {
	    $genoPopId = 4;
	}
	else {
	    $genoPopId = 5;
	}
    }
    else {
	$genoPopId = 9;
    }

    return $genoPopId;
}
