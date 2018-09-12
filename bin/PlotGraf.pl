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
        File Description: script to plot graphs to show distribution of genotype mismatch rates
        Date: 08/15/2016
    ====================================================================================
EOF

use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use strict;
use Data::Dumper;
use GD;
use GD::Text;
use GD::Graph;
use GD::Graph::colour;
use GD::Graph::lines;

#------------------------------------------- Take and set parameters ----------------------------------------#
my $pi = 3.1415926;

# Mean HGMR and cutoff values
my $fsHgmr = 5.1;
my $d2Hgmr = 11.7;
my $d3Hgmr = 17.9;
my $poCutMr = 1.4;
my $fsCutMr = 8.0;
my $d2CutMr = 15.0;
my $d3CutMr = 22.0;

my $idCutAg = 0.16;
my $fsCutAg = 0.36;
my $poCutAg = 0.44;
my $d2CutAg = 0.485;
my $d3CutAg = 0.52;

# Parameters for calculating probability densities of bivariate normal distributions
my $poHgmrK = 10.0;
my ($poMux, $poSx, $poMuy, $poSy, $poCor) = ( 0.037292, 0.073727, 39.477589, 1.468133, 0.104462);
my ($fsMux, $fsSx, $fsMuy, $fsSy, $fsCor) = ( 4.863680, 1.015553, 33.592671, 2.322399, 0.784207);
my ($d2Mux, $d2Sx, $d2Muy, $d2Sy, $d2Cor) = (11.162160, 1.353599, 47.304737, 1.096585, 0.802548);
my ($d3Mux, $d3Sx, $d3Muy, $d3Sy, $d3Cor) = (17.320116, 1.354330, 51.069507, 0.958527, 0.830226);

my ($inFile, $outFile, $graphType, $params) = GetParameters(@ARGV);

unless ($inFile && $outFile && $graphType) {
    print "\n$disclaim\n";
    print "\nUsage: PlotGraf.pl <input related subject file> <output png file> <graph type> [Options]\n\n";
    print "Note:\n";
    print "    Valid graph types are:\n";
    print "        1 = HGMR histogram\n";
    print "        2 = AGMR histogram\n";
    print "        3 = HGMR + AGMR scatter plot\n\n";
    print "Options:\n";
    print "    -gw     graph width:  Set graph width in pixels\n";
    print "    -gh     graph height: Set graph height in pixels\n";
    print "    -xmax   max x value:  Set maximum HGMR or AGMR on x-axis of the histogram\n";
    print "    -ymax   max y value:  Set maximum number of pairs on y-axis of the histogram\n";
    print "    -dot    size:         Set dot size in pixels on the scatter plot\n";
    print "    -hfd    size:         Set dot size in pixels for HF (half sib + full cousin) pairs\n";
    print "\n";
    exit 0;
}

if (-e $outFile) {
    print "\nERROR: output file $outFile already exists!\n\n";
    exit;
}

my $binSize   = $params->{"bin"} ;
my $maxPair   = $params->{"ymax"};
my $maxMr     = $params->{"xmax"};
my $xStep     = $params->{"xtick"};
my $yStep     = $params->{"ytick"};
my $dotSize   = $params->{"dot"} ? $params->{"dot"} : 2;
my $hfdSize   = $params->{"hfd"} ? $params->{"hfd"} : 2;
my $fillDot = 1;
if ($dotSize < 0) {
    $dotSize *= -1;
    $fillDot = 0;
}
$fillDot = 0 if ($dotSize < 3);

$dotSize = 8 if ($dotSize > 8);
$dotSize = 1 if ($dotSize < 1);
$hfdSize = 8 if ($hfdSize > 8);
$hfdSize = 1 if ($hfdSize < 1);

$maxPair = 0 if ($maxPair !~ /^\d+$/);
$maxMr = 100 if ($maxMr > 100);

if ($graphType == 2) {
    $maxMr   = 100 unless ($maxMr);
    $binSize = 1   unless ($binSize);
    $xStep   = 5   unless ($xStep);
    $binSize = 1   unless ($binSize);
}
else {
    $maxMr = 25 unless ($maxMr);
    $xStep = 1 unless ($xStep);
    $binSize = 0.2 unless ($binSize);
}

my $xMinorTicks = $params->{"minorTicks"};
$xMinorTicks = 5 unless ($xMinorTicks);


### Set window sizes and gaps

my $leftEdge = 80;      # The left edge of the graph
my $graphWidth = 1200;  # Width of the plot (rectangle box only)
my $rightEdge = 20;     # The right edge

my $topEdge = 10;       # The top edge
my $graphHeight = 600;  # Height of the plot
my $bottomEdge = 180;   # The bottom edge

my $axisLabelGap = 20;  # Distance between tick labels and axis labels
my $tickGap = 5;        # Gap between the tick and tick label
my $graphGap = 100;

my $majorTickLen = 6;  # major tick length (in pixels)
my $minorTickLen = 3;  # minor tick length
my $graphTitleHt = 20;  # Height of the graph title

$graphWidth  = $params->{"gw"} if ($params->{"gw"});
$graphHeight = $params->{"gh"} if ($params->{"gh"});
$graphWidth  = 600 if ($graphWidth  < 600);
$graphHeight = 200 if ($graphHeight < 200);

# Font sizes
my $mbFontWidth  = gdMediumBoldFont->width;
my $mbFontHeight = gdMediumBoldFont->height;
my $mlFontWidth  = gdLargeFont->width;
my $mlFontHeight = gdLargeFont->height;

# Width and height of the whole graph
my $imageWidth   = $leftEdge + $graphWidth + $rightEdge;
my $imageHeight  = $topEdge + $mlFontHeight*3 + $graphHeight + $bottomEdge;

my $mlFontWidth  = gdLargeFont->width;
my $mlFontHeight = gdLargeFont->height;

my $mgFontWidth  = gdGiantFont->width;
my $mgFontHeight = gdGiantFont->height;


#------------------------------------------- Set colors  ----------------------------------------------#
my $cgi = new CGI;
my $im = new GD::Image($imageWidth, $imageHeight);

my $white  = $im->colorAllocate(255,255,255);
my $grey   = $im->colorAllocate(180,180,180);
my $red    = $im->colorAllocate(250,0,0);
my $brown  = $im->colorAllocate(220,120,100);
my $cyan   = $im->colorAllocate(0,205,205);
my $green  = $im->colorAllocate(0,155,0);
my $blue   = $im->colorAllocate(0,0,255);
my $black  = $im->colorAllocate(0,  0,  0);
my $yellow = $im->colorAllocate(255,  255,  0);
my $deepYl = $im->colorAllocate(165,  155,  0);
my $gold   = $im->colorAllocate(204, 153,  80);
my $purple = $im->colorAllocate(153, 0,  153);
my $caseC  = $im->colorAllocate(51, 102, 153);
my $dGrey  = $im->colorAllocate(100,100,100);

# Colors of relationship types on the HGMR histogram
my $poColor = $red;
my $fsColor = $blue;
my $d2Color = $green;
my $d3Color = $yellow;
my $hcColor = $white;
my $hfColor = $purple;
my $unColor = $grey;

# Colors of relationship types on the AGMR histogram
my $dpColor = $brown;
my $mtColor = $purple;
my $orColor = $blue;

# Colors of relationship types on the scatter plot
my $d3ScatColor = $deepYl;
my $unGreyVal = 50;
my $unScatColor = $im->colorAllocate($unGreyVal, $unGreyVal, $unGreyVal);

# Colors of contour lines
my $poContour = $im->colorAllocate(150, 150,   0);
my $fsContour = $im->colorAllocate(220, 120, 100);
my $d2Contour = $im->colorAllocate(200, 100, 200);
my $d3Contour = $im->colorAllocate(204, 153,  80);


#----------------------------------------- Read related pairs and plot the grah ------------------------#

my ($dupRef, $relRef, $error) = ReadRelatedPairs($inFile);
if ($error) {
    print "$error\n";
    exit;
}

$fsCutMr = ($fsMux + $d2Mux) / 2;;
$d2CutMr = ($d2Mux + $d3Mux) / 2;;
$d3CutMr = $d3Mux + ($d3Mux - $d2Mux) / 2;

$idCutAg = 23.0;;
$fsCutAg = ($fsMuy + $poMuy) / 2;
$poCutAg = ($poMuy + $d2Muy) / 2;
$d2CutAg = ($d2Muy + $d3Muy) / 2;
$d3CutAg = $d3Muy + ($d3Muy - $d2Muy) / 2;

print "Expected HGMR and AGMR values for different types of relationships:\n";
printf("\tFS HGMR %7.4f\%\n", $fsMux);
printf("\tD2 HGMR %7.4f\%\n", $d2Mux);
printf("\tD3 HGMR %7.4f\%\n", $d3Mux);
printf("\tPO AGMR %7.4f\%\n", $poMuy);
printf("\tFS AGMR %7.4f\%\n", $fsMuy);
printf("\tD2 AGMR %7.4f\%\n", $d2Muy);
printf("\tD3 AGMR %7.4f\%\n", $d3Muy);
print "\n";

my %dupPairs = %$dupRef;
my %relPairs = %$relRef;
my $numDups = keys %dupPairs;
my $numRels = keys %relPairs;
print "Found $numDups pairs of dups, $numRels pairs of relatives.\n";

my $numBins = int($maxMr / $binSize);
my @binPairs = ();
my @mzPairs = ();
my @duPairs = ();
my @poPairs = ();
my @fsPairs = ();
my @hfPairs = ();
my @d2Pairs = ();
my @d3Pairs = ();
my @hcPairs = ();

my @binDupPairs = ();
my @mtPairs = ();
my @dpPairs = ();
my @p0Pairs = (); # PO from dup list
my @p1Pairs = (); # FS
my @p2Pairs = (); # D2
my @p3Pairs = (); # D3

my $hasHF = 0;
foreach my $pair (keys %relPairs) {
    my $pRel = $relPairs{$pair}->{pedRel};
    if ($pRel && $pRel eq "HF") {
	$hasHF = 1;
	last;
    }
}

my $hasPlot = 0;

my $totRelPairs = keys %relPairs;
my $totDupPairs = keys %dupPairs;

if ($totRelPairs + $totDupPairs < 1) {
    print "No related subjects found.\n\n";
    exit;
}

if ($graphType == 1) {  # Plot histogram to show HGMR distribution
    if ($totRelPairs > 0) {
	for my $binNo (0 .. $numBins) {
	    push @binPairs, 0;
	    push @mzPairs, 0;
	    push @duPairs, 0;
	    push @poPairs, 0;
	    push @fsPairs, 0;
	    push @d2Pairs, 0;
	    push @d3Pairs, 0;
	}

	foreach my $pair (keys %relPairs) {
	    my %info = %{$relPairs{$pair}};
	    my $smp1 = $info{smp1};
	    my $smp2 = $info{smp2};
	    my $hgmr = $info{hgmr};
	    my $pRel = $info{pedRel};

	    if ($hgmr <= $maxMr) {
		my $binNo = int ($hgmr/$binSize);
		$binPairs[$binNo]++;
		$mzPairs[$binNo]++ if ($pRel eq "MT");
		$duPairs[$binNo]++ if ($pRel eq "DP");
		$poPairs[$binNo]++ if ($pRel eq "PO");
		$fsPairs[$binNo]++ if ($pRel eq "SB");  # Sibling with one parent miss is very rare, treated as FS
		$hfPairs[$binNo]++ if ($pRel eq "HF");
		$d2Pairs[$binNo]++ if ($pRel eq "HS" || $pRel eq "AV" || $pRel eq "AU" || $pRel eq "GP");
		$d3Pairs[$binNo]++ if ($pRel eq "CS" || $pRel eq "FC" || $pRel eq "HA");
		$hcPairs[$binNo]++ if ($pRel eq "HC");
	    }
	}

	PlotRelativeHistogram();
	$hasPlot = 1;
    }
    else {
	print "No output file was produced because you selected mode 1, which is for non-identical related subjects, but your input data contain 0 non-identical related subjects.\n\n";
    }
}
elsif ($graphType == 2) { # Plot histogram to show AGMR distribution
    if ($totDupPairs > 0) {
	for my $binNo (0 .. $numBins) {
	    push @binDupPairs, 0;
	    push @mtPairs, 0;
	    push @dpPairs, 0;
	    push @p0Pairs, 0;
	    push @p1Pairs, 0;
	    push @p2Pairs, 0;
	    push @p3Pairs, 0;
	}

	foreach my $pair (keys %dupPairs) {
	    my %info = %{$dupPairs{$pair}};
	    my $smp1 = $info{smp1};
	    my $smp2 = $info{smp2};
	    my $agmr = $info{agmr};
	    my $pRel = $info{pedRel};

	    if (($pRel && $pRel =~ /(MT)|(DP)/) || ($agmr <= $maxMr)) {
		my $binNo = int ($agmr/$binSize);
		$binDupPairs[$binNo]++;
		$mtPairs[$binNo]++ if ($pRel eq "MT");
		$dpPairs[$binNo]++ if ($pRel eq "DP");
		$p0Pairs[$binNo]++ if ($pRel eq "PO");
		$p1Pairs[$binNo]++ if ($pRel eq "FS");
		$p2Pairs[$binNo]++ if ($pRel eq "HS" || $pRel eq "AV" || $pRel eq "AU" || $pRel eq "GP");
		$p3Pairs[$binNo]++ if ($pRel eq "CS" || $pRel eq "FC" || $pRel eq "HA");
	    }
	}

	PlotDuplicateHistogram();
	$hasPlot = 1;
    }
    else {
	print "No output file was produced because you selected mode 2, which is for dups, but your input data contain 0 dups.\n\n";
    }
}
elsif ($graphType == 3) {  # Scatter plot to show distribution of both HGMR and AGMR
    if ($totRelPairs > 0) {
	PlotScatterPlot();
	$hasPlot = 1;
    }
    else {
	print "No output file was produced because you selected mode 3, which is for non-identical related subjects, but your input data contain 0 non-identical related subjects.\n\n";
    }
}

if ($hasPlot) {
    print $cgi->header( {-type=>'image/png'});

    open(IMG, ">$outFile") or die $!;
    binmode IMG;

    print IMG $im->png;

    print "Graph saved to $outFile\n\n";
    exit;
}


#
# Histogram to show HGMR distribution
#
sub PlotRelativeHistogram
{
    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yTop    = $topEdge + $mlFontWidth * 3;
    my $yBottom = $yTop + $graphHeight;
    my $color   = $black;

    my $relLblGap = 10;
    my $graphTitle = "Homozygous Genotype Mismatch Rates of Closely Related Subjects";
    $yTop = PlotGraphTitle($graphTitle, $xLeft, $yTop);
    $yTop += $mbFontHeight + $relLblGap;     # Leave space for relationship type labels on top of the graph
    $yBottom = $yTop + $graphHeight;

    my $xSt  = 0;
    my $xMax = $maxMr;

    # Find maximum number of pairs in each bar
    my $yMax = 0;
    my $binNo = 0;
    my $binX = $binNo * $binSize;
    while ($binX < $xMax) {
	my $cnt = $binPairs[$binNo];
	$yMax = $cnt if ($cnt > $yMax);
	$binNo++;
	$binX = $binNo * $binSize;
    }

    $yMax = $maxPair if ($maxPair);

    # Plot axes
    ($xMax, $yMax) = PlotHistoAxes($yTop, $yMax, $xSt, $xMax, 1);

    # Plot the histogram
    my $binWid = $binSize * $graphWidth / $xMax;
    $binNo = 0;
    $binX = $binNo * $binSize;
    my $prevXpos = 0;
    while ($binX < $xMax) {
	my $cnt = $binPairs[$binNo];

	my $mzCnt = $mzPairs[$binNo];
	my $duCnt = $duPairs[$binNo];
	my $poCnt = $poPairs[$binNo];
	my $fsCnt = $fsPairs[$binNo];
	my $d2Cnt = $d2Pairs[$binNo];
	my $d3Cnt = $d3Pairs[$binNo];
	my $hcCnt = $hcPairs[$binNo];  # Half cousins, won't show them on histogram

	my $xPos = $xLeft + $binX * $graphWidth / $xMax + $binWid;
	my $yPos = $yBottom - ($cnt-$hcCnt) * $graphHeight / $yMax;
	$yPos = $yTop if ($yPos < $yTop);

	my $mzPos = $yBottom - $mzCnt * $graphHeight / $yMax;
	my $duPos = $mzPos - $duCnt * $graphHeight / $yMax;
	my $poPos = $duPos - $poCnt * $graphHeight / $yMax;
	my $fsPos = $poPos - $fsCnt * $graphHeight / $yMax;
	my $d2Pos = $fsPos - $d2Cnt * $graphHeight / $yMax;
	my $d3Pos = $d2Pos - $d3Cnt * $graphHeight / $yMax;

	$mzPos = $yTop if ($mzPos < $yTop);
	$duPos = $yTop if ($duPos < $yTop);
	$poPos = $yTop if ($poPos < $yTop);
	$fsPos = $yTop if ($fsPos < $yTop);
	$d2Pos = $yTop if ($d2Pos < $yTop);
	$d3Pos = $yTop if ($d3Pos < $yTop);

	if ($yBottom - $yPos > 0.5) {
	    $prevXpos = $xPos - $binWid unless ($prevXpos);
	    $im->filledRectangle($prevXpos, $yBottom, $xPos, $mzPos, $mtColor) if ($yBottom - $mzPos > 0.5);
	    $im->filledRectangle($prevXpos, $mzPos, $xPos, $duPos, $dpColor) if ($mzPos - $duPos > 0.5);
	    $im->filledRectangle($prevXpos, $duPos, $xPos, $poPos, $poColor) if ($duPos - $poPos > 0.5);
	    $im->filledRectangle($prevXpos, $poPos, $xPos, $fsPos, $fsColor) if ($poPos - $fsPos > 0.5);
	    $im->filledRectangle($prevXpos, $fsPos, $xPos, $d2Pos, $d2Color) if ($fsPos - $d2Pos > 0.5);
	    $im->filledRectangle($prevXpos, $d2Pos, $xPos, $d3Pos, $d3Color) if ($d2Pos - $d3Pos > 0.5);
	    $im->filledRectangle($prevXpos, $d3Pos, $xPos, $yPos,  $unColor) if ($d3Pos - $yPos > 0.5);

	    $im->rectangle($prevXpos, $yBottom, $xPos, $yPos, $color);
	}

	$binNo++;
	$binX = $binNo * $binSize;
	$prevXpos = $xPos;
    }

    # Lable the regions of different relationships on top of the graph
    my $x1 = $xLeft + $poCutMr * $graphWidth / $xMax;
    my $x2 = $xLeft + $fsCutMr * $graphWidth / $xMax;
    my $x3 = $xLeft + $d2CutMr * $graphWidth / $xMax;
    my $x4 = $xLeft + $d3CutMr * $graphWidth / $xMax;
    $im->line($x1, $yTop, $x1, $yBottom, $cyan);
    $im->line($x2, $yTop, $x2, $yBottom, $cyan);
    $im->line($x3, $yTop, $x3, $yBottom, $cyan);
    $im->line($x4, $yTop, $x4, $yBottom, $cyan);
    my $poHead = "PO";
    my $fsHead = "FS";
    my $d2Head = "D2";
    my $d3Head = "D3";
    my $unHead = "UN";

    my $poHeadLen = length($poHead) * $mbFontWidth;
    my $fsHeadLen = length($fsHead) * $mbFontWidth;
    my $d2HeadLen = length($d2Head) * $mbFontWidth;
    my $d3HeadLen = length($d3Head) * $mbFontWidth;
    my $unHeadLen = length($unHead) * $mbFontWidth;

    my $poHeadX = ($xLeft + $x1) / 2.0 - $poHeadLen / 2.0;
    my $fsHeadX = ($x1 + $x2) / 2.0 - $fsHeadLen / 2.0;
    my $d2HeadX = ($x2 + $x3) / 2.0 - $d2HeadLen / 2.0;
    my $d3HeadX = ($x3 + $x4) / 2.0 - $d3HeadLen / 2.0;
    my $unHeadX = ($x4 + $xLeft + $graphWidth) / 2.0 - $unHeadLen / 2.0;

    my $yHeadPos = $yTop - 20;
    $im->string(gdMediumBoldFont, $poHeadX, $yHeadPos, $poHead, $black);
    $im->string(gdMediumBoldFont, $fsHeadX, $yHeadPos, $fsHead, $black);
    $im->string(gdMediumBoldFont, $d2HeadX, $yHeadPos, $d2Head, $black);
    $im->string(gdMediumBoldFont, $d3HeadX, $yHeadPos, $d3Head, $black);
    $im->string(gdMediumBoldFont, $unHeadX, $yHeadPos, $unHead, $black);

    ShowBarLegends($xLeft, $yBottom);
}


#
# Histogram to show AGMR distribution
#
sub PlotDuplicateHistogram
{
    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yTop    = $topEdge + $mlFontWidth * 3;
    my $yBottom = $yTop + $graphHeight;
    my $color   = $black;

    my $relLblGap = 0;
    my $graphTitle = "Distribution of AGMR of Duplicate Samples and Monozygous Twins";
    $yTop = PlotGraphTitle($graphTitle, $xLeft, $yTop);
    $yTop += $mbFontHeight + $relLblGap;
    $yBottom = $yTop + $graphHeight;

    my $xSt = 0;
    my $xMax = $maxMr;

    # Find maximum number of pairs in each bar
    my $yMax = 0;
    my $binNo = 0;
    my $binX = $binNo * $binSize;
    while ($binX < $xMax) {
	my $cnt = $binDupPairs[$binNo];
	$yMax = $cnt if ($cnt > $yMax);
	$binNo++;
	$binX = $binNo * $binSize;
    }
    $yMax = $maxPair if ($maxPair);

    # Plot axes
    ($xMax, $yMax) = PlotHistoAxes($yTop, $yMax, $xSt, $xMax);

    # Plot cutoff lines
    my $x1 = $xLeft + $idCutAg * $graphWidth / $xMax;
    my $x2 = $xLeft + $fsCutAg * $graphWidth / $xMax;
    my $x3 = $xLeft + $poCutAg * $graphWidth / $xMax;
    my $x4 = $xLeft + $d3CutAg * $graphWidth / $xMax;
    $im->line($x1, $yTop, $x1, $yBottom, $cyan);
    $im->line($x2, $yTop, $x2, $yBottom, $cyan);
    $im->line($x3, $yTop, $x3, $yBottom, $cyan);
    $im->line($x4, $yTop, $x4, $yBottom, $cyan);

    my $idHead = "Duplicates or MZ twins";
    my $fsHead = "Full sibling";
    my $poHead = "Parent-offspring";
    my $d2Head = "2nd or 3rd degree";
    my $unHead = "Unrelated";

    if ($graphWidth < 80000) {
	$idHead = "ID";
	$poHead = "PO";
	$fsHead = "FS";
	$d2Head = "D2/D3";
	$unHead = "UN";
    }

    my $idHeadLen = length($idHead) * $mbFontWidth;
    my $fsHeadLen = length($fsHead) * $mbFontWidth;
    my $poHeadLen = length($poHead) * $mbFontWidth;
    my $d2HeadLen = length($d2Head) * $mbFontWidth;
    my $unHeadLen = length($unHead) * $mbFontWidth;

    my $idHeadX = ($xLeft + $x1) / 2.0 - $idHeadLen / 2.0;
    my $fsHeadX = ($x1 + $x2) / 2.0 - $fsHeadLen / 2.0;
    my $poHeadX = ($x2 + $x3) / 2.0 - $poHeadLen / 2.0;
    my $d2HeadX = ($x3 + $x4) / 2.0 - $d2HeadLen / 2.0;
    my $unHeadX = ($x4 + $xLeft + $graphWidth) / 2.0 - $unHeadLen / 2.0;

    my $yHeadPos = $yTop - 20;
    $im->string(gdMediumBoldFont, $idHeadX, $yHeadPos, $idHead, $black);
    $im->string(gdMediumBoldFont, $fsHeadX, $yHeadPos, $fsHead, $black);
    $im->string(gdMediumBoldFont, $poHeadX, $yHeadPos, $poHead, $black);
    $im->string(gdMediumBoldFont, $d2HeadX, $yHeadPos, $d2Head, $black);
    $im->string(gdMediumBoldFont, $unHeadX, $yHeadPos, $unHead, $black);

    # Plot the histogram
    my $binWid = $binSize * $graphWidth / $xMax;
    $binNo = 0;
    $binX = $binNo * $binSize;
    my $prevXpos = 0;
    while ($binX < $xMax) {
	my $cnt = $binDupPairs[$binNo];

	my $dpCnt = $dpPairs[$binNo];
	my $mtCnt = $mtPairs[$binNo];
	my $poCnt = $p0Pairs[$binNo];
	my $fsCnt = $p1Pairs[$binNo];
	my $d2Cnt = $p2Pairs[$binNo];
	my $d3Cnt = $p3Pairs[$binNo];

	my $xPos = $xLeft + $binX * $graphWidth / $xMax + $binWid;
	my $yPos = $yBottom - $cnt * $graphHeight / $yMax;
	$yPos = $yTop if ($yPos < $yTop);

	my $dpPos = $yBottom - $dpCnt * $graphHeight / $yMax;
	my $mtPos = $dpPos - $mtCnt * $graphHeight / $yMax;
	my $poPos = $mtPos - $poCnt * $graphHeight / $yMax;
	my $fsPos = $poPos - $fsCnt * $graphHeight / $yMax;
	my $d2Pos = $fsPos - $d2Cnt * $graphHeight / $yMax;
	my $d3Pos = $d2Pos - $d3Cnt * $graphHeight / $yMax;

	$dpPos = $yTop if ($dpPos < $yTop);
	$mtPos = $yTop if ($mtPos < $yTop);
	$poPos = $yTop if ($poPos < $yTop);
	$fsPos = $yTop if ($fsPos < $yTop);
	$d2Pos = $yTop if ($d2Pos < $yTop);
	$d3Pos = $yTop if ($d3Pos < $yTop);

	if ($yBottom - $yPos > 0.5) {
	    $prevXpos = $xPos - $binWid unless ($prevXpos);
	    $im->filledRectangle($prevXpos, $yBottom, $xPos, $dpPos, $dpColor) if ($yBottom - $dpPos > 0.5);
	    $im->filledRectangle($prevXpos, $dpPos, $xPos, $mtPos, $mtColor) if ($dpPos - $mtPos > 0.5);
	    $im->filledRectangle($prevXpos, $mtPos, $xPos, $poPos, $poColor) if ($mtPos - $poPos > 0.5);
	    $im->filledRectangle($prevXpos, $poPos, $xPos, $fsPos, $fsColor) if ($poPos - $fsPos > 0.5);
	    $im->filledRectangle($prevXpos, $fsPos, $xPos, $d2Pos, $d2Color) if ($fsPos - $d2Pos > 0.5);
	    $im->filledRectangle($prevXpos, $d2Pos, $xPos, $d3Pos, $d3Color) if ($d2Pos - $d3Pos > 0.5);
	    $im->filledRectangle($prevXpos, $d3Pos, $xPos, $yPos, $unColor)  if ($d3Pos - $yPos > 0.5);

	    $im->rectangle($prevXpos, $yBottom, $xPos, $yPos, $color);
	}

	$binNo++;
	$binX = $binNo * $binSize;
	$prevXpos = $xPos;
    }

    ShowBarLegends($xLeft, $yBottom);
}


sub ShowBarLegends()
{
    my ($xLeft, $yBottom) = @_;

    my $dpLgd = "DP";
    my $mtLgd = "MT";
    my $poLgd = "PO";
    my $fsLgd = "FS";
    my $d2Lgd = "D2";
    my $d3Lgd = "D3";
    my $unLgd = "UN";

    my $dpStrLen = length($dpLgd);
    my $mtStrLen = length($mtLgd);
    my $poStrLen = length($poLgd);
    my $fsStrLen = length($fsLgd);
    my $d2StrLen = length($d2Lgd);
    my $d3StrLen = length($d3Lgd);

    my $dpLen = $mbFontWidth * $dpStrLen;
    my $mtLen = $mbFontWidth * $mtStrLen;
    my $poLen = $mbFontWidth * $poStrLen;
    my $fsLen = $mbFontWidth * $fsStrLen;
    my $d2Len = $mbFontWidth * $d2StrLen;
    my $d3Len = $mbFontWidth * $d3StrLen;

    my $lgdGap = 30;
    my $strGap = 7;
    my $lgdWd = 20;
    my $lgdHt = 10;
    my $lgxPos = $xLeft + 15;
    my $lgyPos = $yBottom + 80;

    my $noteStr = "Note: colored bars represent related subjects derived from the pedigree and SSM file";
    $im->string(gdMediumBoldFont, $xLeft, $lgyPos, $noteStr, $black);
    $lgyPos += 25;
    my $stryPos = $lgyPos - 2;

    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $dpColor);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $dpLgd, $black);

    $lgxPos += $dpLen + $lgdWd + $lgdGap;
    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $mtColor);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $mtLgd, $black);

    $lgxPos += $mtLen + $lgdWd + $lgdGap;
    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $poColor);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $poLgd, $black);

    $lgxPos += $poLen + $lgdWd + $lgdGap;
    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $fsColor);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $fsLgd, $black);

    $lgxPos += $fsLen + $lgdWd + $lgdGap;
    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $d2Color);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $d2Lgd, $black);

    $lgxPos += $d2Len + $lgdWd + $lgdGap;
    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $d3Color);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $d3Lgd, $black);

    $lgxPos += $d3Len + $lgdWd + $lgdGap;
    $im->filledRectangle($lgxPos, $lgyPos, $lgxPos+$lgdWd, $lgyPos+$lgdHt, $unColor);
    $im->string(gdMediumBoldFont, $lgxPos+$lgdWd+$strGap, $stryPos, $unLgd, $black);
}

#
# Scatter plot to show distribution of both HGMR and AGMR
#
sub PlotScatterPlot
{
    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yTop    = $topEdge + $mlFontWidth * 3;
    my $yBottom = $topEdge + $graphHeight + $mlFontWidth * 3;

    # Plot axes
    my ($xSt, $xStep, $numXt, $xMinors) = (0, 5, 5, 5);
    my ($ySt, $yStep, $numYt, $yMinors) = (20, 5, 8, 5);

    my $maxHgmr = 28;
    $numXt = 6 if ($maxHgmr > 25);
    $numXt = 7 if ($maxHgmr > 30);
    $numXt = 8 if ($maxHgmr > 35);

#    $im->string(gdMediumBoldFont, $xLeft, $yTop + 50, "HF (half sibling + first cousin)", $black);

    my $xRange = $xStep * $numXt;
    my $yRange = $yStep * $numYt;

    my $xLeft = $leftEdge;
    my $yTop = $yBottom - $graphHeight;
    my $yPos = $yTop + 50;

    my $graphTitle = "Distribution of HGMR and AGMR of Pairs of Closely Related Subjects";
    $yTop = PlotGraphTitle($graphTitle, $xLeft, $yTop);
    $yBottom = $yTop + $graphHeight;

    PlotScatterAxes($yTop, $xSt, $xStep, $numXt, $xMinors, $ySt, $yStep, $numYt, $yMinors);

    # Draw the scatter plot
    foreach my $pair (keys %relPairs) {
	my %info = %{$relPairs{$pair}};
	my $smp1 = $info{smp1};
	my $smp2 = $info{smp2};
	my $hgmr = $info{hgmr};
	my $agmr = $info{agmr};
	my $pRel = $info{pedRel};

	my $color = $unScatColor;
	if    ($pRel eq "PO") { $color = $poColor; }
	elsif ($pRel eq "DP") { $color = $mtColor; }
	elsif ($pRel eq "MT") { $color = $mtColor; }
	elsif ($pRel eq "FS") { $color = $fsColor; }
	elsif ($pRel eq "HF") { $color = $hfColor; }
	elsif ($pRel eq "CS") { $color = $d3ScatColor; }
	elsif ($pRel eq "FC") { $color = $d3ScatColor; }
	elsif ($pRel eq "HA") { $color = $d3ScatColor; }
	elsif ($pRel eq "HS") { $color = $d2Color; }
	elsif ($pRel eq "AU") { $color = $d2Color; }
	elsif ($pRel eq "AV") { $color = $d2Color; }
	elsif ($pRel eq "GP") { $color = $d2Color; }
	elsif ($pRel eq "HC") { $color = $hcColor; }

	my $xPos = $xLeft + ($hgmr - $xSt) * $graphWidth / $xRange;
	my $yPos = $yBottom - ($agmr - $ySt) * $graphHeight / $yRange;
	if ($xPos + $dotSize  < $xLeft + $graphWidth && $yPos > $yTop + $dotSize) {
	    if ($pRel ne "HF") {
		PlotOneDot($xPos, $yPos, $dotSize, $color);
	    }
	    else {
		PlotOneDot($xPos, $yPos, $hfdSize, $color);
	    }
	}
    }

    if ($hasHF) {
	my $dx = 10;
	my $dy = 15;
	my $gapx = 15;
	my $gapy = 6;
	$im->arc($xLeft + $dx, $yTop + $dy, $hfdSize, $hfdSize, 0, 360, $hfColor);
	$im->fillToBorder($xLeft + $dx, $yTop + $dy, $hfColor, $hfColor);
	$im->string(gdMediumBoldFont, $xLeft + $dx + $gapx, $yTop + $dy - $gapy, "HF (half sibling + first cousin)", $black);
    }

    # Plot contour lines
    my $pop95 = 0.0553867;  # The area of P > $pop95 contains 95% of the pairs
    my $fsp95 = 0.0051547;
    my $d2p95 = 0.0011624;
    my $d3p95 = 0.0003346;
    my $pop99 = 0.0095202;
    my $fsp99 = 0.0009973;
    my $d2p99 = 0.0000356;
    my $d3p99 = 0.0000032;

    PlotPoContour($yTop, $pop95, $xSt, $ySt, $xRange, $yRange, $poContour);
    PlotIsoProbCircle($yTop, "FS", $fsp95, $xSt, $ySt, $xRange, $yRange, $fsContour);
    PlotIsoProbCircle($yTop, "D2", $d3p95, $xSt, $ySt, $xRange, $yRange, $d2Contour);
    PlotIsoProbCircle($yTop, "D3", $d3p95, $xSt, $ySt, $xRange, $yRange, $d3Contour);

    # Plot legends
    my ($poLgdx, $fsLgdx, $d2Lgdx, $d3Lgdx) = (.5, 10, 14, 18);
    my ($poLgdy, $fsLgdy, $d2Lgdy, $d3Lgdy) = (22, 22, 24, 26);
    my $poxPos = $xLeft + ($poLgdx - $xSt) * $graphWidth / $xRange;
    my $poyPos = $yBottom - ($poLgdy - $ySt) * $graphHeight / $yRange;
    my $fsxPos = $xLeft + ($fsLgdx - $xSt) * $graphWidth / $xRange;
    my $fsyPos = $yBottom - ($fsLgdy - $ySt) * $graphHeight / $yRange;
    my $d2xPos = $xLeft + ($d2Lgdx - $xSt) * $graphWidth / $xRange;
    my $d2yPos = $yBottom - ($d2Lgdy - $ySt) * $graphHeight / $yRange;
    my $d3xPos = $xLeft + ($d3Lgdx - $xSt) * $graphWidth / $xRange;
    my $d3yPos = $yBottom - ($d3Lgdy - $ySt) * $graphHeight / $yRange;

    my $poLgd = "Parent-offspring";
    my $fsLgd = "Full sibling";
    my $d2Lgd = "2nd degree";
    my $d3Lgd = "3rd degree";

    my $poStrLen = length($poLgd);
    my $fsStrLen = length($fsLgd);
    my $d2StrLen = length($d2Lgd);
    my $d3StrLen = length($d3Lgd);

    my $poLen = $mbFontWidth * $poStrLen;
    my $fsLen = $mbFontWidth * $fsStrLen;
    my $d2Len = $mbFontWidth * $d2StrLen;
    my $d3Len = $mbFontWidth * $d3StrLen;

    $im->string(gdMediumBoldFont, $poxPos, $poyPos, $poLgd, $black);
    $im->string(gdMediumBoldFont, $fsxPos, $fsyPos, $fsLgd, $black);
    $im->string(gdMediumBoldFont, $d2xPos, $d2yPos, $d2Lgd, $black);
    $im->string(gdMediumBoldFont, $d3xPos, $d3yPos, $d3Lgd, $black);

    DrawArrow($yTop, $poLgdx, $poLgdy, $poStrLen, $poMux, $poMuy, 0.4, $xSt, $xRange, $ySt, $yRange, $poContour);
    DrawArrow($yTop, $fsLgdx, $fsLgdy, $fsStrLen, $fsMux, $fsMuy, 0.4, $xSt, $xRange, $ySt, $yRange, $fsContour);
    DrawArrow($yTop, $d2Lgdx, $d2Lgdy, $d2StrLen, $d2Mux, $d2Muy, 0.3, $xSt, $xRange, $ySt, $yRange, $d2Contour);
    DrawArrow($yTop, $d3Lgdx, $d3Lgdy, $d3StrLen, $d3Mux, $d3Muy, 0.3, $xSt, $xRange, $ySt, $yRange, $poContour);

    my $lgyPos1 = $yBottom + 80;
    my $lgyPos2 = $lgyPos1 + 15;
    my $noteStr1 = "Note: contour line of each relationship type shows the area that is";
    my $noteStr2 = "      expected to contain 95% of the subject pairs of this type";
    $im->string(gdMediumBoldFont, $xLeft, $lgyPos1, $noteStr1, $black);
    $im->string(gdMediumBoldFont, $xLeft, $lgyPos2, $noteStr2, $black);
}

#
# Plot an arrow
#
sub DrawArrow
{
    my ($yTop, $gx, $gy, $gl, $cx, $cy, $r2, $xSt, $xRange, $ySt, $yRange, $color) = @_;

    my $xLeft   = $leftEdge;
    my $yBottom = $yTop + $graphHeight;

    my $x1 = $gx;
    my $y1 = $gy;
    my $x2 = $cx;
    my $y2 = $cy;

    my $x1P = $xLeft + ($x1 - $xSt) * $graphWidth / $xRange + $mbFontWidth * $gl / 2 ;
    my $x2P = $xLeft + ($x2 - $xSt) * $graphWidth / $xRange;
    my $y1P = $yBottom - ($y1 - $ySt) * $graphHeight / $yRange + $mbFontHeight / 2 ;
    my $y2P = $yBottom - ($y2 - $ySt) * $graphHeight / $yRange;

    my $r1 = 0.05;

    my $dx = $x1P - $x2P;
    my $dy = $y1P - $y2P;

    $x1P = $x1P - $dx * $r1;
    $y1P = $y1P - $dy * $r1;
    $x2P = $x2P + $dx * $r2;
    $y2P = $y2P + $dy * $r2;

    $im->line($x1P, $y1P, $x2P, $y2P, $color);

    my $da = 10 * $pi / 180;
    my $aLen = 12;
    my $alpha = atan2($y2P - $y1P, $x2P - $x1P);
    my $ax1 = $x2P - $aLen * cos($alpha + $da);
    my $ay1 = $y2P - $aLen * sin($alpha + $da);
    my $ax2 = $x2P - $aLen * cos($alpha - $da);
    my $ay2 = $y2P - $aLen * sin($alpha - $da);
    $im->line($x2P, $y2P, $ax1, $ay1, $color);
    $im->line($x2P, $y2P, $ax2, $ay2, $color);
    $im->line($ax1, $ay1, $ax2, $ay2, $color);
}

#
# Plot one contour line with same probability density
#
sub PlotIsoProbCircle
{
    my ($yTop, $type, $pval, $xSt, $ySt, $xRange, $yRange, $cirColor) = @_;

    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yBottom = $yTop + $graphHeight;

    my ($minHgmr, $maxHgmr) = (0, 10);
    if ($type =~ /D2/i) {
	($minHgmr, $maxHgmr) = (3, 23);
    }
    elsif ($type =~ /D3/i) {
	($minHgmr, $maxHgmr) = (8, 32);
    }

    my $x = $minHgmr;

    my ($prevx, $prevy1, $prevy2) = (0, 0, 0);
    while ($x < $maxHgmr) {
	my ($hasSol, $agmr1, $agmr2) = GetAgmrGivenHgmr($x, $pval, $type);
	if ($hasSol) {
	    my $xPos = $xLeft + ($x - $xSt) * $graphWidth / $xRange;
	    my $yPos1 = $yBottom - ($agmr1 - $ySt) * $graphHeight / $yRange;
	    my $yPos2 = $yBottom - ($agmr2 - $ySt) * $graphHeight / $yRange;
	    if ($prevx) {
		$im->line($prevx, $prevy1, $xPos, $yPos1, $cirColor);
		$im->line($prevx, $prevy2, $xPos, $yPos2, $cirColor);
	    }
	    else {
		$im->line($xPos, $yPos2, $xPos, $yPos1, $cirColor);
	    }
	    $prevx = $xPos;
	    $prevy1 = $yPos1;
	    $prevy2 = $yPos2;
	}
	elsif ($prevx) {
	    $im->line($prevx, $prevy1, $prevx, $prevy2, $cirColor);
	    $x = $maxHgmr;
	}

	$x += 0.01;
    }
}

#
# Plot one contour line with same probability density for PO pairs
#
sub PlotPoContour
{
    my ($yTop, $pval, $xSt, $ySt, $xRange, $yRange, $cirColor) = @_;

    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yBottom = $yTop + $graphHeight;

    my ($minHgmr, $maxHgmr) = (0, 3);
    my $x = $minHgmr;

    my ($prevx, $prevy1, $prevy2) = (0, 0, 0);
    while ($x < $maxHgmr) {
	my ($hasSol, $agmr1, $agmr2) = GetPoAgmrGivenHgmr($x, $pval);
	if ($hasSol) {
	    my $xPos = $xLeft + ($x - $xSt) * $graphWidth / $xRange;
	    my $yPos1 = $yBottom - ($agmr1 - $ySt) * $graphHeight / $yRange;
	    my $yPos2 = $yBottom - ($agmr2 - $ySt) * $graphHeight / $yRange;
	    if ($prevx) {
		$im->line($prevx, $prevy1, $xPos, $yPos1, $cirColor);
		$im->line($prevx, $prevy2, $xPos, $yPos2, $cirColor);
	    }

	    $prevx = $xPos;
	    $prevy1 = $yPos1;
	    $prevy2 = $yPos2;
	}
	elsif ($prevx) {
	    $im->line($prevx, $prevy1, $prevx, $prevy2, $cirColor);
	    $x = $maxHgmr;
	}

	$x += 0.01;
    }
}

#
# Calculate y value when x and p are known
#
sub GetAgmrGivenHgmr
{
    my ($hgmr, $p, $type) = @_;

    my $isD2 = $type =~ /D2/i ? 1 : 0;
    my $isD3 = $type =~ /D3/i ? 1 : 0;

    my ($rho, $rho2, $sigx, $sigy, $mux, $muy) = (0, 0, 0, 0, 0, 0);

    if ($isD2) {
	$sigx = $d2Sx;
	$sigy = $d2Sy;
	$mux  = $d2Mux;
	$muy  = $d2Muy;
	$rho  = $d2Cor;
    }
    elsif ($isD3) {
	$sigx = $d3Sx;
	$sigy = $d3Sy;
	$mux  = $d3Mux;
	$muy  = $d3Muy;
	$rho  = $d3Cor;
    }
    else {
	$sigx = $fsSx;
	$sigy = $fsSy;
	$mux  = $fsMux;
	$muy  = $fsMuy;
	$rho  = $fsCor;
    }

    my $rho2 = 1 - $rho**2;
    my $c = 1 / (2 * $pi * $sigx * $sigy * sqrt($rho2));

    my $xu1 = $hgmr - $mux;
    my $delta = $rho2 * $sigy**2 * (2 * log($c/$p) -  $xu1**2 / $sigx**2);

    my $hasSol = 0;
    my ($agmr1, $agmr2) = (0, 0);

    if ($delta > 0) {
	$hasSol = 1;
	$agmr1 = $muy + $rho * $xu1 * $sigy / $sigx - sqrt($delta);
	$agmr2 = $muy + $rho * $xu1 * $sigy / $sigx + sqrt($delta);
    }
	
    return ($hasSol, $agmr1, $agmr2);
}

#
# Calculate y value for PO pairs when x and p are known
#
sub GetPoAgmrGivenHgmr
{
    my ($hgmr, $p) = @_;

    my $x = $hgmr;
    my $k = $poHgmrK;

    my $delta = -1*$k*$x - log($p * $poSy * sqrt(2 * $pi) / $k);

    my $hasSol = 0;
    my ($agmr1, $agmr2) = (0, 0);

    if ($delta > 0) {
	$hasSol = 1;
	$agmr1 = $poMuy - $poSy * sqrt(2*$delta);
	$agmr2 = $poMuy + $poSy * sqrt(2*$delta);
    }
	
    return ($hasSol, $agmr1, $agmr2);
}

#
# Plot graph title
#
sub PlotGraphTitle
{
    my ($graphTitle, $xLeft, $yTop) = @_;

    my $graphTotWidth = $graphWidth;

    my $numLetters = length($graphTitle);
    my $strLen = $mgFontWidth * $numLetters;
    my $xPos = $xLeft + ($graphWidth - $strLen) / 2;
    $im->string(gdGiantFont, $xPos, $yTop, $graphTitle, $black);

    $yTop += $mgFontHeight + 30;

    return $yTop;
}

#
# Plot axes for a histogram
#
sub PlotHistoAxes
{
    my ($yTop, $yMax, $xMin, $xMax, $gType) = @_;

    # Plot y-axis
    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yBottom = $yTop + $graphHeight;
    my $yMin = 0;

    $im->line($xLeft, $yTop, $xLeft, $yBottom, $black);
    $im->line($xRight, $yTop, $xRight, $yBottom, $black);

    # Plot ticks
    my ($numDigits, $firstDigit, $unit) = GetFirstDigit($yMax);
    my ($step, $numTicks, $minorTicks) = GetTicks($firstDigit, $unit);

    my $fullVal = $step * $numTicks;
    my $tickLblLen = length($fullVal);

    my $numTicks = int(($fullVal - $yMin) / $step + 0.5);
    my $minorTicks = 5;

    my $xMajor = $xLeft - $majorTickLen;
    my $dy = $graphHeight/$numTicks;
    my $minorDy = $graphHeight/($numTicks*$minorTicks);

    my $tickLblLen = 5;
    my $xString = $xLeft - $mbFontWidth*$tickLblLen - $majorTickLen - $tickGap;

    for my $i (0 .. $numTicks) {
	my $val = $yMin + $step * $i;
	my $valStr = sprintf("%5d", $val);
	my $y = int($yBottom - $dy*$i);
	$im->line($xLeft, $y, $xMajor, $y, $black);
	$im->string(gdMediumBoldFont, $xString, $y-$mbFontHeight/2, $valStr, $black);
	
	# Plot minor ticks
	if ($i < $numTicks) {
	    my $yMinor = $y;
	    my $xMinor = $xLeft - $minorTickLen;
	    for my $j (1 .. $minorTicks-1) {
		$yMinor -= $minorDy;
		$im->line($xLeft, $yMinor, $xMinor, $yMinor, $black);
	    }
	}
    }

    $im->line($xLeft, $yBottom, $xRight, $yBottom, $black);

    # Draw y-axis label
    my $yLabelStr = "Number of Pairs";

    my $lblLen = length($yLabelStr);
    my $yLabel = $yBottom - ($graphHeight - $mbFontWidth*$lblLen)/2;
    $im->stringUp(gdMediumBoldFont, $xString-20, $yLabel, $yLabelStr, $black);

    # Plot x-axis
    $im->line($xLeft, $yBottom, $xRight, $yBottom, $black);
    $im->line($xLeft, $yTop, $xRight, $yTop, $black);

    my $numXticks = int(($xMax - $xMin) / $xStep + 0.5);

    my $minorDx = $graphWidth/($numXticks*$xMinorTicks);
    for my $i (0 .. $numXticks) {
	my $x = $xLeft + $i * $graphWidth / $numXticks;
	$im->line($x, $yBottom, $x, $yBottom+$majorTickLen, $black);

	# Plot minor ticks
	if ($i < $numXticks) {
	    my $xMinor = $x;
	    for my $j (1 .. $xMinorTicks-1) {
		$xMinor += $minorDx;
		$im->line($xMinor, $yBottom, $xMinor, $yBottom+$minorTickLen, $black);
	    }
	}
	
	my $xVal = sprintf("%d", $xMin + $i * $xStep);
	$xVal = sprintf("%2.1f", $xMin + $i * $xStep) if ($xMax < 5);
	$tickLblLen = length($xVal);
	my $xLblPos = $x - $mbFontWidth * $tickLblLen / 2;
	my $yLblPos = $yBottom + 12;
	$im->string(gdMediumBoldFont, $xLblPos, $yLblPos, $xVal, $black);
    }

    my $xAxis = "Homozygous Genotype Mismatch Rate (%)";
    $xAxis = "All Genotype Mismatch Rate (%)" if ($gType == 2);

    my $strLetters = length($xAxis);
    my $strLen = $mbFontWidth * $strLetters;
    my $xPos = $xLeft + ($graphWidth - $strLen) / 2;
    $im->string(gdMediumBoldFont, $xPos, $yBottom + 40, $xAxis, $black);

    my $realXmax = $numXticks * $xStep;

    return ($realXmax, $fullVal);
}

#
# Plot axes for a scatter plot
#
sub PlotScatterAxes
{
    my ($yTop, $xMin, $xStep, $numXmajor, $numXminor, $yMin, $yStep, $numYmajor, $numYminor) = @_;

    # Plot y-axis
    my $xLeft   = $leftEdge;
    my $xRight  = $leftEdge + $graphWidth;
    my $yBottom = $yTop + $graphHeight;

    $im->line($xLeft, $yTop, $xLeft, $yBottom, $black);
    $im->line($xRight, $yTop, $xRight, $yBottom, $black);

    my $xMax = $xMin + $xStep * $numXmajor;
    my $yMax = $yMin + $yStep * $numYmajor;
    my $tickLblLen = 2;

    my $xMajor = $xLeft - $majorTickLen;
    my $dy = $graphHeight/$numYmajor;
    my $minorDy = $graphHeight/($numYmajor*$numYminor);

    my $tickLblLen = 5;
    my $xString = $xLeft - $mbFontWidth*$tickLblLen - $majorTickLen - $tickGap;

    for my $i (0 .. $numYmajor) {
	my $val = $yMin + $yStep * $i;
	my $valStr = sprintf("%5d", $val);
	my $y = int($yBottom - $dy*$i);
	$im->line($xLeft, $y, $xMajor, $y, $black);
	$im->string(gdMediumBoldFont, $xString, $y-$mbFontHeight/2, $valStr, $black);
	
	# Plot minor ticks
	if ($i < $numYmajor) {
	    my $yMinor = $y;
	    my $xMinor = $xLeft - $minorTickLen;
	    for my $j (1 .. $numYminor-1) {
		$yMinor -= $minorDy;
		$im->line($xLeft, $yMinor, $xMinor, $yMinor, $black);
	    }
	}
    }

    $im->line($xLeft, $yBottom, $xRight, $yBottom, $black);

    # Draw y-axis label
    my $yLabelStr = "All Genotype Mismatch Rate (%)";

    my $lblLen = length($yLabelStr);
    my $yLabel = $yBottom - ($graphHeight - $mbFontWidth*$lblLen)/2;
    $im->stringUp(gdMediumBoldFont, $xString-20, $yLabel, $yLabelStr, $black);

    # Plot x-axis
    $im->line($xLeft, $yBottom, $xRight, $yBottom, $black);
    $im->line($xLeft, $yTop, $xRight, $yTop, $black);

    my $minorDx = $graphWidth/($numXmajor*$numXminor);
    for my $i (0 .. $numXmajor) {
	my $x = $xLeft + $i * $graphWidth / $numXmajor;
	$im->line($x, $yBottom, $x, $yBottom+$majorTickLen, $black);

	# Plot minor ticks
	if ($i < $numXmajor) {
	    my $xMinor = $x;
	    for my $j (1 .. $numXminor-1) {
		$xMinor += $minorDx;
		$im->line($xMinor, $yBottom, $xMinor, $yBottom+$minorTickLen, $black);
	    }
	}
	
	my $xVal = sprintf("%d", $xMin + $i * $xStep);
	$xVal = sprintf("%2.1f", $xMin + $i * $xStep) if ($xMax < 5);
	$tickLblLen = length($xVal);
	my $xLblPos = $x - $mbFontWidth * $tickLblLen / 2;
	my $yLblPos = $yBottom + 12;
	$im->string(gdMediumBoldFont, $xLblPos, $yLblPos, $xVal, $black);
    }

    my $xAxis = "Homozygous Genotype Mismatch Rate (%)";
    my $strLetters = length($xAxis);
    my $strLen = $mbFontWidth * $strLetters;
    my $xPos = $xLeft + ($graphWidth - $strLen) / 2;
    $im->string(gdMediumBoldFont, $xPos, $yBottom + 40, $xAxis, $black);
}

#
# Plot a filled circle
#
sub PlotOneDot
{
    my ($x, $y, $size, $color) = @_;

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
# Find the first digit given a number.  Used for determining ticks and tick lengths
#
sub GetFirstDigit
{
    use integer;

    my $chkNum = shift;
    my $num = $chkNum;

    die "Number '$num' is not a positive integer!\n" if ($num < 1);

    my $unit = 1;
    my $numDigits = 1;
    while ($num > 44) {
	$num = $num / 10;
	$unit = $unit * 10;
	$numDigits++;
    }

    $numDigits++ if ($num > 9);

    $num++ if ($chkNum > $num * $unit);

    return ($numDigits, $num, $unit);
}

#
# Calculate ticks and tick lengths
#
sub GetTicks
{
    my ($firstDigit, $unit) = @_;

    my $numTicks = 0;
    my $minorTicks = 5;
    my $step = 0;

    if ($firstDigit <= 5) {
	$step = $unit;
	$numTicks = 5;
    }
    elsif ($firstDigit <= 10) {
	$step = $unit;
	$numTicks = 10;
    }
    elsif ($firstDigit <= 15) {
	$step = $unit;
	$numTicks = 15;
    }
    elsif ($firstDigit <= 20) {
	$step = $unit * 5;
	$numTicks = 4;
	$minorTicks = 2;
    }
    elsif ($firstDigit <= 25) {
	$step = $unit * 5;
	$numTicks = 5;
	$minorTicks = 10;
    }
    elsif ($firstDigit <= 30) {
	$step = $unit * 5;
	$numTicks = 6;
    }
    elsif ($firstDigit < 50) {
	$step = $unit * 10;
	$numTicks = 5;
	$minorTicks = 10;
    }

    if ($unit == 1) {
	#$minorTicks = 1;
	$numTicks = 5 if ($firstDigit < 4);
    }

    return ($step, $numTicks, $minorTicks);
}

#
# Read parameters from command line
#
sub GetParameters
{
    my @args = @_;

    my ($inFile, $outFile) = ("", "");
    my $graphType = 0;
    my %params = ();

    if (@args > 2) {
	$inFile = $args[0];
	$outFile = $args[1];
	$graphType = $args[2];

	for my $i (3 .. $#args) {
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

    if ($inFile =~ /^-/ || $outFile =~ /^-/ || $graphType =~ /^-/) {
	($inFile, $outFile) = ("", "");
    }
    else {
	if ($graphType && $graphType ne "1" && $graphType ne "2" && $graphType ne "3") {
	    print "ERROR: Invalid graph type \"$graphType\".  Expected 1, 2 or 3.\n";
	}
    }

    return ($inFile, $outFile, $graphType, \%params);
}

#
# Read related pairs from the input file
#
sub ReadRelatedPairs
{
    my $relFile = shift;

    my %dups = ();
    my %rels = ();
    my $hasMultiDs = 0;
    my $error = "";

    open FILE, $relFile or die "Couldn't open $relFile!\n";
    my $header = <FILE>;
    while ($header !~ /\S/ || $header =~ /^\#/) {
	if    ($header =~ /FS HGMR: (\d+\.\d+)/) {
	    $fsMux = $1;
	}
	elsif ($header =~ /D2 HGMR: (\d+\.\d+)/) {
	    $d2Mux = $1;
	}
	elsif ($header =~ /D3 HGMR: (\d+\.\d+)/) {
	    $d3Mux = $1;
	}
	elsif ($header =~ /PO AGMR: (\d+\.\d+)/) {
	    $poMuy = $1;
	}
	elsif ($header =~ /FS AGMR: (\d+\.\d+)/) {
	    $fsMuy = $1;
	}
	elsif ($header =~ /D2 AGMR: (\d+\.\d+)/) {
	    $d2Muy = $1;
	}
	elsif ($header =~ /D3 AGMR: (\d+\.\d+)/) {
	    $d3Muy = $1;
	}

	$header = <FILE>;
    }

    chomp $header;
    if ($header =~ /^DS1\tsample1\tDS2\tsample2\tsubject1\tsubject2\tsex1\tsex2\tHG match\tHG miss\tHGMR\tAG match\tAG miss\tAGMR\tgeno relation\tped relation\tp_value/) {
	$hasMultiDs = 1;
    }
    elsif ($header !~ /^sample1\tsample2\tsubject1\tsubject2\tsex1\tsex2\tHG match\tHG miss\tHGMR\tAG match\tAG miss\tAGMR\tgeno relation\tped relation\tp_value/) {
	$error = "Invalid file for related pairs.  Expected the following columns:\n";
	$error .= "sample1\nsample2\nsubject1\nsubject2\nsex1\nsex2\n";
	$error .= "HG match\nHG miss\nHGMR\nAG match\nAG miss AGMR\ngeno relation\nped relation\np_value\n";

	return (\%dups, \%rels, $error);
    }

    while(<FILE>) {
	chomp;
	my ($ds1, $smp1, $ds2, $smp2, $sbj1, $sbj2, $sex1, $sex2, $hgSame, $hgDiff, $hgmr, $agSame, $agMiss, $agmr, $genoRel, $pedRel, $pval) = ();
	if ($hasMultiDs) {
	    ($ds1, $smp1, $ds2, $smp2, $sbj1, $sbj2, $sex1, $sex2, $hgSame, $hgDiff, $hgmr, $agSame, $agMiss, $agmr, $genoRel, $pedRel, $pval)
	= split /\t/, $_;
	    $smp1 = $ds1 . "___" . $smp1;
	    $smp2 = $ds2 . "___" . $smp2;
	}
	else {
	    ($smp1, $smp2, $sbj1, $sbj2, $sex1, $sex2, $hgSame, $hgDiff, $hgmr, $agSame, $agMiss, $agmr, $genoRel, $pedRel, $pval)
	= split /\t/, $_;
	}

	my $pair = "$smp1\t$smp2";
	my %info = (smp1 => $smp1, smp2 => $smp2, sbj1 => $sbj1, sex1 => $sex1, sex2 => $sex2, hgmr => $hgmr, agmr => $agmr,
		    genoRel => $genoRel, pedRel => $pedRel, pval => $pval);
	if ($genoRel eq "ID" || ($pedRel && $pedRel =~ /(DP)|(MT)/)) {
	    $dups{$pair} = \%info;
	}
	
	unless ($genoRel eq "ID") {
	    $rels{$pair} = \%info;
	}
    }
    close FILE;

    return (\%dups, \%rels, $error);
}
