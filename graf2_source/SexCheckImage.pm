# $Id$
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
# This software/database is a "United States Government Work" under the
# terms of the United States Copyright Act. It was written as part of
# the author's official duties as a United States Government employee and
# thus cannot be copyrighted. This software/database is freely available
# to the public for use. The National Library of Medicine and the U.S.
# Government have not placed any restriction on its use or reproduction.
#
# Although all reasonable efforts have been taken to ensure the accuracy
# and reliability of the software and data, the NLM and the U.S.
# Government do not and cannot warrant the performance or results that
# may be obtained by using this software or data. The NLM and the U.S.
# Government disclaim all warranties, express or implied, including all
# warranties of performance, merchantability or fitness for any particular
# purpose.
#
# Please cite the author in any work or product based on this material.
#
# ===========================================================================
package SexCheckImage;

use strict;
use warnings;
use CGI qw(:standard);
use GD;
use GD::Text;
use GD::Graph;
use GD::Graph::colour;
use GraphYAxisTicks;

sub new
{
    my ($class) = @_;

    my $self = {
        leftEdge => 80,      # The left edge of the graph
        graphWidth => 1000,  # Width of the plot (rectangle box only)
        graphHeight => 600,  # Height of the plot

        rightEdge => 20,     # The right edge
        topEdge => 10,       # The top edge
        bottomEdge => 100,   # The bottom edge

        axisLabelGap => 20,  # Distance between tick labels and axis labels
        tickGap => 5,        # Gap between the tick and tick label
        graphGap => 100,

        majorTickLen => 6,   # major tick length (in pixels)
        minorTickLen => 3,   # minor tick length
        graphTitleHt => 20,  # Height of the graph title
        xMinorTicks => 10,   # number of minor ticks on x-axis

        maxBinSmps => 0,     # max number of samples display on y-axis
    };

    bless ($self, $class);

    return $self;
}

sub SetImageSize
{
    my $self = shift;

    # Width and height of the whole graph
    $self->{titleHeight} = gdLargeFont->height * 5;
    $self->{imageWidth}  = $self->{leftEdge} + $self->{graphWidth} + $self->{rightEdge},
    $self->{imageHeight} = $self->{topEdge} + $self->{titleHeight} + $self->{graphHeight} + $self->{bottomEdge},
}

sub SetParameters
{
    my ($self, $params) = @_;

    foreach my $key (keys %$params) {
        if ($key && defined $self->{$key}) {
            $self->{$key} = $params->{$key};
        }
    }
}

sub PlotGraph
{
    my ($self, $sexChkRes, $pngFile) = @_;

    my $mbFontWidth  = gdMediumBoldFont->width;
    my $mbFontHeight = gdMediumBoldFont->height;
    my $mlFontWidth  = gdLargeFont->width;
    my $mlFontHeight = gdLargeFont->height;

    $self->SetImageSize();

    my $graphWidth  = $self->{graphWidth};
    my $graphHeight = $self->{graphHeight};

    my $img = new GD::Image($self->{imageWidth}, $self->{imageHeight});

    my $greyNum = 150;

    my $black  = $img->colorAllocate(  0,   0,   0);
    my $white  = $img->colorAllocate(255, 255, 255);

    my $maleColor   = $img->colorAllocate(  0,   0, 255); # Blue
    my $femaleColor = $img->colorAllocate(250,   0,   0); # Red
    my $cutoffColor = $img->colorAllocate(  0, 255, 255); # Cyan
    my $otherColor  = $img->colorAllocate($greyNum, $greyNum, $greyNum); # Grey

    my $xLeft   = $self->{leftEdge};
    my $yTop    = $self->{topEdge}  + $self->{titleHeight};
    my $xRight  = $xLeft + $graphWidth;
    my $yBottom = $yTop  + $graphHeight;

    my $yMax = $self->{maxBinSmps} ? $self->{maxBinSmps} :  $sexChkRes->{maxBinSmps};
    my $xMax = 100;

    my @g12Smps = @{$sexChkRes->{g12Smps}};
    my @g12Males = @{$sexChkRes->{g12Males}};
    my @g12Females = @{$sexChkRes->{g12Females}};

    open(IMG, ">$pngFile") or die $!;

    binmode IMG;

    $img->filledRectangle(0, $self->{imageHeight}, $self->{imageWidth}, 0, $white);

    my $graphTitle = "Comparison of subject phenotype and genotype sexes (determined by GRAF)";

    my $titleLen = length($graphTitle);
    my $xTitle = $xLeft + ($self->{graphWidth} - $titleLen*$mlFontWidth) / 2;
    my $yTitle = $yTop - $mlFontHeight * 4;

    $img->string(gdLargeFont, $xTitle, $yTitle, $graphTitle, $black);

    my $lgdTop = 0;
    ($yMax, $lgdTop) = $self->PlotAxes($img, $xLeft, $xRight, $yTop, $yBottom, $yMax, $black);
    $lgdTop += 20;

    # Show 'unknown sex' legend only if there are samples with phenotype sexes unknown
    my $showOther = $sexChkRes->{numUnkSexes} ? 1 : 0;
    $self->PlotLegends($img, $xLeft, $lgdTop, $black, $maleColor, $femaleColor, $otherColor, $showOther);

    my $xMaleCut    = $xLeft + $sexChkRes->{maxMaleHetroRate}   * $graphWidth * 1.0 / $xMax;
    my $xFemaleCut1 = $xLeft + $sexChkRes->{minFemaleHetroRate} * $graphWidth * 1.0 / $xMax;
    my $xFemaleCut2 = $xLeft + $sexChkRes->{maxFemaleHetroRate} * $graphWidth * 1.0 / $xMax;
    $img->line($xMaleCut,   $yTop, $xMaleCut,   $yBottom, $cutoffColor);
    $img->line($xFemaleCut1, $yTop, $xFemaleCut1, $yBottom, $cutoffColor);
    $img->line($xFemaleCut2, $yTop, $xFemaleCut2, $yBottom, $cutoffColor);

    my $ySexLbl    = $yTop - 20;
    my $maleLen    = $mbFontWidth * length("Male");
    my $femaleLen  = $mbFontWidth * length("Female");
    my $xMaleLbl   = $xLeft + ($xMaleCut - $xLeft - $maleLen) / 2;
    my $xFemaleLbl = $xFemaleCut1 + ($xFemaleCut2 - $xFemaleCut1 - $femaleLen) / 2;
    $img->string(gdMediumBoldFont, $xMaleLbl,   $ySexLbl, "Male",   $black);
    $img->string(gdMediumBoldFont, $xFemaleLbl, $ySexLbl, "Female", $black);

    # Show min and max FP SNPs with genotypes per subject
    my $snpLgdEdge = 105; # Room to show "FP SNPs/Smp"
    $snpLgdEdge = 85 if ($graphWidth < 401);

    my $snpLgdx = $xLeft + $graphWidth - $snpLgdEdge;
    my $snpLgdy = $yTop + 10;
    my $snpLgdGap = $mbFontHeight * 1.5;

    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, "SNPs/Sample:", $black);

    $snpLgdx += $mbFontWidth;
    $snpLgdy += $snpLgdGap;
    my $minSnpText = sprintf(" Min %5d", $sexChkRes->{minXsnps});
    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $minSnpText, $black);

    $snpLgdy += $snpLgdGap;
    my $maxSnpText = sprintf(" Max %5d", $sexChkRes->{maxXsnps});
    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $maxSnpText, $black);

    $snpLgdy += $snpLgdGap;
    my $meanSnpText = sprintf("Mean %5.0f", $sexChkRes->{meanXsnps});
    $img->string(gdMediumBoldFont, $snpLgdx, $snpLgdy, $meanSnpText, $black);

    my $binSize = $sexChkRes->{binSize};
    my $binWid  = $binSize * $graphWidth / $xMax;
    my $binNo = 0;
    my $binX = $binNo * $binSize;
    my $prevXpos = $xLeft;
    while ($binX < $xMax) {
        my $cnt = $g12Smps[$binNo];
        my $males = $g12Males[$binNo];
        my $females = $g12Females[$binNo];
        my $others = $cnt - $males - $females;

        my $xPos = $xLeft + $binX * $graphWidth / $xMax + $binWid;
        my $yPos = $yBottom - $cnt * $graphHeight / $yMax;

        my $yPos1 = $yBottom - $males * $graphHeight / $yMax;
        my $yPos2 = $yBottom - ($males+$females) * $graphHeight / $yMax;

        my $color1 = $maleColor;
        my $color2 = $femaleColor;

        # Plot females on bottom in the area expected to be males
        if ($binX < 20) {
            $yPos1 = $yBottom - $females * $graphHeight / $yMax;
            $yPos2 = $yBottom - ($males+$females) * $graphHeight / $yMax;
            $color1 = $femaleColor;
            $color2 = $maleColor;
        }

        $yPos  = $yTop if ($yPos  < $yTop);
        $yPos1 = $yTop if ($yPos1 < $yTop);
        $yPos2 = $yTop if ($yPos2 < $yTop);

        $img->filledRectangle($prevXpos, $yBottom, $xPos, $yPos1, $color1)  if ($yBottom - $yPos1 > 0.5);
        $img->filledRectangle($prevXpos, $yPos1, $xPos, $yPos2, $color2)    if ($yPos1 - $yPos2 > 0.5);
        $img->filledRectangle($prevXpos, $yPos2, $xPos, $yPos, $otherColor) if ($yPos2 - $yPos > 0.5 && $others > 0);

        $img->rectangle($prevXpos, $yBottom, $xPos, $yPos, $black) if ($yBottom - $yPos > 0.5);

        $binNo++;
        $binX = $binNo * $binSize;
        $prevXpos = $xPos;
    }

    print IMG $img->png;

    print "Histogram saved to $pngFile\n\n";
}

sub PlotAxes
{
    my ($self, $img, $xLeft, $xRight, $yTop, $yBottom, $yMax, $lineColor) = @_;

    my $mbFontWidth  = gdMediumBoldFont->width;
    my $mbFontHeight = gdMediumBoldFont->height;
    my $graphHeight = $self->{graphHeight};
    my $graphWidth  = $self->{graphWidth};

    die "\nERROR calling PlotAxes: yMax value is 0 or NULL.\n" if (!$yMax || $yMax < 0);

    my $yMin = 0;

    # Plot y-axis
    $img->line($xLeft, $yTop, $xLeft, $yBottom, $lineColor);
    $img->line($xRight, $yTop, $xRight, $yBottom, $lineColor);

    # Plot ticks
    my $axisTicks = new GraphYAxisTicks();
    my $isManualScale = $self->{maxBinSmps} ? 1 : 0;
    my ($step, $numMajorTicks, $numMinorTicks) = $axisTicks->GetTicks($yMax, $isManualScale);

    # Use user defined ticks if there are any
    $step = $self->{yStep} if ($self->{yStep});

    my $fullVal = $step * $numMajorTicks;
    my $tickLblLen = 5;

    my $xMajor = $xLeft - $self->{majorTickLen};
    my $dy = $graphHeight/$numMajorTicks;
    my $minorDy = $graphHeight/($numMajorTicks*$numMinorTicks);

    my $xString = $xLeft - $mbFontWidth*$tickLblLen - $self->{majorTickLen} - $self->{tickGap};

    for my $i (0 .. $numMajorTicks) {
        my $val = $yMin + $step * $i;
        my $valStr = sprintf("%5d", $val);
        my $y = int($yBottom - $dy*$i);
        $img->line($xLeft, $y, $xMajor, $y, $lineColor);
        $img->string(gdMediumBoldFont, $xString, $y-$mbFontHeight/2, $valStr, $lineColor);

        # Plot minor ticks
        if ($i < $numMajorTicks) {
            my $yMinor = $y;
            my $xMinor = $xLeft - $self->{minorTickLen};
            for my $j (1 .. $numMinorTicks-1) {
                $yMinor -= $minorDy;
                $img->line($xLeft, $yMinor, $xMinor, $yMinor, $lineColor);
            }
        }
    }

    $img->line($xLeft, $yBottom, $xRight, $yBottom, $lineColor);

    # Draw y-axis label
    my $yLabelStr = "Number of Samples";

    my $lblLen = length($yLabelStr);

    my $yLabel = $yBottom - ($graphHeight - $mbFontWidth*$lblLen)/2;
    $img->stringUp(gdMediumBoldFont, $xString-15, $yLabel, $yLabelStr, $lineColor);

    # Plot x-axis
    my $xMin = 0;
    my $xMax  = 100;
    my $xStep = 10;
    $img->line($xLeft, $yBottom, $xRight, $yBottom, $lineColor);
    $img->line($xLeft, $yTop, $xRight, $yTop, $lineColor);

    my $numXticks = int(($xMax - $xMin) / $xStep + 0.5);

    my $minorDx = $graphWidth/($numXticks*$self->{xMinorTicks});
    for my $i (0 .. $numXticks) {
        my $x = $xLeft + $i * $graphWidth / $numXticks;
        $img->line($x, $yBottom, $x, $yBottom+$self->{majorTickLen}, $lineColor);

        # Plot minor ticks
        if ($i < $numXticks) {
            my $xMinor = $x;
            for my $j (1 .. $self->{xMinorTicks}-1) {
                $xMinor += $minorDx;
                $img->line($xMinor, $yBottom, $xMinor, $yBottom+$self->{minorTickLen}, $lineColor);
            }
        }

        my $xVal = sprintf("%d", $xMin + $i * $xStep);
        $xVal = sprintf("%2.1f", $xMin + $i * $xStep) if ($xMax < 5);
        $tickLblLen = length($xVal);
        my $xLblPos = $x - $mbFontWidth * $tickLblLen / 2;
        my $yLblPos = $yBottom + 12;
        $img->string(gdMediumBoldFont, $xLblPos, $yLblPos, $xVal, $lineColor);
    }

    my $xAxis = "Percentage of heterozygous genotypes";

    my $strLetters = length($xAxis);
    my $strLen = $mbFontWidth * $strLetters;
    my $xPos = $xLeft + ($graphWidth - $strLen) / 2;
    my $yPosXlbl = $yBottom + 40;
    $img->string(gdMediumBoldFont, $xPos, $yPosXlbl, $xAxis, $lineColor);

    return ($fullVal, $yPosXlbl);
}

sub PlotLegends
{
    my ($self, $img, $xLeft, $yTop, $lineColor, $maleColor, $femaleColor, $otherColor, $showOther) = @_;

    my $lgdGap1 = 10;
    my $lgdGap2 = 50;
    my $lgdGap3 = 65;

    my $lgdLbl = "Phenotype sex:";
    $img->string(gdMediumBoldFont, $xLeft, $yTop, $lgdLbl, $lineColor);

    my $lgdSide = 10;
    my $lgdLen = (length($lgdLbl)+5) * gdMediumBoldFont->width;
    my $xPos = $xLeft + $lgdLen;
    $img->filledRectangle($xPos, $yTop, $xPos + $lgdSide, $yTop+$lgdSide, $maleColor);
    $xPos += $lgdSide + $lgdGap1;
    $img->string(gdMediumBoldFont, $xPos, $yTop, "Male", $lineColor);

    $xPos += $lgdGap2;
    $img->filledRectangle($xPos, $yTop, $xPos + $lgdSide, $yTop+$lgdSide, $femaleColor);
    $xPos += $lgdSide + $lgdGap1;
    $img->string(gdMediumBoldFont, $xPos, $yTop, "Female", $lineColor);

    if ($showOther) {
        $xPos += $lgdGap3;
        $img->filledRectangle($xPos, $yTop, $xPos + $lgdSide, $yTop+$lgdSide, $otherColor);
        $xPos += $lgdSide + $lgdGap1;
        $img->string(gdMediumBoldFont, $xPos, $yTop, "No info", $lineColor);
    }
}

1;
