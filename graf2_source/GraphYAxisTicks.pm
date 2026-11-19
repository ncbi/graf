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
#####################################################################################################
#
# The purpose of this package is to select numbers of major and minor ticks and tick intervals
# so that they can be used to plot data on the y-axis by default.

# The following criteria are used when selecting these values:

# 1. Display all values, leaving some space on top (Don't let the plotted value touch the top edge)
# 2. Do not leave too much space, making sure no more than ~2/3 vertical spaces is not empty
# 3. Use reasonable tick intervals like 5, 10, 1000, etc.
	# 4  Use reasonable number of major and minor ticks
#
#####################################################################################################

package GraphYAxisTicks;

use strict;
use warnings;

sub new
{
    my ($class) = @_;

    # Settig of ticks and tick intervals based on the value to plot

    # Column 1: maximum y-value to plot
    # Column 2: value represented by each major tick
    # Column 3: number of major ticks
    # Column 4: number of minor ticks with each major tick
    my @valueTickStrs = (
    " 5,  1,  5,  5",
    " 8,  1,  8,  2",
    "10,  1, 10,  2",
    "15,  5,  3,  5",
    "20,  5,  4,  5",
    "25,  5,  5,  5",
    "30,  5,  6,  5",
    "40,  5,  8,  5",
    "50, 10,  5, 10" );

    my %valueTicks = ();
    for my $tickStr (@valueTickStrs) {
        my @rowVals = split /\,/, $tickStr;
        my $yValue     = $rowVals[0] + 0;
        my $tickIntv   = $rowVals[1] + 0;
        my $majorTicks = $rowVals[2] + 0;
        my $minorTicks = $rowVals[3] + 0;

        my %ticks = (tickIntv => $tickIntv, majorTicks => $majorTicks, minorTicks => $minorTicks);

        $valueTicks{$yValue} = \%ticks;
    }

    my $self = {
        maxFirstDigit => 49,
        valueTicks => \%valueTicks
    };

    bless ($self, $class);

    return $self;
}

sub SetValueTicks
{
    my ($self, $val, $intv, $numMajors, $numMinors) = @_;

    $self->{$val}->{tickIntv} = $intv;
    $self->{$val}->{majorTicks} = $numMajors;
    $self->{$val}->{minorTicks} = $numMinors;
}

sub ShowValueTicks
{
    my $self = shift;

    print "\nTick settings for different maximum y values to be plotted\n\n";

    my $line = sprintf("%7s  %8s  %10s  %10s", "Value", "Interval", "Major ticks", "Minor ticks");
    print "$line\n";

    my %valTicks = %{$self->{valueTicks}};

    my $minVal = 1;
    foreach my $value (sort {$a <=> $b} keys %valTicks) {
        my %tickInfo = %{$valTicks{$value}};
        my $intv = $tickInfo{tickIntv};
        my $majorTicks = $tickInfo{majorTicks};
        my $minorTicks = $tickInfo{minorTicks};

        my $maxVal = $value - 1;

        $line = sprintf("%2d - %2d  %8d  %10d  %10d", $minVal, $maxVal, $intv, $majorTicks, $minorTicks);
        print "$line\n";

        $minVal = $value;
    }
    print "\n";
}


#
# Given an integer, determine the most significant digit (one or two digits), and number of digits
#
# For example, 21573: set $sigDigit = 21, $unit = 1000, $numDigits = 5
#
sub GetFirstDigit
{
    use integer;
    my ($self, $yVal) = @_;

    my $sigDigit = $yVal;

    die "\nERROR: Number '$yVal' is not a positive integer!\n" if ($yVal < 1);

    my $unit = 1;
    while ($sigDigit > $self->{maxFirstDigit}) {
        $sigDigit = $sigDigit / 10;
        $unit = $unit * 10;
    }

    $sigDigit++ if ($yVal > $sigDigit * $unit);

    return ($sigDigit, $unit);
}

#
# Given significant digit from GetFirstDigit, set number of ticks and tick intervals
# for plotting data on a graph
#
sub GetTicks
{
    my ($self, $yVal, $isManual) = @_;

    die "\nERROR: Number '$yVal' is smaller than 1!\n" if ($yVal < 1);

    my ($firstDigit, $unit) = $self->GetFirstDigit($yVal);

    my $majorTicks = 0;
    my $minorTicks = 0;
    my $interval = 0;

    my %valTicks = %{$self->{valueTicks}};

    foreach my $value (sort {$a <=> $b} keys %valTicks) {
        my %tickInfo = %{$valTicks{$value}};

        # If the scale is manually set by the user, leave no space on top.
	        # For example, if the max y-value is 20, by default, the script sets scale to 25,
        # leaving some space on top. However, if the user specifies the scale as 20,
        # then the script leaves no space on top.
	        if ($firstDigit < $value || ($isManual && $firstDigit == $value)) {
            $interval   = $tickInfo{tickIntv} * $unit;
            $majorTicks = $tickInfo{majorTicks};
            $minorTicks = $tickInfo{minorTicks};

            last;
        }
    }

    unless ($interval) {
        die "\nERROR: failed to determine ticks for first digit $firstDigit!\n";
    }

    return ($interval, $majorTicks, $minorTicks);
}

1;

