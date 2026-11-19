#!/usr/bin/perl5.16

use strict;
use GraphYAxisTicks;
use Test::More;


my @testValues = (
"     0,      0,  0,  0",
"     1,      1,  5,  5",
"     2,      1,  5,  5",
"     5,      1,  8,  2",
"     6,      1,  8,  2",
"    10,      5,  3,  5",
"    12,      5,  3,  5",
"    15,      5,  4,  5",
"    24,      5,  5,  5",
"    32,      5,  8,  5",
"    47,     10,  5, 10",
"    50,     10,  8,  2",
"    72,     10, 10,  2",
"   100,     50,  3,  5",
"   101,     50,  3,  5",
"  1325,    500,  3,  5",
"  1800,    500,  4,  5",
" 61576,  10000,  8,  2",
"386471,  50000,  8,  5",
"520487, 100000,  8,  2",
"  -100,      0,  5,  0",
);

sub Test_GetTicks_WithValue
{
    my ($axisTicks, $testVal, $expIntv, $expMajor, $expMinor) = @_;
    my ($intv, $numMajorTicks, $numMinorTicks) = ();

    eval { ($intv, $numMajorTicks, $numMinorTicks) = $axisTicks->GetTicks($testVal); };
    my $hasErr = $@ ? "Invalid y value" : "Y value OK";

    if ($expIntv) {
	is($intv, $expIntv, "correct interval");
	is($numMajorTicks, $expMajor, "correct major ticks");
	is($numMinorTicks, $expMinor, "correct minor ticks");
    }
    else {
	is("Invalid y value", $hasErr, "correct error report for invalid y value $testVal");
    }
}

for my $valStr (@testValues) {
	my $axisTicks = new GraphYAxisTicks();
    # $axisTicks->ShowValueTicks();
	my ($testVal, $expIntv, $expMajor, $expMinor) = split /\,/, $valStr;
	$testVal += 0;
	$expIntv += 0;
	$expMajor += 0;
	$expMinor += 0;

	my @args = ($axisTicks, $testVal, $expIntv, $expMajor, $expMinor);
	subtest "testing $valStr" => sub { Test_GetTicks_WithValue(@args)};
}

done_testing;

