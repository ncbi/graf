#!/opt/perl/5.26.3/bin/perl
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
        File Description: script to convert SNP IDs to RS IDs in bim file, based on
                          chromosome and position values
        Date: 04/17/2020
    ====================================================================================
EOF

use strict;
use Getopt::Long;
use Cwd 'abs_path';

if (@ARGV < 1 || $ARGV[0] =~ /^-/) {
    print "\n$disclaim\n\n";
    print "This script checks chromosome and position for each SNPs in the PLINK .bim file to determine\n";
    print "if it is a GRAF fingerprint SNPs. Assuming genome build GRCh 37 or 38, the script does the following:\n\n";
    print "    if the SNP is a GRAF fingerprint SNP then\n";
    print "        if the SNP ID is an RS ID then\n";
    print "            if the SNP ID is different from that expected based on chromosome and position then\n";
    print "                replace the SNP ID as \"RS_ID_conflict\"\n";
    print "                report it in the WARNING message\n";
    print "        else\n";
    print "            replace the SNP ID with rs ID\n\n";
    print "If output file is specified, it saves the new results to the output file. Otherwise it updates";
    print "the input file and saves the original .bim file to <plink_set_basename>_orig.bim\n\n";
    print "Usage: SetRsIdsInBimFile.pl <plink_bim_file> [output_bim_file]\n";
    print "\n";
    exit 0;
}

my %opt = ();

my $ret = GetOptions(\%opt, "fpSnpFile|fpsnps=s");

#------------------------------------------- Global constants -----------------------------------------------#
my $NUM_ALL_FP_SNPS = 11000;
my $MIN_FP_SNPS_FOR_BUILD = 100; # number of FP SNPs needed to determine genome build
my $NUM_BIM_COLUMNS = 6;


#------------------------------------------- Check input file -----------------------------------------------#
my $bimFile = $ARGV[0];
my $outputFile = $ARGV[1];
$outputFile = "" if ($outputFile =~ /^-/);

my $script = $0;
my $pathToScript = abs_path($script);
my $dataDir  = ".";
if ($pathToScript =~ /(.+)\//) {
    $dataDir = "$1/data";
}

my ($bimErr, $origFile) = QcCheckBimFile($bimFile, $outputFile);
if ($bimErr) {
    print "\n$bimErr!\n\n";
    exit 1;
}

#------------------------------------------- Read FP SNPs from file  ----------------------------------------#
# Global variables to avoid passing long hash tables
my %gb37FpChrPos = ();
my %gb38FpChrPos = ();

my $fpSnpFile = "$dataDir/FP_SNPs.txt";
$fpSnpFile = $opt{fpSnpFile} if ($opt{fpSnpFile});

die "\nERROR: didn't find FP SNP file $fpSnpFile!\n" unless (-e $fpSnpFile);

GetFpSnpChromosomePositions($fpSnpFile);

#------------------------------------------- Determine genome build -----------------------------------------#
my $genoBld = GetGenomeBuild($bimFile);
unless ($genoBld) {
    print "\nERROR: bim file might be either not enough FP SNPs for GRAF, or in a build other than GB37 or GB38.\n";
    exit 2;
}


#------------------------------------------- Set RS IDs -----------------------------------------------------#
SetRsIdsForBimFile($bimFile, $genoBld, $origFile);


#------------------------------------------- Functions ------------------------------------------------------#
sub GetFpSnpChromosomePositions
{
    my $fpSnpFile = shift;

    my %g37Pos = ();
    my %g38Pos = ();

    open FILE, $fpSnpFile or die "Couldn't open $fpSnpFile!\n";
    my $header = <FILE>;
    if ($header !~ /^rs\#\tchromosome\tGB37_position\tGB38_position/) {
	die "\nERROR: $fpSnpFile is not the GRAF FP SNP file!\nheader $header\n";
    }

    my $numSnps = 0;
    while (<FILE>) {
	chomp;
	my ($rs, $chr, $g37, $g38, $a1, $a2) = split /\t/, $_;
	if ($rs && $chr && $g37 && $g38) {
	    $gb37FpChrPos{"$chr\t$g37"} = $rs;
	    $gb38FpChrPos{"$chr\t$g38"} = $rs;
	}

	$numSnps++;
    }

    close FILE;

    die "\nERROR: found $numSnps FP SNPs in $fpSnpFile. Expected $NUM_ALL_FP_SNPS FP SNPs.\n" if ($numSnps != $NUM_ALL_FP_SNPS);

    print "\nFound $numSnps FP SNPs in $fpSnpFile.\n";
}

sub QcCheckBimFile
{
    my ($bimFile, $outputFile) = @_;

    my ($error, $bimDir, $origFile) = ("", "", "");

    if ($bimFile !~ /\.bim$/i) {
	$error = "File $bimFile doesn't have extension .bim";
    }
    elsif (not -e $bimFile) {
	$error = "File $bimFile doesn't exist";
    }
    elsif (not -T $bimFile) {
	$error = "File $bimFile is not a plain text file";
    }
    elsif (!$outputFile) {
	$bimDir = $bimFile;
	if ($bimFile =~ /(.+)\//) {
	    $bimDir = $1;
	}
	else {
	    $bimDir = ".";
	}

	$origFile = $bimFile;
	$origFile =~ s/\.bim$/_orig\.bim/;

	if (not -w $bimDir) {
	    $error = "Directory $bimDir/ is not writeable";
	}
	elsif (-e $origFile) {
	    $error = "File $origFile already exists. Please delete or remove this file from directory '$bimDir'";
	}

	my $headLine = `head -1 $bimFile`;
	my @columns = split /\s+/, $headLine;
	my $numColumns = @columns;
	if ($numColumns != $NUM_BIM_COLUMNS) {
	    $error = "$bimFile is invalid. The first row has $numColumns columns (expected $NUM_BIM_COLUMNS columns)";
	}
    }

    return ($error, $origFile);
}

sub GetGenomeBuild
{
    my $bimFile = shift;

    print "Checking genome build of $bimFile.\n";

    my $bldNo = 0;

    my ($numChkSnps, $numGb37Fps, $numGb38Fps) = (0, 0, 0);

    open FILE, $bimFile or die "\nERROR: Couldn't open file $bimFile!\n";
    while(<FILE>) {
	chomp;

	my ($chr, $snp, $pos) = GetChromosomePosition($_);

	if ($gb37FpChrPos{"$chr\t$pos"}) {
	    $numGb37Fps++;
	}

	if ($gb38FpChrPos{"$chr\t$pos"}) {
	    $numGb38Fps++;
	}

	if ($numGb37Fps >= $MIN_FP_SNPS_FOR_BUILD || $numGb38Fps >= $MIN_FP_SNPS_FOR_BUILD) {
	    if ($numGb37Fps > $numGb38Fps) {
		$bldNo = 37;
	    }
	    else {
		$bldNo = 38;
	    }

	    last;
	}

	$numChkSnps++;

	if ($numChkSnps % 100000 == 0) {
	    print "\tChecked $numChkSnps SNPs.\n";
	}
    }
    close FILE;

    print "Checked total $numChkSnps SNPs. File is with GRCh$bldNo\n";

    return $bldNo;
}

sub SetRsIdsForBimFile
{
    my ($bimFile, $gbNum, $origFile) = @_;

    print "\nSetting RS IDs for file $bimFile.\n";
    my @newFileLines = ();

    my ($numSnps, $numFpSnps, $numSetSnps, $numErrSnps) = (0, 0, 0, 0);

    open FILE, $bimFile or die "\nERROR: Couldn't open file $bimFile!\n";
    while(<FILE>) {
	chomp;
	my $line = $_;

	my ($chr, $snp, $pos) = GetChromosomePosition($line);

	my $fpRsNum = 0;
	if ($gbNum == 37) {
	    $fpRsNum = $gb37FpChrPos{"$chr\t$pos"};
	}
	elsif ($gbNum == 38) {
	    $fpRsNum = $gb37FpChrPos{"$chr\t$pos"};
	}

	if ($fpRsNum) {
	    if ($snp =~ /^rs(\d+)/i) {
		my $fileRs = $1;

		if ($fileRs == $fpRsNum) {
		    push @newFileLines, $line;
		}
		else {
		    $line =~ s/$snp/RS_ID_conflict/;
		    push @newFileLines, $line;
		    $numErrSnps++;
		}
	    }
	    else {
		$line =~ s/$snp/rs$fpRsNum/;
		push @newFileLines, $line;
		$numSetSnps++;
	    }

	    $numFpSnps++;
	}
	else {
	    push @newFileLines, $line;
	}

	$numSnps++;

	if ($numSnps % 100000 == 0) {
	    print "\tChecked $numSnps SNPs. Found $numFpSnps FP SNPs. Set RS IDs for $numSetSnps SNPs.";
	    print " $numErrSnps have RS ID conflicts." if ($numErrSnps > 0);
	    print "\n";
	}
    }
    close FILE;

    print "Checked total $numSnps SNPs. Found $numFpSnps FP SNPs. Set RS IDs for $numSetSnps SNPs.";
    print " $numErrSnps have RS ID conflicts." if ($numErrSnps > 0);
    print "\n";

    if ($numSetSnps + $numErrSnps > 0) {
        unless($outputFile) {
	    system("mv $bimFile $origFile");
	    die "\nERROR: Failed to rename $bimFile to $origFile!\n" if (-e $bimFile || not -e $origFile);
	    $outputFile = $bimFile;
	}

	open FILE, ">$outputFile" or die "\nERROR: couldn't open $outputFile for writing!\n";
	for my $line (@newFileLines) {
	    print FILE "$line\n";
	}
	close FILE;

	print "Done. Set RS IDs for $numSetSnps of $numSnps SNPs.\n";
	if ($outputFile) {
	    print "Results saved to $outputFile.\n\n";;
	}
	else {
	    print "Updated $bimFile.\nOriginal file saved to $origFile.\n\n";;
	}

	if ($numErrSnps > 0) {
	    print "\nWARNING: $numErrSnps FP SNPs have RS ID conflicts between SNP IDs and " .
	          "the RS IDs determined from chromosomes and positions.\n";
	    print "These SNPs are marked as 'RS_ID_conflicts' in the bim file and will be ignored by GRAF.\n";
	}
    }
    else {
	print "\nFile $bimFile already includes SNPs with RS IDs.\n";
    }
    print "\n";

    return ($numSnps, $numFpSnps, $numSetSnps, $numErrSnps);
}


sub GetChromosomePosition
{
    my $line = shift;

    my @cols = split /\s+/, $line;
    my $chrStr = $cols[0];
    my $snp = $cols[1];
    my $pos = $cols[3];

    my $chr = $chrStr;
    $chr = 23 if ($chrStr =~ /^x/i || $chrStr =~ /^chrx/i);

    return ($chr, $snp, $pos);
}
