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
package SexCheckResults;

use strict;
use warnings;

my $SEX_FILE_COLUMNS = 7;

sub new
{
    my ($class) = @_;

    my $self = {
        numAllSmps => 0,
        numFpSmps  => 0,
        binSize => 1,

        minFpSnpsForSexCheck => 20,

        # The cutoff values of heterozygosity rates for males and females
        maxMaleHetroRate   => 8,
        minFemaleHetroRate => 20,
        maxFemaleHetroRate => 60,

        smpSexStats => {},  # sex check results such as pheno/geno sexes, genotype counts for each sample
        summary => [],      # summary info to be saved to output file
    };

    bless ($self, $class);

    return $self;
}

sub SetParameters
{
    my ($self, $params) = @_;

    foreach my $key (keys %$params) {
        if ($key && defined $self->{$key}) {
            $self->{$key} = $params->{$key};
        }
    }

    if ($self->{maxMaleHetroRate} > 20) {
        die "\nERROR: Invalid maximum male heterozygosity rate (%) $self->{maxMaleHetroRate}: expected 0 - 20\n"
    }

    if ($self->{minFemaleHetroRate} < $self->{maxMaleHetroRate} || $self->{maxFemaleHetroRate} > 100) {
        die "\nERROR: Invalid minimum female heterozygosity rate (%) $self->{minFemaleHetroRate}: " .
	        "expected $self->{maxMaleHetroRate} - 100\n";
    }

    if ($self->{maxFemaleHetroRate} < $self->{minFemaleHetroRate} || $self->{maxFemaleHetroRate} > 100) {
        die "\nERROR: Invalid maximum female heterozygosity rate (%) $self->{maxFemaleHetroRate}: " .
	        "expected $self->{minFemaleHetroRate} - 100\n";
    }
}

sub SummarizeResults
{
    my ($self) = @_;

    my $numBins = int(100 / $self->{binSize});

    # Counts of samples with heterozygotic genotypes
    my @g12Smps = ();
    my @g12Males = ();
    my @g12Females = ();
    for my $binNo (0 .. $numBins) {
        push @g12Smps, 0;
        push @g12Males, 0;
        push @g12Females, 0;
    }

    my ($minXsnps, $maxXsnps) = (10000, 0);
    my ($numAllSmps, $numFpSmps, $numMales, $numFemales) = (0, 0, 0, 0);
    my %smpSexStats = %{$self->{smpSexStats}};
    my $totSnps = 0;

    my @errLines = ();       # Samples with pheno and geno sex conflicts
    my @errGenoLines = ();   # Samples with potential errors in genotpes
    my $numValidSexSmps = 0; # Number of samples with pheno and geno sexes both valid
    my $numMatchSexSmps = 0; # Number of samples with pheno and geno sex matches

    foreach my $smpName (keys %smpSexStats) {
        next unless $smpSexStats{$smpName};

        my ($smp, $phenoSex, $rptGenoSex, $sexMatch, $numHetros, $numFpSnps, $hetroRate) = @{$smpSexStats{$smpName}};

        my $saveLine = "$smp\t$phenoSex\t$rptGenoSex\t$numHetros\t$numFpSnps\t$hetroRate";

        my $genoSex = 0;

        if ($numFpSnps >= $self->{minFpSnpsForSexCheck}) {
            $minXsnps = $numFpSnps if ($numFpSnps < $minXsnps);
            $maxXsnps = $numFpSnps if ($numFpSnps > $maxXsnps);
            $totSnps += $numFpSnps;

            my $g12Rate = -1;
            $g12Rate = $numHetros * 100.0 / $numFpSnps;

            my $binNo = int($g12Rate / $self->{binSize});
            $g12Smps[$binNo]++;

            if ($g12Rate < $self->{maxMaleHetroRate}) {
                $genoSex = 1;
            }
            elsif ($g12Rate >= $self->{minFemaleHetroRate} && $g12Rate <= $self->{maxFemaleHetroRate}) {
                $genoSex = 2;
            }
            elsif ($g12Rate > $self->{maxFemaleHetroRate}) {
                push @errGenoLines, $saveLine;
            }

            $saveLine = "$smp\t$phenoSex\t$genoSex\t$numHetros\t$numFpSnps\t$hetroRate";
            if ($phenoSex) {
                if ($phenoSex == 1) {
                    $g12Males[$binNo]++;
                    $numMales++;

                    if ($genoSex == 2) {
                        push @errLines, $saveLine;
                    }
                    elsif ($genoSex == 1) {
                        $numMatchSexSmps++;
                    }

                    $numValidSexSmps++;
                }
                elsif ($phenoSex == 2) {
                    $g12Females[$binNo]++;
                    $numFemales++;

                    if ($genoSex == 1) {
                        push @errLines, $saveLine;
                    }
                    elsif ($genoSex == 2) {
                        $numMatchSexSmps++;
                    }

                    $numValidSexSmps++;
                }
            }

            $numFpSmps++;
        }

        $numAllSmps++;
    }

    my $numUnkSexes = $numFpSmps - $numMales - $numFemales;

    my $meanXsnps = $numFpSmps > 0 ? int($totSnps / $numFpSmps + 0.5) : 0;

    my $maxBinSmps = 0;
    foreach my $binNo (0 .. $numBins) {
        $maxBinSmps = $g12Smps[$binNo] if ($maxBinSmps < $g12Smps[$binNo]);
    }

    my $numErrGenoSmps = @errGenoLines;
    my $numErrSmps = @errLines;

    $self->{g12Smps}    = \@g12Smps;
    $self->{g12Males}   = \@g12Males;
    $self->{g12Females} = \@g12Females;

    $self->{maxBinSmps} = $maxBinSmps;

    $self->{numAllSmps} = $numAllSmps;
    $self->{numFpSmps}  = $numFpSmps;

    $self->{minXsnps}   = $minXsnps;
    $self->{maxXsnps}   = $maxXsnps;
    $self->{meanXsnps}  = $meanXsnps;

    $self->{numMales}    = $numMales;
    $self->{numFemales}  = $numFemales;
    $self->{numUnkSexes} = $numUnkSexes;

    my $colHeader = "Sample\tphenotype_sex\tgenotype_sex\t#heterozygous_SNPs\t#compared_SNPs\theterozygous_rate";

    my @header = ();
    push @header, "Total $numAllSmps samples found. $numFpSmps samples have enough genotypes to determine sexes.";
    push @header, "Phenotype sex counts: $numMales males, $numFemales females, $numUnkSexes unknown sex samples";
    push @header, "Mean FP SNPs/sample: $meanXsnps";
    push @header, "$numValidSexSmps samples have known phenotype sexes and enough genotypes to determine sexes.";
    push @header, "$numMatchSexSmps samples have phenotype sexes confirmed with genotypes.";
    push @header, "\n$numErrSmps samples have phenotype and genotype sex discrepancies.";

    my @outputLines = ();
    for my $line (@header) {
        print "$line\n";
        push @outputLines, $line;
    }

    if ($numErrSmps > 0) {
        push @outputLines, "\n$colHeader";
        for my $line (@errLines) {
            push @outputLines, $line;
        }
    }

    if ($numErrGenoSmps > 0) {
        my $headLine = "\nWARNING: $numErrGenoSmps samples have unexpected heterozygosity rates";
        print "$headLine\n";
        push @outputLines, $headLine;

        push @outputLines, "\n$colHeader";
        for my $line (@errGenoLines) {
            push @outputLines, $line;
        }
    }

    $self->{summary} = \@outputLines;
}

sub SaveSummary
{
    my ($self, $outSumFile) = @_;

    my @outputLines = @{$self->{summary}};

    open FILE, ">$outSumFile" or die "\nERROR: couldn't open $outSumFile for writing.\n\n";
    for my $line (@outputLines) {
        print FILE "$line\n";
    }
    close FILE;

    print "Sex check summary saved to $outSumFile\n\n";
}

sub GetSampleSexResults
{
    my ($self, $sexResFile, $smpSexFile) = @_;

    my $error = "";
    my %smpSexRes = ();

    my $rowNo = 1;
    open FILE, $sexResFile or die "\nERROR: Couldn't open $sexResFile!\n";
    while (<FILE>) {
        chomp;
        next if ($_ !~ /\S/);
        next if ($_ =~ /sample_id\tpheno_sex/);

        my @vals = split /\t/, $_;
        my $numCols = @vals;

        if ($numCols != $SEX_FILE_COLUMNS) {
            $error = "$numCols columns found at Row #$rowNo in $sexResFile. Expected $SEX_FILE_COLUMNS columns.";
            last;
        }

        my $smp = $vals[0];

        if ($smpSexRes{$smp}) {
            $error = "Duplicate sample $smp found in $sexResFile.";
            last;
        }
        else {
            $smpSexRes{$smp} = \@vals;
        }

        # Make sure columns 1 - 5 are all integers
        for my $col (1 .. $#vals-1) {
            my $colNum = $col + 1;
            my $val = $vals[$col];
            $val =~ s/^-//;
            if ($val !~ /^\d+$/) {
                $error = "Non-integer value '$vals[0]' found in $sexResFile, Row #$rowNo column #$colNum.";
                last;
            }
        }

        last if ($error);

        $smpSexRes{$smp} = \@vals;

        $rowNo++;
    }
    close FILE;

    $self->{smpSexStats} = \%smpSexRes;
    $self->{numAllSmps} = keys %smpSexRes;

    print "Read sex check results of $self->{numAllSmps} samples from $sexResFile.\n\n";

    if ($smpSexFile) {
        $self->SetSamplePhenotypeSexes($smpSexFile);
    }

    $self->SummarizeResults();

    return (\%smpSexRes, $error);
}

sub SetSamplePhenotypeSexes
{
    my ($self, $smpSexFile) = @_;

    die "$smpSexFile is not a plain text file.\n" unless (-T $smpSexFile);

    print "Resetting phenotype sexes using data in $smpSexFile\n\n";

    my ($numSmps, $numOverlapSmps, $numExtraSmps, $numInvSexSmps) = (0, 0, 0, 0);
    my ($numNewSexSmps, $numDelSexSmps, $numUpdSexSmps) = (0, 0, 0);

    my %famSmpSexes = ();
    my %smpSexStats = %{$self->{smpSexStats}};
    foreach my $smpName (keys %smpSexStats) {
        my ($smp, $phenoSex, $rptGenoSex, $sexMatch, $numHetros, $numFpSnps, $hetroRate) = @{$smpSexStats{$smpName}};
        $famSmpSexes{$smpName} = $phenoSex;
    }

    my %smpSexes = ();
    open FILE, $smpSexFile or die "\nERROR: couldn't open $smpSexFile.\n";
    while(<FILE>) {
        chomp;
        my ($smpName, $sex) = split /\s+/, $_;
        next unless ($smpName && $sex !~ /^sex/i);

        if ($sex =~ /^(\d)$/ && $sex < 3) {
            if (defined $famSmpSexes{$smpName}) {
                my $famSex = $famSmpSexes{$smpName};

                if ($sex && !$famSex) {
                    $numNewSexSmps++;
                }
                elsif (!$sex && $famSex) {
                    $numDelSexSmps++;
                }
                elsif ($sex && $famSex) {
                    if ($sex != $famSex) {
                        $numUpdSexSmps++;
                    }
                }

                $smpSexes{$smpName} = $sex;
                $numOverlapSmps++;
            }
            else {
                $numExtraSmps++;
            }
        }
        else {
            print "smp $smpName sex '$sex'\n";
            $numInvSexSmps++;
        }

        $numSmps++;
    }
    close FILE;

    print "Read phenotype sexes of $numSmps samples from file.\n";
    print "\t$numOverlapSmps samples are also found in the fam file.\n";
    print "\tWARNING: $numExtraSmps samples are not found in the fam file.\n" if ($numExtraSmps > 0);

    if ($numInvSexSmps > 0) {
        die "\nERROR: $numInvSexSmps samples do not have valid sex values (0, 1 or 2)\n";
    }

    my $numTouchSexSmps = $numNewSexSmps + $numDelSexSmps + $numUpdSexSmps;

    if ($numTouchSexSmps > 0) {
        print "\nIn comparison to the fam file, below are the changes in the sample sex file:\n";
        print "\t$numNewSexSmps samples have sex values changed from unknown to male or female\n" if ($numNewSexSmps > 0);
        print "\t$numDelSexSmps samples have sex values changed from male or female to unknown\n" if ($numDelSexSmps > 0);
        print "\t$numUpdSexSmps samples have sex values changed from male to female or vise versa\n" if ($numUpdSexSmps > 0);
        print "\n";

        foreach my $smpName (keys %smpSexes) {
            if ($self->{smpSexStats}->{$smpName}) {
                $self->{smpSexStats}->{$smpName}->[1] = $smpSexes{$smpName};
            }
        }
    }
    else {
        print "\nSample phenotype sexes of samples have already read from the PLINK fam file.\n\n";
    }
}

1;
