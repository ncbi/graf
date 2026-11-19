use strict;
use warnings FATAL => 'all';
use Test::More;
use SexCheckResults;

my $INPUT_FILE = 'data/SexCheckResTest01.txt';

subtest 'checking creation of SexCheckResults ' , sub {
    my $sexCheckRes = new SexCheckResults();
    my $fileErr = $sexCheckRes->GetSampleSexResults($INPUT_FILE);

    is($fileErr,"");
    is($sexCheckRes->{numAllSmps}, 552);
    is($sexCheckRes->{numFpSmps}, 552);
    is($sexCheckRes->{numMales}, 328);
    is($sexCheckRes->{numFemales}, 224);
    is($sexCheckRes->{numUnkSexes}, 0);
    is($sexCheckRes->{meanXsnps}, 146);
    is($sexCheckRes->{maxBinSmps}, 327);

};

done_testing;