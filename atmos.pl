#!/usr/bin/perl -w

my $fn="atmosph.in";
open(INP, "<", $fn) or die "can not open $fn\n";

# read head lines
my $line;
for(my $i=1;$i<=12;$i++) {
    $line=<INP>;
    print "$line";
}

# read data
my @col;
while(<INP>)
{
    last if /^end/;
    @col=split;
    $col[8]=($col[12]+$col[13])/2.0;
    print join(" ",@col),"\n";
}

# print last 2 lines
print;
$line=<INP>;
print "$line";
