### macs2 receive bigwig and then do the correlation  
### merge all the correlations from wigCorrelation
### Shicheng Guo, Shihcheng.Guo@Gmail.com
### 12/30/2019

use strict;
my @file=glob("*.o*");
my %data;
foreach my $file(@file){
if($file =~ m/(\d*_\w*).bw\.+(\d*_\w*).bw/){
open F,$file;
while(<F>){
my($cor)=split /\s+/;
$data{$1}{$2}=$cor;
}
}
}

my @sam=sort keys %data;
my $header=join("\t",@sam);
print "\t$header\n";
foreach my $sam1(sort keys %data){
    print "$sam1";
        foreach my $sam2(sort keys %data){
        print "\t$data{$sam1}{$sam2}";
        }
        print "\n";
}
