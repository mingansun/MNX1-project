#!/usr/bin/perl
use warnings;
use strict;

=head1 Description

This script is used to get the flanking phastCons profile surrounding specified regions.
1) The regions are specified in BED file. 2) in phastCons dir, the data are for each chr
and named by chr name; 3) flanking length is the flanking length for one side. For example,
5000 means 5000 bp upstream and 5000 bp downstream.


=head1 Usage

perl get_flanking_phastCons.pl <bed file> <phastcons dir> <flanking length>  >  output

=head1 Example

perl get_flanking_phastCons.pl hg19.CSM.bed phastCons_dir 5000 >hg19.CSM.bed.phastCons

=head1 Version

v1.0, Ming-an Sun, May 20, 2014
v1.1, Ming-an Sun, Aug 1, 2021

=cut

die `pod2text $0` unless @ARGV==3;

my ($bed_file, $phastcons_dir, $flank_length) = @ARGV;

$phastcons_dir =~ s/\/$//;

warn "\nReading coordinates from $bed_file ...\n";
my %regions;
my @coords;
my $con = 0;
open(BED, $bed_file)||die"Cannot open $bed_file\n";
while(<BED>){
    $con++;
    print "$con\n" if $con%1000000==0;
    my ($chr, $start, $end, $strand) = split;
    push(@{$regions{$chr}}, [$start, $end, $strand]);
    push(@coords, "$chr\t$start\t$end\t$strand");
}
close BED;

warn "\nCheck the phastCons flanking specified regions ...\n";
my %stat;
my %info;
foreach my $chr (sort keys %regions){
    warn "$chr\t";

    # coordnates for which the phastcons will be extracted
    my %coord;
    foreach my $region (@{$regions{$chr}}){
        my $start  = $region->[0];
        my $end    = $region->[1];
        my $strand = $region->[2];
        my $middle = int(($start+$end)/2);
        my $flank_start = $middle-$flank_length;
        my $flank_end   = $middle+$flank_length;
        for(my $i=$flank_start; $i<=$flank_end; $i++){
            $coord{$chr}->{$i} = '';
        }
    }
    
    # read phastcons scores
    my $phastcons_file;
    if(-e "$phastcons_dir/$chr"){
        $phastcons_file = "$phastcons_dir/$chr";
    }elsif(-e "$phastcons_dir/$chr\.gz"){
        $phastcons_file = "$phastcons_dir/$chr\.gz";
    }else{
        warn "Phastcons file for $chr undetected in $phastcons_dir!\n";
        next;
    }
    my %phastcons;
    &read_phastcon_to_hash(\%phastcons, $phastcons_file, \%coord);
    warn "read2hash\t";

    foreach my $region (@{$regions{$chr}}){
        my @tmp;
        my $start  = $region->[0];
        my $end    = $region->[1];
        my $strand = $region->[2];
        my $middle = int(($region->[0]+$region->[1])/2);
        for(my $idx=0-$flank_length; $idx<=$flank_length; $idx++){
            my $i = $strand eq "-" ? $middle-$idx : $middle+$idx;
            if(defined $phastcons{$chr}->{$i}){
                $stat{$idx}->[0]++;
                $stat{$idx}->[1] += $phastcons{$chr}->{$i};
                push(@tmp, $phastcons{$chr}->{$i});
            }
            else{
                push(@tmp, 'NA');
            }
        }
        $info{"$chr\t$start\t$end\t$strand"} = join("\t", @tmp);

    }
    delete($regions{$chr});
    delete($phastcons{$chr});
    warn "OK\n";
}

warn "\nOutput results\n";

open(OUT1, ">$bed_file\.flankCons")||die"Cannot open $bed_file\.flankCons\n";
foreach my $x (@coords){
    print OUT1 "$x\t$info{$x}\n";
}

close OUT1;

open(OUT2, ">$bed_file\.flankAvgCons")||die"Cannot open $bed_file\.flankAvgCons\n";
for(my $idx=0-$flank_length; $idx<=$flank_length; $idx++){
    print "$idx\t";
    print $stat{$idx}->[0] . "\t";
    print $stat{$idx}->[1] . "\t";
    print OUT2 sprintf("%.5f", $stat{$idx}->[1]/$stat{$idx}->[0]);
	print OUT2 "\n";
}
close OUT2;

warn "\nDone.\n";


############## subroutines ################

## read phastCons score to hash
sub read_phastcon_to_hash{
    my ($hashP, $phastcons_file, $coordHashP) = @_;
    open(IN, $phastcons_file =~ /\.gz$/ ? "zcat $phastcons_file |" : $phastcons_file)||die"Cannot open $phastcons_file\n";
    $/ = "fixedStep";
    while(my $ln = <IN>){
        my @a = split("\n", $ln);
        next unless @a>=2;
        pop @a if $a[-1] =~ /fixedStep/; # remove the last line if necessary
        my $header = shift @a;
        my ($this_chr, $start);
        if($header =~ /chrom=(chr\S+)\s+start=(\d+)/){
            ($this_chr, $start) = ($1, $2);
        }else{
            die "ERROR: cannot parse $header in $phastcons_file\n";
        }
        my $position = $start;
        foreach(@a){
            if(/^([\d\.]+)/){
                if(defined ${$coordHashP}{$this_chr}->{$position}){
                    ${$hashP}{$this_chr}->{$position} = $1;
                    $position++;
                }
            }
        }
    }
    $/ = "\n";
    close IN;
}


__END__

# bed file
chr10	73522	73579	+
chr10	73536	73591	+
chr10	73555	73601	+
chr10	73630	73662	+
chr10	74103	74153	+
chr10	81764	81802	+
# phastcons file of wig format
fixedStep chrom=chr1 start=10918 step=1
0.254
0.253
0.251
0.249
0.247
0.244
0.242
0.239
0.236
