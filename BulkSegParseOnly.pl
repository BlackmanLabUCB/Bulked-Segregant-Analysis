#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;

# Filter the variant pileup file yielded from the merged bam file and output a list of alleles by position
my $varpileupfile = shift;
my $file_count = 2; 
parse_pileup($varpileupfile, $file_count);

exit;

sub parse_pileup {
	my $pileup_file = shift;
	my $number_of_files = shift;
	my $upper_limit = $number_of_files * 100;
	
	#set up the I/O stream
	open (IN, $pileup_file);
	my $root;
	if ($pileup_file =~ /^(\w*).pileup$/) {$root = $1;}
	my $finalpileup_file = "final.pileup";
	open (OUT, ">$finalpileup_file");
	
	while (<IN>) {
		my $line = $_; chomp $line; $line =~ s/\s/\t/ig;
		my @SNPinfo = split ("\t",$line);
		
		#skip if SNP is an indel
		next if ($SNPinfo[2] eq "\*");
		
		#skip if coverage is too low or too high
		my $coverage = $SNPinfo[3];
		if ($coverage < 4 || $coverage > $upper_limit) {
			next;
		}
		
		$SNPinfo[4] =~ s/[^AGCTagct\.,]//ig;		#remove nonstandard chars

		#skip if there is no SNP
		if ($SNPinfo[4] =~ /^[\,.]$/) {
			next;
		}

		my $consbase = $SNPinfo[2];
		$SNPinfo[4] =~ s/[\.,]/$consbase/ig;		#replace consensus with actual base
		$SNPinfo[4] = uc $SNPinfo[4];				#change all to uppercase
		
		#skip if number of bases given is not equal to read coverage (gets rid of indels assoc line)
		if (length $SNPinfo[4] != $coverage) {
			next;
		}
		#skip if SNP reported because all bases different from reference
		if ($SNPinfo[4] =~ /^A+$|^G+$|^T+$|^C+$/) {
			next;
		}
		
		#parse SNP information into hash keyed by base (A,C,G,T) with read count for a base as value
		my @SNPlist = split ("",$SNPinfo[4]);
		my @QUALlist = split ("",$SNPinfo[5]);
		my %SNPs;
		my %QUALs;
		foreach my $base (@SNPlist) {
			$SNPs{$base}++;
			if (not defined $QUALs{$base}) {
				$QUALs{$base} = shift(@QUALlist);
			} else {
				$QUALs{$base} .= shift(@QUALlist);
			}
		}

		#skip if only one or more than 2 alleles segregate at the SNP
		next if scalar keys %SNPs != 2;
		
		#skip if either allele is supported by fewer than two reads of good quality
		#collect alleles
		my $skip = 0;
		my @alleles;
		$alleles[0] = $consbase;
		foreach my $allele (keys %SNPs) {
			$skip = 1 if $SNPs{$allele} == 1;
			my $good_qual = 0;
			print "$SNPinfo[1] $allele $QUALs{$allele}\n";
			my @allele_quals = split ("",$QUALs{$allele});
			foreach my $pos (@allele_quals) {
				$good_qual++ if $pos gt '3';
			}
			$skip = 1 if $good_qual < 2;
			print "not enough high quality SNPs for allele $allele at $SNPinfo[1]\n" if $good_qual<2;
			next if $skip == 1;
			if ($allele ne $consbase && not defined $alleles[1]) {
				$alleles[1] = $allele;
			} elsif ($allele ne $consbase) {
				$alleles[2] = $allele;
			}
		}
		next if $skip == 1;
		
		#print out SNP if it passed all those criteria
		print OUT "$SNPinfo[0]\t$SNPinfo[1]\t" . join ("\t", @alleles) . "\n";
	}
	close IN;
	close OUT;
}
