#!/usr/bin/perl -w
use warnings;
use strict;

#receive input data from user
my @filelist = @ARGV or die "no filenames provided\n usage: perl BulkSegSlidingWindow.pl SNP_catalog_file pileupfile(s) outfile";
my $SNP_catalog = shift @filelist;
my $outfile = pop @filelist;

#create hash of SNP catalog information
my %SNPhash;
my $SNP_hash_ref = \%SNPhash;
input_SNPs($SNP_catalog, $SNP_hash_ref);

#read in allele counts at each SNP from all pileup files
foreach my $pileup (@filelist) {
	print "Starting counts of pileup file $pileup \n";
	count_pileup_SNPs ($pileup, $SNP_hash_ref);
	print "Finished counts of pileup file $pileup \n";
}

#print out allele counts and allele frequency differences
SNPprinter($SNP_hash_ref, \@filelist, $outfile);
	
exit;


sub input_SNPs {
	my $SNP_file = shift;
	my $SNP_hash_reference = shift;
	open (SNP, $SNP_file);
	
	#scroll through SNP catalog and enter SNP info into nested hash data structure
	#keys of outermost hash = scaffolds
	#keys of mid-level hash = positions
	#keys of innermost hash = SNP_1, SNP_2, SNP_1_count_filex, etc
	foreach my $entry (<SNP>) {
		chomp $entry; $entry =~ s/\s/\t/ig;
		my @SNPinfo = split ("\t",$entry);
		next if scalar @SNPinfo > 4;     #skips entries with more than two alleles or the two alleles differ from the reference
		my $scaffnum = $1 if ($SNPinfo[0] =~ /scaffold_(\d*)/);
		
		#extend scaffold and position numbers with zeros to facilitate later sorting
		while (length $scaffnum < 4) {
			$scaffnum = "0" . $scaffnum;
		}
		$SNPinfo[0] = "Scaffold_$scaffnum";
		while (length $SNPinfo[1] < 9) {
			$SNPinfo[1] = "0" . $SNPinfo[1];
		}
		
		${$SNP_hash_reference}{$SNPinfo[0]}{$SNPinfo[1]}{'SNP_1'} = $SNPinfo[2];
		${$SNP_hash_reference}{$SNPinfo[0]}{$SNPinfo[1]}{'SNP_2'} = $SNPinfo[3];
	}
	close SNP;
}

sub count_pileup_SNPs {
	my $pileup_file = shift;
	my $SNP_hash_reference = shift;
	my $file_root = $1 if $pileup_file =~ /(.+)\.pileup/;
	my $SNP_1_count = "snp_1_count_" . $file_root;
	my $SNP_2_count = "snp_2_count_" . $file_root;
	open (PUP, $pileup_file);
	foreach my $entry (<PUP>) {
		chomp $entry; $entry =~ s/\s/\t/ig;
		my @position_info = split ("\t",$entry);
		my $scaffnum = $1 if ($position_info[0] =~ /scaffold_(\d*)/);
		
		#extend scaffold and position numbers with zeros to facilitate later sorting
		while (length $scaffnum < 4) {
			$scaffnum = "0" . $scaffnum;
		}
		$position_info[0] = "Scaffold_$scaffnum";
		while (length $position_info[1] < 9) {
			$position_info[1] = "0" . $position_info[1];
		}

		if (exists ${$SNP_hash_reference}{$position_info[0]}{$position_info[1]}) {
				my $scaff = $position_info[0];
				my $pos = $position_info[1];
				my $consbase = $position_info[2];
				my $SNP_1 = ${$SNP_hash_reference}{$scaff}{$pos}{'SNP_1'};
				my $SNP_2 = ${$SNP_hash_reference}{$scaff}{$pos}{'SNP_2'};
				$position_info[8] =~ s/[^AGCTagct\.,]//ig;		#remove nonstandard chars
				$position_info[8] =~ s/[\.,]/$consbase/ig;		#replace consensus with actual base
				$position_info[8] = uc $position_info[8];		#change all to uppercase
				${$SNP_hash_reference}{$scaff}{$pos}{$SNP_1_count} = $position_info[8] =~ s/$SNP_1/$SNP_1/g;
				${$SNP_hash_reference}{$scaff}{$pos}{$SNP_2_count} = $position_info[8] =~ s/$SNP_2/$SNP_2/g;
				${$SNP_hash_reference}{$scaff}{$pos}{$SNP_1_count} = 0 unless (${$SNP_hash_reference}{$scaff}{$pos}{$SNP_1_count});
				${$SNP_hash_reference}{$scaff}{$pos}{$SNP_2_count} = 0 unless (${$SNP_hash_reference}{$scaff}{$pos}{$SNP_2_count});
			} else {
				next;
			}
		}
	close PUP;
}

sub SNPprinter {
	my $SNP_hash_reference = shift;
	my $file_list_ref = shift;
	my @file_list = @{$file_list_ref};
	my $output_file = shift;
	my %window_hash;
	my @count_keys; 
	open (OUT, ">$output_file");
	
	#This loops stores the keys for the SNP counts in the array @count_keys
	#The counters hash, which will temporarily the number of alleles from each file for a given window, is also initiated.
	foreach my $filename (@file_list) {
		my $file_root = $1 if $filename =~ /(.+)\.pileup/;
		my $SNP_1_count = "snp_1_count_" . $file_root;
		my $SNP_2_count = "snp_2_count_" . $file_root;
		push @count_keys, $SNP_1_count;
		push @count_keys, $SNP_2_count;
	}
		
	#Outer loops scroll through the SNP_hash by position
	foreach my $scaff (sort keys %{$SNP_hash_reference}) {
		foreach my $pos (sort keys %{${$SNP_hash_reference}{$scaff}}) {
			my $flag = 0;
			my @counts_array;
			foreach my $index (@count_keys) {
				if (exists ${$SNP_hash_reference}{$scaff}{$pos}{$index}) {
					$flag++;
				}
			}
			next if ($flag < 1);
	#This loop print the counts for each SNP from each file
			foreach my $index (@count_keys) {
				if (exists ${$SNP_hash_reference}{$scaff}{$pos}{$index}) {
					push (@counts_array, ${$SNP_hash_reference}{$scaff}{$pos}{$index});
				} else {
					push (@counts_array, "0");
				}
			}
			next if (($counts_array[0] == 0 && $counts_array[1] == 0) || ($counts_array[2] == 0 && $counts_array[3] == 0));
			my $freq1; my $freq2; my $total_count; my $diff;
			$total_count = $counts_array[0] + $counts_array[1] + $counts_array[2] + $counts_array[3];
			next if $total_count < 20; 
			next if $total_count > 160; 
			next if ($counts_array[0] + $counts_array[1]) <= 10;
			next if ($counts_array[2] + $counts_array[3]) <= 10;
			next if ((($counts_array[1] < 3) && ($counts_array [3] < 3)) || (($counts_array[0] < 3) && ($counts_array [2] < 3)));
			if ($counts_array[0] >= $counts_array[1]) {
				$freq1 = sprintf("%.2f", $counts_array[0]/($counts_array[0]+$counts_array[1]));
				$freq2 = sprintf("%.2f", $counts_array[2]/($counts_array[2]+$counts_array[3]));
				$diff = $freq1 - $freq2;
			} else {
				$freq1 = sprintf("%.2f", $counts_array[1]/($counts_array[0]+$counts_array[1]));
				$freq2 = sprintf("%.2f", $counts_array[3]/($counts_array[2]+$counts_array[3]));
				$diff = $freq1 - $freq2;
			}
			print OUT "$scaff\t$pos\t$counts_array[0]\t$counts_array[1]\t$counts_array[2]\t$counts_array[3]\t";
			print OUT "$total_count\t$freq1\t$freq2\t$diff\n";
			}
	}
}

