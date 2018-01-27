#!/usr/bin/perl
#Copy (C) 2015-2016  Helsinki University.
#Written by Ji Zhang, MD, PhD and Mirko Rossi, DVM, PhD

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY. See the GNU General Public License for 
#more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Revision notes:
#version: 2.1
#4.10.2017
#bugs fixed

#version: 2.0
#5.5.2015
#adding switch -m to save all the temporary files for diagnostic purpose.
#adding switch -n to allow using a multi-Fasta file of the the allele sequences (nt) as reference (a fix number of 999999999 will be assigned to each expect-d).
#bugs fixed.

#version: 1.0
#the first version.


# Switch:
#   -h: Show help
#   -v: Show progress on the screen. [Default = off].
#   -b: Do not output a file for BAPS.
#   -a: Do not output core genome files (Fasta format and eXtended Multi-Fasta format).
#   -l: Disable hyperlink in HTML output (to reduce the file size) [Default = off].
#   -m: Saving all the temporary files to a folder.
#   -n: Using a multi-Fasta file of the the allele sequences (nt) as reference.

# Option:
#   -o: Run the analysis with exiting wgMLST scheme.
#   -g: Name of the text file listing the genomic sequences need to be analyzed. Each file name should occupy a line.
#   -r: Name of the reference genome sequence (GenBank format; or multi-Fasta format (nt) in case the switch -n is on).
#   -c: Minimum coverage of alignment to define an allele. [Default = 100].
#   -i: Minimum identity percentage to define an allele. [Default = 80].
#   -t: Minimum coverage of alignment to define a truncated allele. [Default = 50].

use strict;
use Getopt::Std;
use Benchmark; 
my $T0 = time();

getopts('g:r:d:c:i:t:v b a h l o m n ');
our ($opt_g, $opt_r, $opt_d, $opt_c, $opt_i, $opt_t, $opt_v, $opt_b, $opt_a, $opt_h, $opt_l, $opt_o, $opt_m, $opt_n);
my $distance_d = $opt_d || "20000";
my $coverage_c = $opt_c || "100";
my $identity_i = $opt_i || "80";
my $coverage_t = $opt_t || "50";

if($opt_h==1){
print ("
Genome Profiler - extraction of allele profiles from genomic sequences.
version: 2.1
4.10.2017

Example commands: 

If you run the data for the first time, you could use one of the genome as reference to built a new scheme (ad hoc mode):
perl GeP.pl -r reference.gbk -g genome_list.txt

Or you could use a multi-Fasta file of the the allele sequences (nt) as reference (a fix number of 999999999 will be assigned to expect-d):
perl GeP.pl -r reference.ffn -g genome_list.txt -n

Or you could run your data with the wgMLST scheme created previously by GeP:
perl GeP.pl -g genome_list.txt -o

 Switch:
   -h: Show help.
   -v: Do not show progress on the screen.
   -b: Do not output a file for BAPS.
   -a: Do not output core genome files (Fasta format and eXtended Multi-Fasta format).
   -l: Disable hyperlink in HTML output (to reduce the file size) [Default = off].
   -m: Saving all the temporary files to a folder.
   -n: Using a multi-FASTA file of the the allele sequences as reference.

 Option:
   -o: Run the analysis with exiting wgMLST scheme.
   -g: Name of the text file listing the genomic sequences need to be analyzed. Each file name should occupy a line.
   -r: Name of the reference (default: GenBank format; or multi-Fasta format (nt) in case the switch -n is on).
   -c: Minimum coverage of alignment to define an allele. [Default = 100].
   -i: Minimum identity percentage to define an allele. [Default = 80].
   -t: Minimum coverage of alignment to define a truncated allele. [Default = 50].
   
");
exit;
}

print ("Making preparations for the analysis. Please wait ...\n") if ($opt_v == 0);
###################
open(GENOME, "<$opt_g") or die "Cannot open genome list file!";
open(OUT, ">genome_list.tmp");
open(REPORT, ">>output.report.tmp");
my $counter = 0;
while(<GENOME>){
	if($_ =~ /^\s/){
		next;
	}
	else{
		$counter++;
		chomp;
		print OUT "$_\n";
	}
}
print REPORT "Number of scanned genomes: ", "$counter<br>";
close GENOME;
close OUT;
close REPORT;
###################
my($rep_list);
if($opt_o==0){
	if($opt_n==1){
		print ("GeP is running in ad hoc mode ...\n") if ($opt_v == 0);
		open(IN, "<$opt_r") or die "Cannot open multi-fasta file!";
		open(DIC, ">>output.gene_list.txt");
		open(LIST, ">>gene_list.tmp");
		open(REPORT, ">>output.report.tmp");
		open(SCHEME, ">scheme_summary.txt");
		$/=">";
		my $counter = 0;
		foreach(<IN>){
			if($_ =~ /^>/){
				next;
			}else{
				my($head, $ntseq, $ntseqlen, $count_nt, $aa_base, $codon, $aaseq, $aaseqlen);
				my(@nt);
				my(%aacode);
				$counter++;
				open(NT, ">nt.gene$counter.fas");
				open(AA, ">aa.gene$counter.fas");
				$head = $_;
				$head =~ s/\n.*//gs;
				$ntseq = $_;
				$ntseq =~ s/.*\n//;
				$ntseq =~ s/>//;
				$ntseq =~ s/\s//gs;
				$ntseq =~ lc($ntseq);
				$ntseqlen = length($ntseq);
				
				%aacode = (
#				The Bacterial, Archaeal and Plant Plastid Code (transl_table=11).
				TTT => "F", TCT => "S", TAT => "Y", TGT => "C",
				TTC => "F", TCC => "S", TAC => "Y", TGC => "C",
				TTA => "L", TCA => "S", TAA => "*", TGA => "*",
				TTG => "L", TCG => "S", TAG => "*", TGG => "W",
				CTT => "L", CCT => "P", CAT => "H", CGT => "R",
				CTC => "L", CCC => "P", CAC => "H", CGC => "R",
				CTA => "L", CCA => "P", CAA => "Q", CGA => "R",
				CTG => "L", CCG => "P", CAG => "Q", CGG => "R",
				ATT => "I", ACT => "T", AAT => "N", AGT => "S",
				ATC => "I", ACC => "T", AAC => "N", AGC => "S",
				ATA => "I", ACA => "T", AAA => "K", AGA => "R",
				ATG => "M", ACG => "T", AAG => "K", AGG => "R",
				GTT => "V", GCT => "A", GAT => "D", GGT => "G",
				GTC => "V", GCC => "A", GAC => "D", GGC => "G",
				GTA => "V", GCA => "A", GAA => "E", GGA => "G",
				GTG => "V", GCG => "A", GAG => "E", GGG => "G",
				);
				
				@nt = split "", $ntseq;
				$count_nt = 0;
				foreach(@nt){
					if($count_nt==3){
						$aa_base = $aacode{$codon};
						$aaseq = $aaseq.$aa_base;
						$count_nt = 1;
						$codon = $_;
					}else{
						$count_nt++;
						$codon = $codon.$_;
					}
				}
				
				$aaseqlen = length($aaseq);
				
				print NT ">gene$counter", "_$ntseqlen", "_1\n", "$ntseq\n";
				print AA ">gene$counter", "_$aaseqlen\n$aaseq\n";
				print DIC "gene$counter\t$head\n";
				print LIST "gene$counter\t", "999999999\n";
				system("makeblastdb -in aa.gene$counter.fas -dbtype prot -logfile makeblastdb.log");
				system("makeblastdb -in nt.gene$counter.fas -dbtype nucl -logfile makeblastdb.log");
				close AA;
				close NT;
			}
		}
		$/="\n";
		print REPORT "Number of the reference genes: ", "$counter<br>";
		print SCHEME "Number of the reference genes: ", "$counter\n";
		close IN;
		close DIC;
		close LIST;
		close REPORT;
		close SCHEME;
		system ("mkdir scheme");
		system("cp gene_list.tmp ./scheme/gene_list.tmp");
	}else{
		print ("GeP is running in ad hoc mode ...\n") if ($opt_v == 0);
		my($counter, $counter2, $poz, $contig_id, $strand, $gene_type, $start, $stop, $locus_tag, $product, $nt, $nt_seq, $len, $aa_seq);
		my(@in, @strands, @starts, @stops, @locus_tags, @products, @aa_seqs, @gene_types);
		open(IN, "<$opt_r") or die "Cannot open reference gbk file!";
		open(OUT, ">>$opt_r.tmp");
		print OUT ("contig_id	start	stop	locus_tag	product	nt_seq	aa_seq");
		print OUT "\n";
		$/="\n//";
		@in = <IN>;
		$/="\n";
		$counter = 0;
		foreach(@in){
			$counter++;
			open(TMP, ">$opt_r.$counter.tmp");
			$_ =~ s/\n//gs;
			$_ =~ s/  +gene  +/\n$&/gs;
			$_ =~ s/  +CDS  +/\n$&/gs;
			$_ =~ s/  +rRNA  +/\n$&/gs;
			$_ =~ s/  +tRNA  +/\n$&/gs;
			$_ =~ s/ORIGIN  +1/\n$&/;
			print TMP "$_";
			close TMP;
			
			@strands=(); @starts=(); @stops=(); @locus_tags=(); @products=(); @aa_seqs=(); @gene_types=();
			open(TMP, "<$opt_r.$counter.tmp");
			while(<TMP>){
				chomp;
				$contig_id = "$opt_r.$counter";
				if($_ =~ /^  +CDS  +/){
					chomp;
					$aa_seq = $_;
					if($aa_seq =~ m/  +\/translation\=\"/){
						$aa_seq =~ s/.*  +\/translation\=\"//;
						$aa_seq =~ s/\".*//;
						$aa_seq =~ s/ //g;
						$poz = $_;
						$poz =~ s/\>//g;
						$poz =~ s/\<//g;
						if($poz =~ m/  +complement\([0-9]+\.\.[0-9]+\)/){
							$strand = "minus";
						}
						else{
							$strand = "plus";
						}
						$poz =~ /[0-9]+\.\.[0-9]+/;
						$poz = $&;
						$start = $poz;
						$stop = $poz;
						$stop =~ s/[0-9]+\.\.//;
						$start =~ s/\.\.[0-9]+//;
						$gene_type = "CDS";
						
						if($_ =~ m/  \/locus_tag\=\"/){
							$locus_tag = $_;
							$locus_tag =~ s/.*  \/locus_tag\=\"//;
							$locus_tag =~ s/\".*//;
						}
						elsif($_ =~ m/[0-9]+\.peg\.[0-9]+/){
							$locus_tag = $_;
							$locus_tag =~ m/[0-9]+\.peg\.[0-9]+/;
							$locus_tag = $&;
							$locus_tag =~ s/.*peg/peg/;
						}
						else{
							$locus_tag = "NA";
						}
						
						$product = $_;
						$product =~ s/.*  +\/product\=\"//;
						$product =~ s/\".*//;
						$product =~ s/ +/ /g;
						
						push @strands, $strand;
						push @gene_types, $gene_type;
						push @starts, $start;
						push @stops, $stop;
						push @locus_tags, $locus_tag;
						push @aa_seqs, $aa_seq;
						push @products, $product;
					}
					else{
						next;
					}
				}
				elsif($_ =~ /^ORIGIN  +1/){
					$nt = $_;
					$nt =~ s/$&//;
					$nt =~ s/[0-9]+//gs;
					$nt =~ s/ //gs;
					$nt =~ s/\/\///;
				}
				else{
					next;
				}
			}
			close TMP;
			system("rm -f $opt_r.$counter.tmp");
			
			$counter2 = -1;
			foreach(@starts){
				$counter2++;
				$start = $_;
				$stop = $stops["$counter2"];
				$gene_type = $gene_types["$counter2"];
				$locus_tag = $locus_tags["$counter2"];
				$product = $products["$counter2"];
				$aa_seq = $aa_seqs["$counter2"];
				$strand = $strands["$counter2"];
				$start = $start-1;
				$len = $stop-$start;
				$nt_seq = substr $nt, $start, $len;
				if($strand eq "minus"){
					$nt_seq = reverse($nt_seq);
					$nt_seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
					$start = $start+1;
				}
				else{
					$start = $start+1;
				}
				
				if($gene_type eq "CDS"){
					print OUT "$contig_id\t$start\t$stop\t$locus_tag\t$product\t$nt_seq\t$aa_seq\n";
				}
				else{
					next;
				}
			}
		}
		close IN;
		close OUT;
		
		open(REF, "<$opt_r.tmp");
		open(LIST, ">>gene_list.tmp");
		open(LIST2, ">>output.gene_list.txt");
		open(REPORT, ">>output.report.tmp");
		open(SCHEME, ">scheme_summary.txt");
		open(FFN, ">$opt_r.ffn.tmp");
		open(FNA, ">$opt_r.fna.tmp");
		print FNA ">$opt_r.fna\n";
		my($counter, $contig_id_new, $contig_id_old, $poz1, $poz2, $locus_tag, $product, $annotation, $nt, $aa, $ntlength, $aalength,
				$start_new, $start_old, $end_new, $end_old, $distance);
		my(@in);
		$counter = -1;
		while(<REF>){
			$counter++;
			chomp;
			if($counter == 0){
				next;
				}
			else{
				@in = split (/\t/,$_);
				$contig_id_new = $in[0];
				$poz1 = $in[1];
				$poz2 = $in[2];
				$locus_tag = $in[3];
				$product = $in[4];
				$annotation = "$product"."; "."$locus_tag";
				$nt = $in[5];
				$aa = $in[6];
				$ntlength = length($nt);
				$aalength = length($aa);
				
				if($poz1<$poz2){
					$start_new = $poz1;
					$end_new = $poz2;
				}
				else{
					$start_new = $poz2;
					$end_new = $poz1;
				}
				
				if($contig_id_new eq $contig_id_old){
					$distance = $start_new - $end_old;
					if($distance<0){
						$distance = 0;
					}
					$distance =$distance + 10;
				}
				else{
					$distance = -1;
				}
				
				$contig_id_old = $contig_id_new;
				$start_old = $start_new;
				$end_old = $end_new;
				
				open(NT, ">nt.gene$counter.fas");
				print NT ">gene", "$counter", "_$ntlength", "_1", "\n", "$nt", "\n";
				print FFN ">gene", "$counter", "_$ntlength", "_1", "\n", "$nt", "\n";
				print FNA "$nt", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
				close NT;
				
				open(AA, ">aa.gene$counter.fas");
				print AA ">gene", "$counter", "_$aalength", "\n", "$aa", "\n";
				close AA;
				
				print LIST "gene$counter\t", "$distance\n";
				print LIST2 "gene$counter\t", "$annotation\n";
				
				system ("makeblastdb -in aa.gene$counter.fas -dbtype prot -logfile makeblastdb.log");
				system ("makeblastdb -in nt.gene$counter.fas -dbtype nucl -logfile makeblastdb.log");
			}
		}
		print REPORT "Number of the reference genes: ", "$counter<br>";
		print SCHEME "Number of the reference genes: ", "$counter\n";
		close LIST;
		close LIST2;
		close REPORT;
		close SCHEME;
		close REF;
		###################checking for multi-copy genes in the reference genome
		print ("checking for multi-copy genes in the reference genome ...\n") if ($opt_v == 0);
		system ("makeblastdb -in $opt_r.ffn.tmp -dbtype nucl -logfile makeblastdb.log");
		open(SCHEME, ">>scheme_summary.txt");
		my($sseqid, $sseqid_prev, $identity, $allele_len, $s_start, $s_end, $q_start, $q_end,
			$s_len, $coverage, $gene_id, $gene_id_prev, $rep_list, $rep_sum, $counter);
		my(@blastn, @blastn1, @rep_list);
		@blastn = readpipe("blastn -query $opt_r.fna.tmp -db $opt_r.ffn.tmp -num_threads 2 -dust no -penalty -3 -reward 2 -max_target_seqs 99999 -word_size 11 -outfmt 6");
		foreach(@blastn){
			chomp;
			@blastn1 = split (/\t/,$_);
			$gene_id = $blastn1[1];
			$gene_id =~ s/_.*//g;
			$allele_len = $blastn1[1];
			$allele_len =~ /_[0-9]+_/;
			$allele_len = $&;
			$allele_len =~ s/_//g;
			$identity = $blastn1[2];
			$q_start = $blastn1[6];
			$q_end = $blastn1[7];
			$s_start = $blastn1[8];
			$s_end = $blastn1[9];
			$s_len = abs($s_end-$s_start)+1;
			$coverage = $s_len/$allele_len*100;
			if($coverage>=$coverage_t&&$identity>=$identity_i){
				if($gene_id eq $gene_id_prev){
					$counter++;
					unless($counter>1){
						push @rep_list, $gene_id;
					}
					$gene_id_prev = $gene_id;
				}
				else{
					$counter = 0;
					$gene_id_prev = $gene_id;
				}
			}
			else{
				next;
			}
		}
		$rep_list = join (",", @rep_list);
		$rep_sum = @rep_list;
		print ("$rep_sum possible multi-copy genes ($rep_list) in the reference genome are recorded.\n") if ($opt_v == 0);
		$rep_list = "$rep_list".",";
		print SCHEME "multi-copy gene = ", "$rep_list";
		close SCHEME;
		system ("mkdir scheme");
		system("cp gene_list.tmp ./scheme/gene_list.tmp");
	}
}
else{
	print ("GeP is running with existing wgMLST scheme ...\n") if ($opt_v == 0);
	system("cp ./scheme/gene_list.tmp ./");
	system("mv ./scheme/aa.* ./");
	system("mv ./scheme/nt.* ./");
	system("mv ./scheme/scheme_summary.txt ./");
	system("mv ./scheme/output.gene_list.txt ./");
	open(IN, "<scheme_summary.txt");
	my($line, $rep_list, $loci_sum, $loci_counter);
	$line=0;
	foreach(<IN>){
		chomp;
		$line++;
		if($line==1){
			$loci_sum = $_;
			$loci_sum =~ /[0-9]+/;
			$loci_sum = $&;
			$loci_counter = 0;
			while($loci_counter<=$loci_sum){
				$loci_counter++;
				system("makeblastdb -in aa.gene$loci_counter.fas -dbtype prot -logfile makeblastdb.log");
				system("makeblastdb -in nt.gene$loci_counter.fas -dbtype nucl -logfile makeblastdb.log");
			}
			open(REPORT, ">>output.report.tmp");
			print REPORT "$_<br>";
			close REPORT;
		}
		elsif($line==2){
			$rep_list = $_;
			$rep_list =~ s/multi-copy gene \= //;
		}
		else{}
	}
	close IN;
}
###################
my ($genome, $seq, $counter, $spacer);
open(GENOME, "<genome_list.tmp");
while(<GENOME>){
	chomp;
	$genome = $_;
	open(SEQIN, "<$genome") or die "Cannot open genome sequence $genome!";
	open(OUT, ">>$genome.combined.fas");
	$/=">";
	print OUT ">", "$genome\n";
	while(<SEQIN>){
		$spacer = ();
		$counter = 0;
		while($counter<$distance_d){
			$spacer = "$spacer"."N";
			$counter++;
		}
		
		if($_ =~ /^>/){
			next;
		}
		else{
			my $seq = $_;
			$seq =~ s/.*\n//;
			$seq =~ s/>//;
			$seq =~ s/\n//g;
			print OUT "$seq", "\n$spacer\n";
		}
	}
	close OUT;
	close SEQIN;
	$/="\n";
}
close GENOME;
###################
my ($genome, $gene, $nt_db, $allele_new, $line_counter, $q_end_prev, $q_start_prev, 
	$allele_len, $allele, $identity, $aln_len, $mismatch, $gap, $q_start, $q_end, $s_start, $s_end, 
	$s_end2, $q_end2, $s_len, $q_len, $allele_counter, $copy, $d, $allele_old, $allele_counter2, 
	$strand, $coverage, $poz, $poz_edge, $distance, $line_OK, $distance_OK, $poz_OK, $allele_OK, $strand_OK, $coverage_OK, $identity_OK, $q_len_OK,
	$edge1_start, $edge1_len, $edge2_start, $edge2_len, $edge2_end,
	$poz_edge, $poz_edge_OK, $q_len_prev, $total_len, $coding_len, $go_on, $dis_index);
open(GENOME, "<genome_list.tmp");
open(RESULTS, ">>output.txt");
open(BLAST, ">>blast_output.tmp");
print RESULTS "genome_id\t", "gene_id\t", "allele\t", "prev_locus_start\t", "prev_locus_end\t", "hit_start\t", "hit_end\t", "expected_d\t", "actual_d\t", "program\t", "copy\n";
while(<GENOME>){
	chomp;
	$genome = $_;
	$q_start_prev = "NA";
	$q_end_prev = "NA";
	open(GENE, "<gene_list.tmp");
	while(<GENE>){
		chomp;
		$gene = $_;
		$gene =~ s/\t.*//;
		$d = $_;
		$d =~ s/.*\t//;
		my ($perfectn_hit, $perfectn_found, $perfectn_poz, $perfectn_d, $copy, $allele_old);
		my (@blastn, @blastn1, @blastx, @blastx1, @distances, @sorted, @pozs, @pozs_sorted, @lines, @lines_1st, @dises_index); 
		my (%poz_hash_1st, %poz_edge_hash_1st, %distance_hash_1st, %allele_hash_1st, %strand_hash_1st, %coverage_hash_1st, %identity_hash_1st, %q_len_hash_1st,
		%poz_hash, %poz_edge_hash, %distance_hash, %allele_hash, %strand_hash, %coverage_hash, %identity_hash, %q_len_hash);
		
		$allele_old=(); $allele_counter2=0; $go_on = "yes";$perfectn_found = "no";
		print ("blastn: ", "$genome", " vs. ", "$gene", " ...\n") if ($opt_v == 0);
		print BLAST "blastn: $genome vs. $gene\n", "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n";
		@blastn = readpipe("blastn -query $genome.combined.fas -db nt.$gene.fas -outfmt 6 -dust no -penalty -3 -reward 2 -word_size 11");
		open (NT, "<nt.$gene.fas");
		$/=undef;
		$nt_db = <NT>;
		$allele_counter = ($nt_db =~ s/>/>/g);
		$allele_new = $allele_counter+1;
		close NT;
		$/="\n";
		$copy = 0; $line_counter = 0;
		
		foreach(@blastn){
			print BLAST "$_";
			chomp;
			$line_counter++;
			@blastn1 = split (/\t/,$_);
			
			$allele_len = $blastn1[1];
			$allele_len =~ m/[0-9]+_[0-9]+$/;
			$allele_len = $&;
			$allele_len =~ s/_.*//g;
			
			$allele = $blastn1[1];
			$allele =~ s/.*_//g;
			
			$identity = $blastn1[2];
			$aln_len = $blastn1[3];
			$mismatch = $blastn1[4];
			$gap = $blastn1[5];
			$q_start = $blastn1[6];
			$q_end = $blastn1[7];
			$s_start = $blastn1[8];
			$s_end = $blastn1[9];
			
			if($s_end>$s_start){
				$strand = "plus";
			}
			else{
				$strand = "minus";
			}
			
			$s_len = abs($s_end-$s_start)+1;
			$q_len = abs($q_end-$q_start)+1;
			
			$coverage = $s_len/$allele_len;
			$coverage = $coverage*100;
			
			if($strand eq "minus"){
				$edge2_end = $q_end + $s_end -1;
				$edge1_start = $q_start - $allele_len + $s_start;
			}
			else{
				$edge1_start = $q_start - $s_start + 1;
				$edge2_end = $q_end + $allele_len - $s_end;
			}
			
			$poz = "$q_start"."\t"."$q_end";
			$poz_edge = "$edge1_start"."\t"."$edge2_end";
			
			if($coverage>=$coverage_t&&$identity>=$identity_i){
				if($allele ne $allele_old){
					$allele_old = $allele;
					$allele_counter2++;
				}
				else{
					$allele_old = $allele;
				}
			
				@pozs=(); @pozs_sorted=();
				push @pozs, $q_start;
				push @pozs, $q_end;
				push @pozs, $q_start_prev;
				push @pozs, $q_end_prev;
				@pozs_sorted = sort{$a<=>$b}@pozs;
				$distance = $pozs_sorted[2] - $pozs_sorted[1];
				$total_len = $pozs_sorted[3] - $pozs_sorted[0] + 1;
				$coding_len = $q_len + $q_len_prev;
				unless($total_len>$coding_len){
					$distance = 0;
				}
				$dis_index = abs($distance);
				
				if($allele_counter2 == 1){
					$copy++;
					push @dises_index, $dis_index;
					$poz_hash_1st{$dis_index} = $poz;
					$poz_edge_hash_1st{$dis_index} = $poz_edge;
					$distance_hash_1st{$dis_index} = $distance;
					$allele_hash_1st{$dis_index} = $allele;
					$strand_hash_1st{$dis_index} = $strand;
					$coverage_hash_1st{$dis_index} = $coverage;
					$identity_hash_1st{$dis_index} = $identity;
					$q_len_hash_1st{$dis_index} = $q_len;
				}
				else{};
				
				push @lines, $line_counter;
				$poz_hash{$line_counter} = $poz;
				$poz_edge_hash{$line_counter} = $poz_edge;
				$allele_hash{$line_counter} = $allele;
				$strand_hash{$line_counter} = $strand;
				$coverage_hash{$line_counter} = $coverage;
				$identity_hash{$line_counter} = $identity;
				$distance_hash{$line_counter} = $distance;
				$q_len_hash{$line_counter} = $q_len;
			}
			else{
				next;
			}
		}
		
		if($copy>0){
			foreach(@lines){
				chomp;
				$line_OK = $_;
				$poz_OK = $poz_hash{$line_OK};
				$allele_OK = $allele_hash{$line_OK};
				$coverage_OK = $coverage_hash{$line_OK};
				$identity_OK = $identity_hash{$line_OK};
				$distance_OK = $distance_hash{$line_OK};
				$q_len_OK = $q_len_hash{$line_OK};
				
				if($perfectn_found eq "no"&&$coverage_OK==100&&$identity_OK==100){
					$perfectn_hit = $allele_OK;
					$perfectn_found = "yes";
					$perfectn_poz = $poz_OK;
					$perfectn_d = $distance_OK;
				}
				else{}
				
				if($rep_list =~ /$gene,/){
					if($go_on eq "yes"&&$coverage_OK==100&&$identity_OK==100&&$distance_OK<=$d){
						print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $poz_OK;
						$q_end_prev =~ s/.*\t//;
						$q_start_prev = $poz_OK;
						$q_start_prev =~ s/\t.*//;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					else{}
				}
				else{
					if($go_on eq "yes"&&$copy==1&&$coverage_OK==100&&$identity_OK==100){
						print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $poz_OK;
						$q_end_prev =~ s/.*\t//;
						$q_start_prev = $poz_OK;
						$q_start_prev =~ s/\t.*//;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					elsif($go_on eq "yes"&&$copy>1&&$coverage_OK==100&&$identity_OK==100&&$distance_OK<=$d){
						print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $poz_OK;
						$q_end_prev =~ s/.*\t//;
						$q_start_prev = $poz_OK;
						$q_start_prev =~ s/\t.*//;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					else{}
				}
			}
			
			if($go_on eq "yes"){
				chomp;
				@sorted=sort{$a<=>$b}@dises_index;
				$distance_OK = $sorted[0];
				$poz_OK = $poz_hash_1st{$distance_OK};
				$poz_edge_OK = $poz_edge_hash_1st{$distance_OK};
				$allele_OK = $allele_hash_1st{$distance_OK};
				$strand_OK = $strand_hash_1st{$distance_OK};
				$coverage_OK = $coverage_hash_1st{$distance_OK};
				$identity_OK = $identity_hash_1st{$distance_OK};
				$q_len_OK = $q_len_hash_1st{$distance_OK};
				
				my($extract_start, $extract_end, $new_allele_len, $seqin, $extract,
						$extract_edge_start, $extract_edge_end, $extract_edge_len, $extract_edge,
						$extract_start_codon, $new_allele_len_codon, $extract_codon, $extract_codon_len, $codon_start, $extra_codon,
						$last_codon_start, $last_codon, $stoped, $extra_stop);
				open(SEQIN, "<$genome.combined.fas");
				$extract_start = $poz_OK;
				$extract_start =~ s/\t.*//;
				$extract_end = $poz_OK;
				$extract_end =~ s/.*\t//;
				$new_allele_len = $extract_end-$extract_start+1;
				$extract_start = $extract_start-1;
				
				$extract_edge_start = $poz_edge_OK;
				$extract_edge_start =~ s/\t.*//;
				$extract_edge_end = $poz_edge_OK;
				$extract_edge_end =~ s/.*\t//;
				$extract_edge_len = $extract_edge_end-$extract_edge_start + 1;
				$extract_edge_start = $extract_edge_start-1;
				
				$/=undef;
				$seqin = <SEQIN>;
				$seqin =~ s/>.*\n//;
				$seqin =~ s/\s//g;
				$extract = substr $seqin, $extract_start, $new_allele_len;
				$extract_edge = substr $seqin, $extract_edge_start, $extract_edge_len;
				$/="\n";
				
				if($strand_OK eq "minus"){
					$extract = reverse($extract);
					$extract =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
					$extract_edge = reverse($extract_edge);
					$extract_edge =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
					
					$extract_start_codon = $extract_start - 3;
					$new_allele_len_codon = $new_allele_len + 3;
					$extract_codon = substr $seqin, $extract_start_codon, $new_allele_len_codon;
					$extract_codon = reverse($extract_codon);
					$extract_codon =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
				}
				else{
					$extract_start_codon = $extract_start;
					$new_allele_len_codon = $new_allele_len + 3;
					$extract_codon = substr $seqin, $extract_start_codon, $new_allele_len_codon;
				}
				
				$extract_start = $extract_start+1;
				$extract =~ s/\s//g;
				$extract_edge =~ s/\s//g;
				$extract_codon =~ s/\s//g;
				
				$last_codon_start = $new_allele_len - 3;
				$last_codon = substr $extract, $last_codon_start, 3;
				
				$extract_codon_len = length($extract_codon);
				$codon_start = $extract_codon_len - 3;
				$extra_codon = substr $extract_codon, $codon_start, 3;
				
				if($last_codon=~/taa/i or $last_codon=~/tag/i or $last_codon=~/tga/i){
					$stoped = "yes";
				}
				else{
					$stoped = "no";
				}
				
				if($extra_codon=~/taa/i or $extra_codon=~/tag/i or $extra_codon=~/tga/i){
					$extra_stop = "yes";
				}
				else{
					$extra_stop = "no";
				}
				
				if($copy==1&&$coverage_OK==100&&$identity_OK==100){
					print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
					$q_end_prev = $poz_OK;
					$q_end_prev =~ s/.*\t//;
					$q_start_prev = $poz_OK;
					$q_start_prev =~ s/\t.*//;
					$q_len_prev = $q_len_OK;
					$go_on = "no";
				}
				elsif($copy==1&&$coverage_OK>=$coverage_c&&$identity_OK != 100){
					if($go_on eq "yes"&&$extract =~ /[RYSWKMBDHVNX]/i){
						print RESULTS "$genome\t", "$gene\t", "n\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $extract_end;
						$q_start_prev = $extract_start;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					elsif($go_on eq "yes"&&$extract_edge =~ /N+/i){
						print RESULTS "$genome\t", "$gene\t", "N\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $extract_end;
						$q_start_prev = $extract_start;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					elsif($go_on eq "yes"){
						if($extra_stop eq "yes" && $stoped eq "no"){
							open(NTDB, ">>nt.$gene.fas");
							print NTDB "\n>", "$gene", "_", "$new_allele_len_codon", "_","$allele_new\n", "$extract_codon";
							system ("makeblastdb -in nt.$gene.fas -dbtype nucl -logfile makeblastdb.log");
							print RESULTS "$genome\t", "$gene\t", "$allele_new\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
							if($strand_OK==0){
							$q_end_prev = $extract_end;
							$q_start_prev = $extract_start - 3;
							$q_len_prev = $q_len_OK;
							}
							else{
							$q_end_prev = $extract_end +3;
							$q_start_prev = $extract_start;
							$q_len_prev = $q_len_OK;
							}
							close NTDB;
							$go_on = "no";
						}
						else{
							open(NTDB, ">>nt.$gene.fas");
							print NTDB "\n>", "$gene", "_", "$new_allele_len", "_","$allele_new\n", "$extract";
							system ("makeblastdb -in nt.$gene.fas -dbtype nucl -logfile makeblastdb.log");
							print RESULTS "$genome\t", "$gene\t", "$allele_new\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
							$q_end_prev = $extract_end;
							$q_start_prev = $extract_start;
							$q_len_prev = $q_len_OK;
							close NTDB;
							$go_on = "no";
						}
					}
					else{};
				}
				elsif($copy>1&&$coverage_OK==100&&$identity_OK==100&&$distance_OK<=$d){
					print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
					$q_end_prev = $poz_OK;
					$q_end_prev =~ s/.*\t//;
					$q_start_prev = $poz_OK;
					$q_start_prev =~ s/\t.*//;
					$q_len_prev = $q_len_OK;
					$go_on = "no";
				}
				elsif($copy>1&&$coverage_OK>=$coverage_c&&$distance_OK<=$d){
					if($go_on eq "yes"&&$extract =~ /[RYSWKMBDHVNX]/i){
						print RESULTS "$genome\t", "$gene\t", "n\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $extract_end;
						$q_start_prev = $extract_start;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					elsif($go_on eq "yes"&&$extract_edge =~ /N+/i){
						print RESULTS "$genome\t", "$gene\t", "N\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
						$q_end_prev = $extract_end;
						$q_start_prev = $extract_start;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					elsif($go_on eq "yes"){
						if($extra_stop eq "yes" && $stoped eq "no"){
							open(NTDB, ">>nt.$gene.fas");
							print NTDB "\n>", "$gene", "_", "$new_allele_len_codon", "_","$allele_new\n", "$extract_codon";
							system ("makeblastdb -in nt.$gene.fas -dbtype nucl -logfile makeblastdb.log");
							print RESULTS "$genome\t", "$gene\t", "$allele_new\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
							if($strand_OK==0){
							$q_end_prev = $extract_end;
							$q_start_prev = $extract_start - 3;
							$q_len_prev = $q_len_OK;
							}
							else{
							$q_end_prev = $extract_end +3;
							$q_start_prev = $extract_start;
							$q_len_prev = $q_len_OK;
							}
							close NTDB;
							$go_on = "no";
						}
						else{
							open(NTDB, ">>nt.$gene.fas");
							print NTDB "\n>", "$gene", "_", "$new_allele_len", "_","$allele_new\n", "$extract";
							system ("makeblastdb -in nt.$gene.fas -dbtype nucl -logfile makeblastdb.log");
							print RESULTS "$genome\t", "$gene\t", "$allele_new\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastn\t", "$copy\n";
							$q_end_prev = $extract_end;
							$q_start_prev = $extract_start;
							$q_len_prev = $q_len_OK;
							close NTDB;
							$go_on = "no";
						}
					}
					else{};
				}
				else{
					my($line_counter, $copy_x, $extract_start, $extract_end, $new_allele_len, $seqin, $extract,
						$extract_edge_start, $extract_edge_end, $extract_edge_len, $extract_edge, $dis_index);
					my(@blastx, @blastx1, @distances, @sorted, @pozs, @pozs_sorted, @lines, @lines_1st, @dises_index); 
					my(%poz_hash_1st, %poz_edge_hash_1st, %distance_hash_1st, %allele_hash_1st, %strand_hash_1st, %coverage_hash_1st, %identity_hash_1st,
						%poz_hash, %poz_edge_hash, %distance_hash, %allele_hash, %strand_hash, %coverage_hash, %identity_hash);
					print ("blastx: ", "$genome", " vs. ", "$gene", " ...\n") if ($opt_v == 0);
					print BLAST "blastx: $genome vs. $gene\n", "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n";
					@blastx = readpipe("blastx -query $genome.combined.fas -db aa.$gene.fas -gapopen 11 -seg no -outfmt 6");
					$copy_x=0;
					if(@blastx){
						foreach (@blastx){
							print BLAST "$_";
							chomp;
							$line_counter++;
							@blastx1 = split (/\t/,$_);
							$allele_len = $blastx1[1];
							$allele_len =~ s/.*_//g;
							
							$identity = $blastx1[2];
							$aln_len = $blastx1[3];
							$mismatch = $blastx1[4];
							$gap = $blastx1[5];
							$q_start = $blastx1[6];
							$q_end = $blastx1[7];
							$s_start = $blastx1[8];
							$s_end = $blastx1[9];
							
							if($q_end>$q_start){
								$strand = "plus";
								}
							else{
								$strand = "minus";
								$q_end2 = $q_end;
								$q_end = $q_start;
								$q_start = $q_end2;
							}
							
							$s_len = abs($s_end-$s_start)+1;
							$q_len = abs($q_end-$q_start)+1;
							
							$coverage = $s_len/$allele_len;
							$coverage = $coverage*100;
							
							$allele_len = $allele_len*3 + 3;
							
							unless($s_start==1){
								$s_start = 3*($s_start-1)+1;
							}
							$s_end = 3*$s_end;
							
							if($strand eq "minus"){
								$edge1_start = $q_start - $allele_len + $s_end;
								$edge2_end = $q_end + $s_start - 1;
							}
							else{
								$edge1_start = $q_start - $s_start + 1;
								$edge2_end = $q_end + $allele_len - $s_end;
							}
							
							$poz = "$q_start"."\t"."$q_end";
							$poz_edge = "$edge1_start"."\t"."$edge2_end";
							
							if($identity>=$identity_i&&$coverage>=$coverage_t){
								$copy_x++;
								if($coverage>=$coverage_c){
									$allele = $allele_new;
								}
								else{
									$allele = "T";
								}
								
								@pozs=(); @pozs_sorted=();
								push @pozs, $q_start;
								push @pozs, $q_end;
								push @pozs, $q_start_prev;
								push @pozs, $q_end_prev;
								@pozs_sorted = sort{$a<=>$b}@pozs;
								$distance = $pozs_sorted[2] - $pozs_sorted[1];
								$total_len = $pozs_sorted[3] - $pozs_sorted[0] + 1;
								$coding_len = $q_len + $q_len_prev;
								unless($total_len>$coding_len){
									$distance = 0;
								}
								$dis_index = abs($distance);
								
								push @dises_index, $dis_index;
								$poz_hash{$dis_index} = $poz;
								$poz_edge_hash{$dis_index} = $poz_edge;
								$allele_hash{$dis_index} = $allele;
								$strand_hash{$dis_index} = $strand;
								$coverage_hash{$dis_index} = $coverage;
								$identity_hash{$dis_index} = $identity;
								$distance_hash{$dis_index} = $distance;
								$q_len_hash{$dis_index} = $q_len;
							}
							else{}
						}
						
						if(@dises_index){
							@sorted=sort{$a<=>$b}@dises_index;
							$distance_OK = $sorted[0];
							$poz_OK = $poz_hash{$distance_OK};
							$poz_edge_OK = $poz_edge_hash{$distance_OK};
							$allele_OK = $allele_hash{$distance_OK};
							$coverage_OK = $coverage_hash{$distance_OK};
							$identity_OK = $identity_hash{$distance_OK};
							$strand_OK = $strand_hash{$distance_OK};
							$distance_OK = $distance_hash{$distance_OK};
							$q_len_OK = $q_len_hash{$distance_OK};
							
							my($extract_start, $extract_end, $new_allele_len, $seqin, $extract,
									$extract_edge_start, $extract_edge_end, $extract_edge_len, $extract_edge,
									$extract_start_codon, $new_allele_len_codon, $extract_codon, $extract_codon_len, $codon_start, $extra_codon);
							open(SEQIN, "<$genome.combined.fas");
							$extract_start = $poz_OK;
							$extract_start =~ s/\t.*//;
							$extract_end = $poz_OK;
							$extract_end =~ s/.*\t//;
							$new_allele_len = $extract_end-$extract_start+1;
							$extract_start = $extract_start-1;
							
							$extract_edge_start = $poz_edge_OK;
							$extract_edge_start =~ s/\t.*//;
							$extract_edge_end = $poz_edge_OK;
							$extract_edge_end =~ s/.*\t//;
							$extract_edge_len = $extract_edge_end-$extract_edge_start + 1;
							$extract_edge_start = $extract_edge_start-1;
							
							$/=undef;
							$seqin = <SEQIN>;
							$seqin =~ s/>.*\n//;
							$seqin =~ s/\s//g;
							$extract = substr $seqin, $extract_start, $new_allele_len;
							$extract_edge = substr $seqin, $extract_edge_start, $extract_edge_len;
							$/="\n";
							
							if($strand_OK eq "minus"){
								$extract = reverse($extract);
								$extract =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
								$extract_edge = reverse($extract_edge);
								$extract_edge =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
								
								$extract_start_codon = $extract_start - 3;
								$new_allele_len_codon = $new_allele_len + 3;
								$extract_codon = substr $seqin, $extract_start_codon, $new_allele_len_codon;
								$extract_codon = reverse($extract_codon);
								$extract_codon =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
							}
							else{
								$extract_start_codon = $extract_start;
								$new_allele_len_codon = $new_allele_len + 3;
								$extract_codon = substr $seqin, $extract_start_codon, $new_allele_len_codon;
							}
							
							$extract_start = $extract_start+1;
							$extract =~ s/\s//g;
							$extract_edge =~ s/\s//g;
							$extract_codon =~ s/\s//g;
							
							$extract_codon_len = length($extract_codon);
							$codon_start = $extract_codon_len - 3;
							$extra_codon = substr $extract_codon, $codon_start, 3;
							
							if($go_on eq "yes"&&$perfectn_found eq "yes"&&$copy_x==1){
								print RESULTS "$genome\t", "$gene\t", "$perfectn_hit\t", "$q_start_prev\t", "$q_end_prev\t", "$perfectn_poz\t", "$d\t", "$perfectn_d\t", "blastn\t", "$copy\n";
								$q_end_prev = $perfectn_poz;
								$q_end_prev =~ s/.*\t//;
								$q_start_prev = $perfectn_poz;
								$q_start_prev =~ s/\t.*//;
								$q_len_prev = $q_len_OK;
								$go_on = "no";
							}
							elsif($go_on eq "yes"&&$copy_x>1&&$distance_OK>$d){
								print RESULTS "$genome\t", "$gene\t", "D\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
								$go_on = "no";
								$q_end_prev = "NA";
								$q_start_prev = "NA";
								$q_len_prev = $q_len_OK;
							}
							elsif($go_on eq "yes"){
								if($allele_OK eq "T"){
									print RESULTS "$genome\t", "$gene\t", "T\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
									$q_end_prev = $extract_end;
									$q_start_prev = $extract_start;
									$q_len_prev = $q_len_OK;
									$go_on = "no";
								}
								elsif($allele_OK==$allele_new){
									if($go_on eq "yes"&&$extract =~ /[RYSWKMBDHVNX]/i){
										print RESULTS "$genome\t", "$gene\t", "n\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
										$q_end_prev = $extract_end;
										$q_start_prev = $extract_start;
										$q_len_prev = $q_len_OK;
										$go_on = "no";
									}
									elsif($go_on eq "yes"&&$extract_edge =~ /N+/i){
										print RESULTS "$genome\t", "$gene\t", "N\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
										$q_end_prev = $extract_end;
										$q_start_prev = $extract_start;
										$q_len_prev = $q_len_OK;
										$go_on = "no";
									}
									elsif($go_on eq "yes"){
										if($extra_codon =~ /[RYSWKMBDHVNX]/i){
											print RESULTS "$genome\t", "$gene\t", "N\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
											$q_end_prev = $extract_end;
											$q_start_prev = $extract_start;
											$q_len_prev = $q_len_OK;
											$go_on = "no";
										}
										else{
											open(NTDB, ">>nt.$gene.fas");
											print NTDB "\n>", "$gene", "_", "$new_allele_len_codon", "_","$allele_new\n", "$extract_codon";
											system ("makeblastdb -in nt.$gene.fas -dbtype nucl -logfile makeblastdb.log");
											print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
											if($strand_OK==0){
												$q_end_prev = $extract_end;
												$q_start_prev = $extract_start - 3;
												$q_len_prev = $q_len_OK;
											}
											else{
												$q_end_prev = $extract_end +3;
												$q_start_prev = $extract_start;
												$q_len_prev = $q_len_OK;
											}
											close NTDB;
											$go_on = "no";
										}
									}
									else{}
								}
								else{}
							}
							else{}
						}
						else{
							print RESULTS "$genome\t", "$gene\t", "M\n";
							$go_on = "no";
						}
					}
					else{
						print RESULTS "$genome\t", "$gene\t", "M\n";
						$go_on = "no";
					}
				}
			}
			else{}
		}
		else{
			my($line_counter, $copy_x, $extract_start, $extract_end, $new_allele_len, $seqin, $extract,
				$extract_edge_start, $extract_edge_end, $extract_edge_len, $extract_edge, $dis_index);
			my(@blastx, @blastx1, @distances, @sorted, @pozs, @pozs_sorted, @lines, @lines_1st, @dises_index); 
			my(%poz_hash_1st, %poz_edge_hash_1st, %distance_hash_1st, %allele_hash_1st, %strand_hash_1st, %coverage_hash_1st, %identity_hash_1st,
				%poz_hash, %poz_edge_hash, %distance_hash, %allele_hash, %strand_hash, %coverage_hash, %identity_hash);
			print ("blastx: ", "$genome", " vs. ", "$gene", " ...\n") if ($opt_v == 0);
			print BLAST "blastx: $genome vs. $gene\n", "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n";
			@blastx = readpipe("blastx -query $genome.combined.fas -db aa.$gene.fas -gapopen 11 -seg no -outfmt 6");
			$copy_x=0;
			if(@blastx){
				foreach (@blastx){
					print BLAST "$_";
					chomp;
					$line_counter++;
					@blastx1 = split (/\t/,$_);
					$allele_len = $blastx1[1];
					$allele_len =~ s/.*_//g;
					
					$identity = $blastx1[2];
					$aln_len = $blastx1[3];
					$mismatch = $blastx1[4];
					$gap = $blastx1[5];
					$q_start = $blastx1[6];
					$q_end = $blastx1[7];
					$s_start = $blastx1[8];
					$s_end = $blastx1[9];
					
					if($q_end>$q_start){
						$strand = "plus";
						}
					else{
						$strand = "minus";
						$q_end2 = $q_end;
						$q_end = $q_start;
						$q_start = $q_end2;
					}
					
					$s_len = abs($s_end-$s_start)+1;
					$q_len = abs($q_end-$q_start)+1;
					
					$coverage = $s_len/$allele_len;
					$coverage = $coverage*100;
					
					$allele_len = $allele_len*3 + 3;
					
					unless($s_start==1){
						$s_start = 3*($s_start-1)+1;
					}
					$s_end = 3*$s_end;
					
					if($strand eq "minus"){
						$edge1_start = $q_start - $allele_len + $s_end;
						$edge2_end = $q_end + $s_start - 1;
					}
					else{
						$edge1_start = $q_start - $s_start + 1;
						$edge2_end = $q_end + $allele_len - $s_end;
					}
					
					$poz = "$q_start"."\t"."$q_end";
					$poz_edge = "$edge1_start"."\t"."$edge2_end";
					
					if($identity>=$identity_i&&$coverage>=$coverage_t){
						$copy_x++;
						if($coverage>=$coverage_c){
							$allele = $allele_new;
						}
						else{
							$allele = "T";
						}
						
						@pozs=(); @pozs_sorted=();
						push @pozs, $q_start;
						push @pozs, $q_end;
						push @pozs, $q_start_prev;
						push @pozs, $q_end_prev;
						@pozs_sorted = sort{$a<=>$b}@pozs;
						$distance = $pozs_sorted[2] - $pozs_sorted[1];
						$total_len = $pozs_sorted[3] - $pozs_sorted[0] + 1;
						$coding_len = $q_len + $q_len_prev;
						unless($total_len>$coding_len){
							$distance = 0;
						}
						$dis_index = abs($distance);
						
						push @dises_index, $dis_index;
						$poz_hash{$dis_index} = $poz;
						$poz_edge_hash{$dis_index} = $poz_edge;
						$allele_hash{$dis_index} = $allele;
						$strand_hash{$dis_index} = $strand;
						$coverage_hash{$dis_index} = $coverage;
						$identity_hash{$dis_index} = $identity;
						$distance_hash{$dis_index} = $distance;
						$q_len_hash{$dis_index} = $q_len;
					}
					else{}
				}
				
				if(@dises_index){
					@sorted=sort{$a<=>$b}@dises_index;
					$distance_OK = $sorted[0];
					$poz_OK = $poz_hash{$distance_OK};
					$poz_edge_OK = $poz_edge_hash{$distance_OK};
					$allele_OK = $allele_hash{$distance_OK};
					$coverage_OK = $coverage_hash{$distance_OK};
					$identity_OK = $identity_hash{$distance_OK};
					$strand_OK = $strand_hash{$distance_OK};
					$distance_OK = $distance_hash{$distance_OK};
					$q_len_OK = $q_len_hash{$distance_OK};
					
					my($extract_start, $extract_end, $new_allele_len, $seqin, $extract,
							$extract_edge_start, $extract_edge_end, $extract_edge_len, $extract_edge,
							$extract_start_codon, $new_allele_len_codon, $extract_codon, $extract_codon_len, $codon_start, $extra_codon);
					open(SEQIN, "<$genome.combined.fas");
					$extract_start = $poz_OK;
					$extract_start =~ s/\t.*//;
					$extract_end = $poz_OK;
					$extract_end =~ s/.*\t//;
					$new_allele_len = $extract_end-$extract_start+1;
					$extract_start = $extract_start-1;
					
					$extract_edge_start = $poz_edge_OK;
					$extract_edge_start =~ s/\t.*//;
					$extract_edge_end = $poz_edge_OK;
					$extract_edge_end =~ s/.*\t//;
					$extract_edge_len = $extract_edge_end-$extract_edge_start + 1;
					$extract_edge_start = $extract_edge_start-1;
					
					$/=undef;
					$seqin = <SEQIN>;
					$seqin =~ s/>.*\n//;
					$seqin =~ s/\s//g;
					$extract = substr $seqin, $extract_start, $new_allele_len;
					$extract_edge = substr $seqin, $extract_edge_start, $extract_edge_len;
					$/="\n";
					
					if($strand_OK eq "minus"){
						$extract = reverse($extract);
						$extract =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
						$extract_edge = reverse($extract_edge);
						$extract_edge =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
						
						$extract_start_codon = $extract_start - 3;
						$new_allele_len_codon = $new_allele_len + 3;
						$extract_codon = substr $seqin, $extract_start_codon, $new_allele_len_codon;
						$extract_codon = reverse($extract_codon);
						$extract_codon =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
					}
					else{
						$extract_start_codon = $extract_start;
						$new_allele_len_codon = $new_allele_len + 3;
						$extract_codon = substr $seqin, $extract_start_codon, $new_allele_len_codon;
					}
					
					$extract_start = $extract_start+1;
					$extract =~ s/\s//g;
					$extract_edge =~ s/\s//g;
					$extract_codon =~ s/\s//g;
					
					$extract_codon_len = length($extract_codon);
					$codon_start = $extract_codon_len - 3;
					$extra_codon = substr $extract_codon, $codon_start, 3;
					
					if($go_on eq "yes"&&$perfectn_found eq "yes"&&$copy_x==1){
						print RESULTS "$genome\t", "$gene\t", "$perfectn_hit\t", "$q_start_prev\t", "$q_end_prev\t", "$perfectn_poz\t", "$d\t", "$perfectn_d\t", "blastn\t", "$copy\n";
						$q_end_prev = $perfectn_poz;
						$q_end_prev =~ s/.*\t//;
						$q_start_prev = $perfectn_poz;
						$q_start_prev =~ s/\t.*//;
						$q_len_prev = $q_len_OK;
						$go_on = "no";
					}
					elsif($go_on eq "yes"&&$copy_x>1&&$distance_OK>$d){
						print RESULTS "$genome\t", "$gene\t", "D\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
						$go_on = "no";
						$q_end_prev = "NA";
						$q_start_prev = "NA";
						$q_len_prev = $q_len_OK;
					}
					elsif($go_on eq "yes"){
						if($allele_OK eq "T"){
							print RESULTS "$genome\t", "$gene\t", "T\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
							$q_end_prev = $extract_end;
							$q_start_prev = $extract_start;
							$q_len_prev = $q_len_OK;
							$go_on = "no";
						}
						elsif($allele_OK==$allele_new){
							if($go_on eq "yes"&&$extract =~ /[RYSWKMBDHVNX]/i){
								print RESULTS "$genome\t", "$gene\t", "n\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
								$q_end_prev = $extract_end;
								$q_start_prev = $extract_start;
								$q_len_prev = $q_len_OK;
								$go_on = "no";
							}
							elsif($go_on eq "yes"&&$extract_edge =~ /N+/i){
								print RESULTS "$genome\t", "$gene\t", "N\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
								$q_end_prev = $extract_end;
								$q_start_prev = $extract_start;
								$q_len_prev = $q_len_OK;
								$go_on = "no";
							}
							elsif($go_on eq "yes"){
								if($extra_codon =~ /[RYSWKMBDHVNX]/i){
									print RESULTS "$genome\t", "$gene\t", "N\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
									$q_end_prev = $extract_end;
									$q_start_prev = $extract_start;
									$q_len_prev = $q_len_OK;
									$go_on = "no";
								}
								else{
									open(NTDB, ">>nt.$gene.fas");
									print NTDB "\n>", "$gene", "_", "$new_allele_len_codon", "_","$allele_new\n", "$extract_codon";
									system ("makeblastdb -in nt.$gene.fas -dbtype nucl -logfile makeblastdb.log");
									print RESULTS "$genome\t", "$gene\t", "$allele_OK\t", "$q_start_prev\t", "$q_end_prev\t", "$poz_OK\t", "$d\t", "$distance_OK\t", "blastx\t", "$copy_x\n";
									if($strand_OK==0){
										$q_end_prev = $extract_end;
										$q_start_prev = $extract_start - 3;
										$q_len_prev = $q_len_OK;
									}
									else{
										$q_end_prev = $extract_end +3;
										$q_start_prev = $extract_start;
										$q_len_prev = $q_len_OK;
									}
									close NTDB;
									$go_on = "no";
								}
							}
							else{}
						}
						else{}
					}
					else{}
				}
				else{
					print RESULTS "$genome\t", "$gene\t", "M\n";
					$go_on = "no";
				}
			}
			else{
				print RESULTS "$genome\t", "$gene\t", "M\n";
				$go_on = "no";
			}
		}
	}
	close GENE;
}

close GENOME;
close RESULTS;
close BLAST;
system ("rm -f makeblastdb.log");
###################
my ($inline, $in, $genome, $genes, $gene, $genome_previous, $allele, $coordinates, $copy, $line_counter, $gene_list,
		 $first_col, $second_col, $number_of_column, $allele_sum1, $allele_sum2, $duplicates);
my (@in, @gene_list);

open(IN, "<gene_list.tmp");
open(OUT, ">gene_list2.tmp");
while(<IN>){
	chomp;
	$_ =~ s/\t.*//;
	print OUT "$_\n";
}
close IN; close OUT;
system("rm -f gene_list.tmp");
system("mv gene_list2.tmp gene_list.tmp");

open(IN, "<output.txt");
open(OUT, ">allele_matrix.tmp");
open(GENE, "<gene_list.tmp");
@gene_list = <GENE>;

$genes = join ("", @gene_list);
$genes =~ s/\n/\t/gs;

$line_counter = 0;

while(<IN>){
	chomp;
	$inline = $_;
	@in = split (/\t/,$inline);
	
	$line_counter++;
	if($line_counter == 1){
		print OUT "genome\t", "$genes";
	}
	else{
		$genome_previous = $genome;
		
		$genome = $in[0];
		$gene = $in[1];
		$allele = $in[2];
		
		if("$genome" eq "$genome_previous"){
			print OUT "\t", "$allele";
			}
		else{
			print OUT "\n", "$genome", "\t", "$allele";
		}
	}
}
close IN;
close OUT;
###################
my ($inline, $col_sum, $col);
my (@column);

open(IN, "<allele_matrix.tmp");
open(OUT, ">>allele_matrix.transposed.tmp");

while (<IN>){
	if($_ =~ m/^\s/){
	next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' allele_matrix.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
close IN;
close OUT;
###################
open(IN, "<allele_matrix.transposed.tmp");
open(OUT, ">allele_matrix.sorted.tmp");
open(OUT2, ">allele_matrix.cleaned.tmp");
open(REPORT, ">>output.report.tmp");
my ($inline, $in, $genome, $gene, $gene_previous, $allele, $coordinates, $copy, $counter,
		 $first_col, $second_col, $number_of_column, $allele_sum1, $allele_sum2, $duplicates, $shared);
my (@in, @in2);

my $line_counter = 0;
my $na = 0;
my $dup = 0;
my $same = 0;
my $dif = 0;
while (<IN>){
	$line_counter++;
	chomp;
	if($line_counter == 1){
		next;
	}
	else{
		$inline = $_;
		$inline =~ s/^\s//;
		$genome = $inline;
		$genome =~ s/\t.*//g;
		$inline =~ s/$genome\t//;
		$inline =~ s/\t$//;
		@in = split (/\t/,$inline);
		$number_of_column = @in;
		$allele_sum1 = $in[0] * $number_of_column;
		
		if($inline =~ /^\s/){
			next;
		}
		else{
			if($inline =~ /[mnt]/i){
				$na++;
				print OUT "$inline\t", "4NA", "\n";
			}
			elsif($inline =~ /d/i){
				$dup++;
				print OUT "$inline\t", "2dup", "\n";
			}
			else{
				my($element1, $element2, $judgement);
				$judgement = "same";
				foreach(@in){
					$element1 = $_;
					foreach(@in){
						$element2 = $_;
						if($element1 != $element2){
							$judgement = "dif";
						}
						else{
							next;
						}
					}
				}
				if($judgement eq "same"){
					$same++;
					print OUT "$inline\t", "3same", "\n";
					print OUT2 "$inline\n";
				}
				else{
					$dif++;
					print OUT "$inline\t", "1dif", "\n";
					print OUT2 "$inline\n";
				}
			}
		}
	}
}
$shared = $same + $dif;
system("paste gene_list.tmp allele_matrix.sorted.tmp >allele_matrix.sorted2.tmp");
print REPORT "Number of shared-loci: ", "$shared<br>",
						"Number of shared-loci that are identical: ", "$same<br>",
						"Number of shared-loci that are polymorphic: ", "$dif<br>",
						"Number of loci that is excluded because of gene duplication: ", "$dup<br>",
						"Number of loci that is excluded because of incomplete information (missing, truncation or containing nucleotide ambiguity): ", "$na<br>";
close IN;
close OUT;
close OUT2;
close REPORT;
###################
open(IN, "<allele_matrix.sorted2.tmp");
open(OUT, ">>allele_matrix.sorted.transposed.tmp");
my ($inline, $in, $genome, $gene, $gene_previous, $allele, $coordinates, $copy, $counter,
		 $first_col, $second_col, $number_of_column, $allele_sum1, $allele_sum2, $duplicates, $shared);
my (@column, @genome);

while (<IN>){
	if($_ =~ m/^\s/){
	next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' allele_matrix.sorted2.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
open(LIST, "<genome_list.tmp");
open(LIST2, ">genome_list2.tmp");
print LIST2 "\n";
while(<LIST>){
	chomp;
	print LIST2 "$_\n";
}
system("paste genome_list2.tmp allele_matrix.sorted.transposed.tmp >output.allele_matrix.txt");
close IN;
close OUT;
close LIST;
close LIST2;
###################
open(IN, "<output.allele_matrix.txt");
my ($isolate_key, $profile_value, $gene_list, $type_list, $gene_key, $type_value, $key, $counter,
		$gene, $gene_counter, $type, $allele1, $allele2, $isolate1, $isolate2, $profiles1, $profiles2);
my (@infile, @genes, @isolates, @types, @profile1, @profile2);
my (%gene_hash, %isolate_hash);

while(<IN>){
	chomp;
	if($_ =~ /^\t/){
		next;
	}
	else{
		$_ =~ s/\t$//;
		$isolate_key = $_;
		$isolate_key =~ s/\t.*//g;
		push @isolates, $isolate_key;
		$profile_value =$_;
		$profile_value =~ s/$isolate_key\t//;
				$profile_value =~ s/\n//;
		$isolate_hash{$isolate_key} = $profile_value;
	}
}

open(IN, "<output.allele_matrix.txt");
@infile = <IN>;
$gene_list = $infile[0];
$gene_list =~ s/^\t//;
$gene_list =~ s/\t\n//;
$type_list = $infile[-1];
$type_list =~ s/^\t//;
$type_list =~ s/\t\n//;
@genes = split('\t', $gene_list);
@types = split('\t', $type_list);

$counter = 0;
foreach(@genes){
	$gene_key = $_;
	$type_value = @types["$counter"];
	$gene_hash{$gene_key} = $type_value;
	$counter++;
}

foreach(@isolates){
	$isolate1 = $_;
	$profiles1 = $isolate_hash{"$isolate1"};
	$gene_counter = 0;
	
	foreach(@genes){
		$gene = $_;
		$type = $types["$gene_counter"];
		
		$gene_counter++;
		if($type eq "1dif"){
			$gene_counter = $gene_counter - 1;
			@profile1 = split ('\t', $profiles1);
			$allele1 = $profile1["$gene_counter"];
			
			foreach(@isolates){
				$isolate2 = $_;
				$profiles2 = $isolate_hash{"$isolate2"};
				@profile2 = split ('\t', $profiles2);
				$allele2 = $profile2["$gene_counter"];
				
				if($allele1 == $allele2){
					next;
				}
				else{
					open (PROFILE, ">>pairwise.$isolate1.$isolate2.txt");
					print PROFILE "$gene\t", "$allele1\t", "$allele2\n";
					close PROFILE;
				}
			}
			$gene_counter++;
		}
		else{
			next;
		}
	}
}
###################
my ($inline, $col_sum, $col);
my (@column);

open(IN, "<allele_matrix.cleaned.tmp");
open(OUT, ">>allele_matrix.cleaned.transposed.tmp");

while (<IN>){
	if($_ =~ m/^\s/){
	next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' allele_matrix.cleaned.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
close IN;
close OUT;
###################
open(IN_1, "<allele_matrix.cleaned.transposed.tmp");
open(GENOME, "<genome_list.tmp");
open(OUT, ">>dif.tmp");
open(OUT2, ">>dis.tmp");
open(BAPS, ">output.BAPS.txt");
open(STRU, ">Structure.tmp");

my ($isolate_1, $isolate_2, $allele_counter, $dif_counter, $allele1, $allele2, $distance,
		$inline, $outline, $genome_name, $line_counter, $genome_counter, $line_counter, $in);
my (@allele_1, @allele_2, @inline);

$line_counter = 0;
foreach (<IN_1>){
	chomp;
	$line_counter++;
	$in = $_;
	$in =~ s/\t$//;
	print BAPS "$in\t", "$line_counter\n";
	print STRU "$in\n";
	
	$isolate_1 = $_;
	@allele_1 = split ("\t",$isolate_1);
	
	open(IN_2, "<allele_matrix.cleaned.transposed.tmp");
	foreach(<IN_2>){
		chomp;
		$allele_counter = 0;
		$dif_counter = 0;
		
		$isolate_2 = $_;
		@allele_2 = split ("\t",$isolate_2);
		
		foreach(@allele_2){
			$allele2 = $_;
			$allele1 = @allele_1[$allele_counter];
			$allele_counter++;

			if ($allele1!=$allele2){
			$dif_counter++;
			}
			else{next;
			}
		}
		$distance = $dif_counter/$allele_counter;
		print OUT "$dif_counter\t";
		print OUT2 "$distance\t";
		close IN_2;
	}
		print OUT "\n";
		print OUT2 "\n";
}

close IN_1;
close OUT;
close BAPS;
close STRU;

system("paste genome_list.tmp Structure.tmp >output.Structure.txt");
system("paste genome_list.tmp dis.tmp >dis.txt.tmp");
system("paste genome_list.tmp dif.tmp >dif.txt.tmp");

open DIF, "<dif.txt.tmp";
open DIS, "<dis.txt.tmp";
open DIF2, ">output.difference_matrix.txt";
open DIS2, ">output.Splitstree.nex";

open (GENOME, "<genome_list.tmp");
$genome_counter = 0;
foreach(<GENOME>){
	$genome_counter++;
}
close GENOME;

print DIS2 
"#NEXUS\n\n",
"BEGIN taxa;", "\n",
"   DIMENSIONS ntax=$genome_counter;\n",
"TAXLABELS\n";

print DIF2 "\t";
open GENOME, "<genome_list.tmp";
foreach(<GENOME>){
	if($_ =~ /^\s/){
		next;
	}
	else{
		chomp;
		print DIS2 "   $_\n";
		print DIF2 "$_\t";
	}
}

print DIF2 "\n";
print DIS2
";\n",
"END;\n\n",
"BEGIN distances;\n",
"   DIMENSIONS ntax=", "$genome_counter;\n",
"   FORMAT\n",
"      triangle=LOWER\n",
"      diagonal\n",
"      labels\n",
"      missing=?\n",
"   ;\n",
"MATRIX\n";

$line_counter = 1;
foreach (<DIS>){
	chomp;
	$line_counter++;
	$inline = $_;
	@inline = split ("\t",$inline);
	$genome_name = $inline[0];
	
	splice @inline, "$line_counter";
	$outline = join("\t", @inline);
	print DIS2 "$outline\n";
#	print DIS2 "@inline\n";
}
print DIS2
";\n",
"END;\n";

$line_counter = 1;
foreach (<DIF>){
	chomp;
	$line_counter++;
	$inline = $_;
	@inline = split ("\t",$inline);
	$genome_name = $inline[0];
	
	splice @inline, "$line_counter";
	$outline = join("\t", @inline);
	print DIF2 "$outline\n";
}
close DIF;
close DIS;
close DIF2;
close DIS2;
###################
open (IN, "<output.difference_matrix.txt");
open (HTML,'>>difference_matrix.tmp');
my ($isolate1, $isolate2, $list_isolate, $line_counter, $inline, $x_counter, $value, $x, $y, $isolates_number, $counter, $start);
my (@infile, @isolates, @values);

@infile = <IN>;
$list_isolate = $infile[0];
$list_isolate =~ s/^\t//;
$list_isolate =~ s/\n//;
@isolates = split ("\t", $list_isolate);

print HTML <<END_OF_HTML;
<HTML><HEAD><TITLE>Difference matrix</TITLE></HEAD><BODY><div><H2><a name="top">Difference matrix</a></H2></div><TABLE BORDER=1><TR>
END_OF_HTML

print HTML "<td><br></td>";

foreach(@isolates){
	print HTML "<td>$_</td>";
}

print HTML "</TR>";

$line_counter = 0;
foreach(@infile){
	chomp;
	if($_ =~ /^\t/){
		next;
	}
	else{
		$inline = $_;
		@values = split ("\t", $inline);
		$y = shift @values;
		
		$line_counter++;
		print HTML "<TR>", "<td>$y</td>";
		
		$x_counter = 0;
		foreach(@values){
			$value = $_;
			$x = $isolates["$x_counter"];
			$x_counter++;
			print HTML "<TH><a href=\"#$x.$y\">$value</a></TH>"
		}
		
		print HTML "</TR>";
	}
}
print HTML "</TABLE>";


open (GENE, "<output.gene_list.txt");
my ($alias, $gene, $allele1, $allele2);
my (@inline);
my (%gene_hash);

while(<GENE>){
	chomp;
	$alias = $_;
	$alias =~ s/\t.*//;
	$gene = $_;
	$gene =~ s/.*\t//;
	$gene_hash{$alias} = $gene;
}
close GENE;

$isolates_number = @isolates;
$counter = 0;
foreach(@isolates){
	$isolate1 = $_;
	$counter++;
	$start = $counter;
	while($start<$isolates_number){
		$isolate2 = $isolates["$start"];
		$start++;
		open (PAIR, "<pairwise.$isolate1.$isolate2.txt");

if($opt_l == 0){
print HTML <<END_OF_HTML;
<H3><a name="$isolate1.$isolate2">$isolate1 VS. $isolate2</a></H3><TABLE BORDER=1><TR><TH>alias</TH><TH>$isolate1</TH><TH>$isolate2</TH><TH>product</TH></TR>
END_OF_HTML
}
else{
print HTML <<END_OF_HTML;
<H3><a name="$isolate1.$isolate2">$isolate1 vs. $isolate2</a></H3>
END_OF_HTML
}

		foreach(<PAIR>){
			chomp;
			@inline = split ("\t", $_);
			$alias = $inline[0];
			$gene = $gene_hash{$alias};
			$allele1 = $inline[1];
			$allele2 = $inline[2];

if($opt_l == 0){
print HTML <<END_OF_HTML;
<TR><Td><a href="../aln_$T0/$alias.aln.txt">$alias</a></Td><Th>$allele1</Th><Th>$allele2</Th><Td>$gene</Td></TR>
END_OF_HTML
}
else{
	print HTML "$alias: ", "$allele1 ", "$allele2; ", "$gene", "<br>";
}

		}
		close PAIR;

print HTML <<END_OF_HTML;
</TABLE><div><a href="#top">top</a></div>
END_OF_HTML

	}
}

print HTML <<END_OF_HTML;
</BODY></HTML>
END_OF_HTML

close HTML;
###################
open (GENE, "<gene_list.tmp");
my ($gene, $seqfile, $allele_counter);
while(<GENE>){
	chomp;
	$gene = $_;
	$seqfile = "nt."."$gene".".fas";
	open (SEQ, "<$seqfile");
	open (OUT, ">$gene.allele.fas");
	$allele_counter = 0;
	while(<SEQ>){
		chomp;
		if($_ =~ />/){
			$_ =~ s/_[0-9]+_/_/;
			$allele_counter ++;
			print OUT "$_\n";
		}
		else{
			print OUT "$_\n";
		}
	}
	close OUT;
	close SEQ;
	
	if($allele_counter>1){
		system ("mafft --clustalout --quiet $gene.allele.fas >$gene.aln.txt");
		system ("mafft --quiet $gene.allele.fas >$gene.aln.fas") if ($opt_a==0);
		system ("rm -f $gene.allele.fas");
	}
	else{
		system ("mv $gene.allele.fas $gene.aln.fas");
	}
}
close GENE;
###################
unless($opt_a==1){
	my($inline, $line_counter, $gene, $seqname, $allele, $allele_db, $seq, $isolate);
	my(@isolates, @isolates2, @alleles);
	open(IN, "<allele_matrix.transposed.tmp");
	open(OUT, ">allele_profiles.tmp");
	$line_counter = 0;
	while(<IN>){
		chomp;
		$inline = $_;
		$line_counter++;
		
		if($line_counter == 1){
			print OUT "$inline\n";
		}
		elsif($inline=~ /\t[MNTD]\t/i){
			next;
		}
		else{
			print OUT "$inline\n";
		}
	}
	close IN;
	close OUT;
	
	open(IN, "<allele_profiles.tmp");
	open(OUT, ">clonalframe.dat");
	$line_counter = 0;
	while(<IN>){
		$inline = $_;
		$inline =~ s/\t\n//;
		$line_counter++;
		
		if($line_counter == 1){
			@isolates = split ("\t", $inline);
			shift @isolates;
			next;
		}
		elsif($inline=~/\t\t/){
			next;
		}
		else{
			@isolates2 = @isolates;
			@alleles = split ("\t", $inline);
			$gene = shift @alleles;
			print OUT "#", "$gene\n";
			
			foreach(@alleles){
				$isolate = shift @isolates2;
				open(ISOLATE, ">>$isolate.core.fas.tmp");
				$allele = $_;
				
				open (DB, "<$gene.aln.fas");
				$/=">";
				
				foreach(<DB>){
					$seqname = $_;
					$seqname =~ s/\n.*//gs;
					$allele_db = $seqname;
					$allele_db =~ s/.*_//g;
					$allele_db =~ s/new//;
					$seq = $_;
					$seq =~ s/$seqname//;
					$seq =~ s/\n//g;
					$seq =~ s/>//;
					
					if($allele_db == $allele){
						print OUT ">", "$isolate\n", "$seq\n";
						print ISOLATE "$seq\n";
					}
					else{
						next;
					}
				}
				
				close DB;
				close ISOLATE;
				$/="\n";
			}
			
		}
		print OUT "=", "\n";
	}
	close OUT;
	open(OUT, ">core_genomes.fas");
	foreach(@isolates){
		$isolate = $_;
		$/="";
		open(IN, "<$isolate.core.fas.tmp");
		foreach(<IN>){
			print OUT ">$isolate.core\n$_";
			system("rm -f $isolate.core.fas.tmp");
		}
		close IN;
	}
	close OUT;
}
###################
open(REPORT, ">>output.report.tmp");
my $T1 = time();
my $time = $T1-$T0;
my $hour = $time/3600;
my $hour_rounded = sprintf("%.1f", $hour);
print REPORT "Start time: $T0<br>", "Total running time: $hour_rounded hours ($time seconds)"; 
close REPORT;

system ("mkdir tmp_$T0");
system ("mkdir combined_contigs_$T0");
system ("mkdir aln_$T0");
system ("mkdir output_$T0");

if($opt_b == 0){
	system ("mv output.BAPS.txt ./output_$T0/BAPS.txt");
}
else{
	system ("rm -f output.BAPS.txt");
}

unless($opt_a == 1){
	system ("mv core_genomes.fas ./output_$T0/");
	system ("mv clonalframe.dat ./output_$T0/");
	system ("rm -f gene*.aln.fas");
}
system ("mv output.allele_matrix.txt ./tmp_$T0/allele_matrix.txt");
system ("mv output.difference_matrix.txt ./tmp_$T0/difference_matrix.txt");
system ("cat output.report.tmp difference_matrix.tmp >difference_matrix.html");
system ("mv output.Splitstree.nex ./output_$T0/Splitstree.nex");
system ("mv output.Structure.txt ./output_$T0/allele_profiles.txt");
system ("cp output.gene_list.txt ./scheme/output.gene_list.txt");
system ("mv output.gene_list.txt ./output_$T0/gene_list.txt");
system ("mv output.txt ./output_$T0/");
system ("mv difference_matrix.html ./output_$T0/");

system ("mv *.combined.fas ./combined_contigs_$T0/");
system ("mv *.tmp ./tmp_$T0/");
system ("mv *.aln.* ./aln_$T0/");
system ("rm -f *.*hr");
system ("rm -f *.*in");
system ("rm -f *.*sq");

#system ("rm -f pairwise.*.txt");
my @files = readpipe("ls");
foreach(<@files>){
	chomp;
	my $filename = $_;
	if($filename =~ m/pairwise\./){
		system "rm $filename";
	}
}


system ("rm -f *.allele.fas");
system ("mv aa.* ./scheme");
system ("mv nt.* ./scheme");
system ("mv scheme_summary.txt ./scheme");

unless($opt_m == 1){
	system ("rm -rf tmp_$T0");
}



