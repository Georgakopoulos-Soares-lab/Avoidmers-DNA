import re,os,sys,glob
from Bio import SeqIO

input_file="chm13v2.0.fa"

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
length=8
bin_size=50_000
datafile=open("number_of_unique_kmers_per_bin_hg38_binsize_"+str(bin_size)+"_kmer_length_"+str(length)+".txt","w")

for fasta in fasta_sequences:
	name, fasta_sequence = fasta.id, str(fasta.seq)
	for bin in range(0, len(fasta_sequence)-bin_size+1, bin_size):
		kmers=set()
		sequence=fasta_sequence[bin:bin+bin_size]
		for i in range(0, len(sequence)-length+1, length):
			kmer=sequence[i:i+length]
			kmers.add(kmer)
		#print(len(kmers),name,bin)
		datafile.write(name+'\t'+str(bin)+'\t'+str(bin+bin_size)+'\t'+str(len(kmers))+'\n')
datafile.close()

