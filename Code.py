# Combine fastq files
def combine_fastq(input_files, output_file):
    with open(output_file, 'w') as outfile:
        for input_file in input_files:
            with open(input_file, 'r') as infile:
                outfile.write(infile.read())
# Different libraries of the same read (R1 or R2)
input_files = ['L1.fastq', 'L2.fastq', 'L3.fastq']
output_file = 'combined R1.fastq'
combine_fastq(input_files, output_file)


#Set paired reads

from Bio import SeqIO
def pair_reads(r1_file, r2_file, output_file):
    # Open the output file for writing
    with open(output_file, 'w') as out_handle:
        # Open the R1 and R2 files for reading
        with open(r1_file, 'r') as r1_handle, open(r2_file, 'r') as r2_handle:
            # Iterate through paired records from R1 and R2
            for record_r1, record_r2 in zip(SeqIO.parse(r1_handle, 'fastq'), SeqIO.parse(r2_handle, 'fastq')):
                # Check that the read identifiers match
                if record_r1.id.split()[0] == record_r2.id.split()[0]:
                    # Write the paired records to the output file
                    SeqIO.write([record_r1, record_r2], out_handle, 'fastq')
                else:
                    print(f"Warning: Read identifiers do not match for {record_r1.id} and {record_r2.id}")

r1_file = 'combined R1.fastq'
r2_file = 'combined R2.fastq'
output_file = 'paired_reads.fastq'
pair_reads(r1_file, r2_file, output_file)

#Quality filtering

from Bio import SeqIO
count = 0
for rec in SeqIO.parse("paired_reads.fastq", "fastq"):
    count += 1
print("%i reads" % count)

good_reads = (
    rec
    for rec in SeqIO.parse("paired_reads.fastq", "fastq")
    if min(rec.letter_annotations["phred_quality"]) >= 20
)
count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)
