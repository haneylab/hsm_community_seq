from subprocess import call


samps = [line.rstrip() for line in open("/Users/haneylab/compbio/260319_CommunityTroubleshoot/reads/samples.txt", 'r')]

for sample in samps:
    call(['cutadapt', '-j', '4', '-a', 'AATGATACGGCGACCACCGAGATCTACACGCTNNNNNNNNNNNNTATGGTAATTGTGTGYCAGCMGCCGCGGTAA',
          '-A', 'CAAGCAGAAGACGGCATACGAGATAGTCAGCCAGCCGGACTACNVGGGTWTCTAAT', f'reads/{sample}/{sample}_R1.fastq',
          f'reads/{sample}/{sample}_R2.fastq', '-o', f'trimmed_reads/{sample}_trimmed_R1.fastq', '-p',
          f'trimmed_reads/{sample}_trimmed_R2.fastq'])
