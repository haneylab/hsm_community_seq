import sys
from subprocess import call


inFile = open(sys.argv[1], 'r')
bc = {}
with open(sys.argv[2], 'r') as o:
    line = o.readline()
    while line:
        if line[0] != ">":
            raise ValueError("Fasta files must start with '>'")
        title = line[1:].rstrip()
        line = o.readline()
        bc[title] = line.rstrip()
        line = o.readline()


header = ''
data = ''


outFile = {}

for i, line in enumerate(inFile):
    if line[0] == "@" and i % 4 == 0:
        if header != '':
            barcode = header.split('+')[-1].split('/')[0]
            if barcode not in bc.values():
                pass
            else:
                if barcode not in outFile:
                    outFile[barcode] = open(barcode + "_R2.fastq", 'a')

                outFile[barcode].write(header + "\n" + data)

        header = line.strip()
        data = ''
    else:
        data += line

barcode = header.split('+')[-1].split('/')[0]
if barcode not in bc.values():
    pass
else:
    if barcode not in outFile:
        outFile[barcode] = open(barcode + "_R2.fastq", 'a')

    outFile[barcode].write(header + "\n" + data)


for samp in bc.keys():
    # call(["mkdir", "{}".format(samp)])
    call(["mv", "{}_R2.fastq".format(bc[samp]), "{}/{}_R2.fastq".format(samp, samp)])
