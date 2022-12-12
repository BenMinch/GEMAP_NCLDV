import os, sys, re, shlex, subprocess
import pandas as pd
#Filters out the internal standards
#cores = sys.argv[1]
input_dir = sys.argv[1]
outdir = sys.argv[2]
ref_genome= sys.argv[3]

#check if SRR is in the input directory		
for n in os.listdir(input_dir):
	if n.endswith(".fulltrim.gz")or(".fulltrim"):
		forward = os.path.join(input_dir, n)
		reverse= re.sub('_1', '_2',forward)
		#get out of loop if forward equals reverse		
		if forward==reverse:
			print('skip')
		else:
			outpath= os.path.join(outdir, n)
			outpath2= re.sub('gz.fulltrim.gz','.sam', outpath)
			cmd = 'bbmap/bbmap.sh in1='+ forward+ ' in2='+ reverse + ' out='+ outpath2 + ' ref='+ ref_genome + ' -k=14 -minid=0.97 -Xmx30g'
			subprocess.call(cmd, shell=True)
#convert output to bam
print('Converting files')
for i in os.listdir(outdir):
	if i.endswith('.sam'):
		read= os.path.join(outdir,i)
		samtools= 'samtools view -S -b ' + read + ' > ' + re.sub('.sam','.bam',read)
		subprocess.call(samtools, shell=True)

#remove the sams
remove= 'rm '+ outdir+ '/*.sam'
subprocess.call(remove, shell=True)
# sort the bams
print('Sorting bams')
for i in os.listdir(outdir):
	if i.endswith('.bam'):
		read= os.path.join(outdir,i)
		samtools_sort= 'samtools sort '+ read+ ' -o ' +re.sub('.bam', '.bam.sorted', read)
		subprocess.call(samtools_sort, shell=True)
#remove the bams
remove= 'rm '+ outdir+ '/*.bam'
subprocess.call(remove, shell=True)

#Run coverm must have this program installed  (conda environment)
print('Running Coverm')

for i in os.listdir(outdir):
	if i.endswith('.sorted'):
		sample= os.path.join(outdir,i)
		coverm= 'coverm contig --methods rpkm count --bam-files '+ sample + ' > ' + re.sub('.bam.sorted', '.coverm', sample)
		subprocess.call(coverm, shell=True)
remove_sort= 'rm ' + outdir+ '/*.sorted'
subprocess.call(remove_sort, shell=True)

