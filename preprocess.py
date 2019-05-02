"""
Set up pipeline for preprocessing samples
Included:
SureCall Agent Pipeline:
	Trimmer
	Align
	Barcode
	*note may include fgbio pipeline as alternative in future
Family Size adjustment
BaseQuality Recalibration

Options
=======
Required
--------
-r	--reference	Reference file
-I	--Input	Input. Depending on which step is selected can have different one.
						bottom directory: process from fastq
						list of files: post process from aligned bam file
						Single file: post process single aligned file
Optional
--------
-s	--start	Where in pipeline to start.
				"fq/fastq" -> start from fastqs
				"bam" -> start from bam
-k	--knownsnps	[Default depends on reference file input] Database files. Multiple inputs allowed
	--subset	[Default finds all options in the fastq_files folder] Only run on a subset of samples within wd 	
	--fastqfolder	[Default fastq_files] location of fastq_files containing folder.
	--familysize [Default: 2] Collapse reads into family size of x
	--keepone	[Default: True] Keep family size of one?
-t	--barcodetag	[Default: XI] tag for barcode in bam files. Agilent uses XI
						*note if implementing fgbio, may change location
	--keeptempfiles [Default: True] Keep temporary files made in between
"""

import sys
import os
from subprocess import call
import re
import argparse
import glob

def fqStart(subset):
"""
determines fq files to use and runs pipeline on them
"""
# find which fq files to use.
	def writelsf():
		with open(sampname + "_fq_preprocess.lsf", "w+") as lsf:
			lsfHeader(lsf, sampname, 24, 6, 64)
			Agilent(lsf, sampname)
			if args.familysize > 1:
				famSize(lsf, sampname, ".umi.sorted.bam")
				BQSR(lsf, ".umi.sorted.debar.bam")
			BQSR(lsf, ".umi.sorted.bam")

if subset is "all":
	files = glob(args.Input[0] + "/" + args.fastqfolder+"/"+"*R1*")
	for f in files:
		sampname = f.split("_R1")[0]
		writelsf()
else:
	with open(args.files[0], 'r+') as files:
		for f in files:
			sampname = f.rstrip['\n']
			writelsf()

def bamStart(runtype):
"""
starts from aligned files and then further processes them. checks to see what needs to be run.
"""
	def writelsf():
		with open(sampname + "_bam_preprocess.lsf", "w+") as lsf:
			lsfHeader(lsf, sampname, 24, 6, 64)
			famSize(lsf, sampname)
			BQSR(lsf, ".debar")
			BQSR(lsf, "")
	if runtype is "single":
		sampname = args.Input[0].re


def Agilent(lsf, sampname):
"""
writes agilent stuff to lsf file
"""
	lsf.write("./agent.sh" + " -d " + args.Input[0] + " -s " + sampname + " -r " + args.reference[0] + " -b " + args.bed[0] + " -fq " + fqfolder)

def BQSR(lsf, sampname, extra):
"""
writes BQSR stuff to lsf file
"""	
	if !isInstance(knownsnps, str):
    	for i in range(0, len(knownsnps)):
        	knownsnps = ",".join(knownsnps)
	lsf.write("./BQSR.sh + -b " + sampname + extra + " -k " + knownsnps + " -o " + sampname + " -r " + args.reference[0])

def famSize(lsf, sampname, extra):
"""
writes family size stuff to lsf file
"""
	lsf.write("./setFamMin.sh -b " + sampname + extra + " -a " + args.barcodetag[0] + " -n " + args.familysize[0])

def lsfHeader(lsf, sampname, time, threads, memory):
	if time < 6:
		length = "short"
	elif time < 24:
		length = "medium"
	else:
		length = "long"
	lsf.write("#BSUB -J " + sampname + '_preprocess.lsf\n')
	lsf.write("#BSUB -W " + time + ":00\n")
	lsf.write("#BSUB -o /rsrch3/home/PanCanRsch/jhuang14/logs" + '\n')
	lsf.write("#BSUB -e /rsrch3/home/PanCanRsch/jhuang14/logs" + '\n')
	lsf.write("#BSUB -cwd /rsrch3/home/PanCanRsch/jhuang14/MolBarcode" + '\n')
	lsf.write("#BSUB -q " + length + '\n')
	lsf.write("#BSUB -u jhuang14@mdanderson.org" + '\n')
	lsf.write("#BSUB -n " + threads + '\n')
	lsf.write("#BSUB -M " + memory + '\n')
	lsf.write("#BSUB -R rusage[mem=" + memory + "]\n")
	lsf.write("#BSUB -N" + '\n')
	lsf.write("#BSUB -B" + '\n')
	lsf.write("module load BWA")
	lsf.write("module load gatk")
	lsf.write("module load samtools")


if __name__ == "__main__":
	## set up argparser
	parser = argparse.ArgumentParser(description='Set up pipeline for preprocessing samples')
	parser.add_argument('-r', '--reference', help="Reference File", nargs=1, required=True, type=str)
	parser.add_argument('-I', '--Input', help="Inputs. Depending on which step is selected can have different one.", nargs=1, required=True, type=str)
	parser.add_argument('-k', '--knownsnps', help="knownsnps", nargs="?", required=False, type=str)
	parser.add_argument('-s', '--start', help="run type", nargs=1, required=False, type=str, default="NA")
	parser.add_argument('-b', '--bed', help="design bed file. Used in Agilent pipeline", nargs=1, required=False, type=str, default="NA")
	parser.add_argument('--subset', help="[Default finds all options in the fastq_files folder] Only run on a subset of samples within wd", nargs="?", required=False, type=str, default="all")
	parser.add_argument('--fastqfolder', help="[Default fastq_files] location of fastq_files containing folder.", nargs=1, required=False, type=str, default="fastq_files")
	parser.add_argument('--familysize', help="[Default: 2] Collapse reads into family size of x", nargs=1, required=False, type=int, default=2)
	parser.add_argument('--keepone', help="[Default: True] Collapse reads into family size of x", nargs=1, required=False, type=bool, default=True)
	parser.add_argument('--barcodetag', help="[Default: 2] Collapse reads into family size of x", nargs=1, required=False, type=str, default="XI")
	parser.add_argument('--keeptempfiles', help="[Default: True] Keep temporary files made in between", nargs=1, required=False, type=bool, default=True)
	args = parser.parse_args()

	## run checks on some inputs

	# check Input type to determine pipeline to run 
	if args.start[0] is "NA":
		if os.path.isfile(args.Input[0]):
			start = "bam"
			if args.Input[0].endswith(".bam"):
				runtype = "single"
			else:
				runtype = "list"
		elif os.path.isdir(args.Input(0)): 
			start = "fq"
		else:
			print("can't find Input")
			exit()
	elif args.start[0] in list("fq", "fastq", "bam"):
		start = args.start[0]
		if args.Input[0].endswith(".bam"):
			runtype = "single"
		else:
			runtype = "list"

	# check knownsnps
	if args.knownsnps is None:
		if "hg19" in args.reference[0]:
			knownsnps = "dbsnp_138.hg19.vcf.gz"
		elif "hg38" in args.reference[0]:
			knownsnps = "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
		else:
			print("Please use hg19 or hg38")
			exit()
	else:
		knownsnps = args.knownsnps

	if start is "fq" or start is "fastq":
		fqStart()
	elif start is "bam":
		bamStart()
