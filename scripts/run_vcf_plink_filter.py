import os
import subprocess

# declare path to sub-directories/folder
HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, os.pardir, 'input')
OUTPUT = os.path.join(HERE, os.pardir, 'output')

# Open files for stdout and stderr
out = open('path.txt', 'w')
err = open('path.err', 'w')

# create output directory if it doesn't exist already
if not os.path.isdir(OUTPUT):
	os.makedirs(OUTPUT)


def run_vcftools_filter(file, MAX_ALLELE, MIN_ALLELE, MISS, MAF, QUAL):
	"""

	:param file: compressed input vcf file
	:param MAX_ALLELE: includes sites with number of alleles less than or equal to two ( bi-allelic)
	:param MIN_ALLELE: includes sites with number of alleles greater than or equal to two ( bi-allelic)
	:param MISS: minimum missing data threshold
	:param MAF: set minor allele frequency threshold
	:param QUAL: includes sites with quality value above 30 threshold
	:return: filtered vcf.gz file
	"""

	# create output directory if it doesn't exist
	if not os.path.isdir(os.path.join(OUTPUT, 'filter')):
		os.makedirs(os.path.join(OUTPUT, 'filter'))

	command = ' '.join(
		['vcftools --gzvcf', file[0], '--max-alleles', MAX_ALLELE, '--min-alleles', MIN_ALLELE, '--max-missing', MISS,
		 '--minQ', QUAL, '--maf', MAF,
		 '--recode', '--recode-INFO-all', '--stdout', '|', 'bgzip', '>',
		 os.path.join(OUTPUT, 'filter', 'acs_mini_project.filtered.vcf.gz')])
	subprocess.run(command, stdout=out, stderr=err, shell=True)


def export_vcf_to_plink(vcf_filter, plink):
	"""

	:param vcf_filter: compressed input vcf filtered file
	:param plink: path to the plink exe
	:return: output three files (.pgen, .pvar, .psam). .pgen contains binary genomic data file. .pvar contains variant information with a header line
	and last header line contains seven columns (Chromosome, variant ID, reference allele, alternate alleles, phred-scaled quaity score, filter).
	.psam contains sample information file individual ID, sex
	"""
	command = ' '.join(
		[plink, '--vcf', vcf_filter, '--out', os.path.join(OUTPUT, 'filter', 'acs_mini_project.filtered')])
	subprocess.run(command, stdout=out, stderr=err, shell=True)


def calculate_LD_R2_threshold(PGEN, PVAR, PSAM, plink, R2):
	"""

	:param PGEN: binary genomic data file
	:param PVAR: variant information
	:param PSAM: sample information
	:param plink: path to the plink exe
	:param R2: pairwise threshold
	:return: will prouduce two files (purne.in and purne.out). All variants within a window of 1kb and r2 greater than 0.1 are removed and
	placed in the purne.out file.
	"""
	command = ' '.join(
		[plink, '--pgen', PGEN, '--pvar', PVAR, '--psam', PSAM, '--set-missing-var-ids @:#', '--indep-pairwise 1 kb 1',
		 R2, '--out', os.path.join(OUTPUT, 'filter', 'acs_mini_project.filtered.r2')])
	subprocess.run(command, stdout=out, stderr=err, shell=True)


def prune_ld_run_pca(PGEN, PVAR, PSAM, plink, PRUNE_OUT):
	"""

	:param PGEN: binary genomic data file
	:param PVAR: variant information
	:param PSAM: sample information
	:param plink: path to the plink exe
	:param PRUNE_OUT:  excluded list variants
	:return:
	"""
	command = ' '.join(
		[plink, '--pgen', PGEN, '--pvar', PVAR, '--psam', PSAM, '--exclude', PRUNE_OUT, '--variance-standardize',
		 '--pca', '--make-bpgen',
		 '--out', os.path.join(OUTPUT, 'filter', 'acs_mini_project.filtered.r2')])
	subprocess.run(command, stdout=out, stderr=err, shell=True)


if __name__ == '__main__':
	vcf_file = [os.path.join(INPUT, files) for files in os.listdir(INPUT) if files.endswith('bgz')][0]
	MAX_ALLELE = str(2)
	MIN_ALLELE = str(2)
	MISS = str(0.9)
	MAF = str(0.01)
	QUAL = str(30)
	run_vcftools_filter(vcf_file, MAX_ALLELE, MIN_ALLELE, MISS, MAF, QUAL)

	vcf_filter = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
				  if i.endswith('gz')][0]
	plink = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
			 if i.endswith('plink2')][0]
	export_vcf_to_plink(vcf_filter, plink)
	PSAM = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
			if i.endswith('psam')][0]
	PVAR = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
			if i.endswith('pvar')][0]
	PGEN = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
			if i.endswith('pgen')][0]
	R2 = str(0.1)
	calculate_LD_R2_threshold(PGEN, PVAR, PSAM, plink, R2)
	PRUNE_OUT = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
				 if i.endswith('prune.out')][0]
	prune_ld_run_pca(PGEN, PVAR, PSAM, plink, PRUNE_OUT)
