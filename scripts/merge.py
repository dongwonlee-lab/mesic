from datetime import datetime
import os
import sys
import subprocess
import logging

from check_samples import check_samples

def format_files(vcf, samfn, varlist, outvcf):
    logging.info('####################')
    logging.info('VCF file:\t%s' % vcf)
    logging.info('Sample file:\t%s' % samfn)
    logging.info('Variant file:\t%s' % varlist)
    logging.info('TEMP file:\t%s' % outvcf)
    
    if samfn == None or samfn == "":
        try: 
            subprocess.run(' '.join(['bcftools', 'view', '-Oz', '-o', outvcf, '-i', 'ID=@%s' % varlist, vcf]), shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            sys.exit(1)
            
    else:
        try:
            subprocess.run(' '.join(['bcftools', 'view', '-Oz', '-o', outvcf, '-S', samfn, '--force-samples', '-i', 'ID=@%s' % varlist, vcf]), shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            sys.exit(1)
            
    try:
        subprocess.run(' '.join(['tabix', '-fp', 'vcf', outvcf]), shell=True, check=True)
    except subprocess.CallProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        sys.exit(1)


def run(chrom, filelist_fn, output, tempdir, snpsonly, python_lib):
    start = datetime.now()
    
    ext = output.split('.')[-1]
    if ext != 'gz':
        logging.error('Please double check that your output ends with ".gz". Current file extension is %s' % ext)
        sys.exit(1)
        
    filelists = {}
    count = 0
    f1 = open(filelist_fn, 'r')
    for line in f1.readlines():
        fn = line.strip().split(',')
        vcf = fn[0]
        variant = fn[1]
        sample = fn[2]
        base = os.path.basename(os.path.splitext(vcf)[0])
        temp = '%s/%s_%s.vcf.gz' % (tempdir, base, count)
        filelists[count] = [vcf, variant, sample, temp]
        count = count +1
    f1.close()
    
    logging.info('###########################')
    logging.info('#### Merging VCF files ####')
    logging.info('###########################')
    logging.info('### Arguments')
    logging.info('Chromosome:\t%s' % chrom)
    logging.info('File lists:\t%s' % filelist_fn)
    logging.info('Output:\t%s' % output)
    logging.info('TEMP directory\t%s' % tempdir)
    
    check_samples(filelists, python_lib)
    
    logging.info('####################')
    logging.info('Formatting files...')
    
    for idx in filelists:
        vcf = filelists[idx][0]
        variant = filelists[idx][1]
        sample = filelists[idx][2]
        temp = filelists[idx][3]        
        format_files(vcf, sample, variant, temp)
               
    mergelist = '%s/mergelist-%s.txt' % (tempdir, chrom)
    f1 = open(mergelist, 'w')
    for idx in filelists:
        temp = filelists[idx][3]
        f1.write('%s\n' % temp)    
    f1.close()
           
    logging.info('####################')
    logging.info('Merging cohorts')
    
    if snpsonly:    
        tempmerge = '%s/merged-%s.temp.vcf.gz' % (tempdir, chrom)
        try:
            subprocess.run(' '.join(['bcftools', 'merge', '-m' ,'none', '-Oz', '-o', tempmerge, '-l', mergelist]), shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            sys.exit(1)

        logging.info('####################')
        logging.info('Filtering for SNPs only')
        try:
            f = 'TYPE="snp"'
            subprocess.run(' '.join(['bcftools', 'view', '-Oz', '-o', output, '-i', "'%s'" % f, tempmerge]), shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            sys.exit(1)
    else:
        try:
            subprocess.run(' '.join(['bcftools', 'merge', '-m' ,'none', '-Oz', '-o', output, '-l', mergelist]), shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            sys.exit(1)
    logging.info('####################')
    logging.info('Merging done')
    logging.info('Merged VCF saved to:\n\t%s' % output)
    logging.info('Runtime: ' + str(datetime.now()-start) + '\n')
