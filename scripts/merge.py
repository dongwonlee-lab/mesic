import subprocess
from datetime import datetime

def run(chrom, afn, bfn, varlist, output,
        asamfn, bsamfn, tempdir):
    start = datetime.now()
    
    print('###########################')
    print('##### Variant Overlap #####')
    print('###########################')
    print('### Arguments')
    print('Chromosome\t%s' % chrom)
    print('VCF A\t\t\t%s' % afn)
    print('VCF B\t\t\t%s' % bfn)
    print('Variant list\t%s' % varlist)
    print('Sample list A\t%s' % asamfn)
    print('Sample list B\t%s' % bsamfn)
    print('Output\t\t\t%s' % output)
    print('TEMP directory\t%s' % tempdir)
    
    amerge = '%s/cohortA.extract.vcf.gz' % tempdir
    bmerge = '%s/cohortB.extract.vcf.gz' % tempdir
    
    print('####################')
    print('Formatting files...')
    print('\tCohort A: %s' % afn)
    if asamfn == None:
        subprocess.run(['bcftools', 'view', '-Oz', '-o',  amerge, '-i', 'ID=@%s' % varlist, afn])
    else:
        subprocess.run(['bcftools', 'view', '-Oz', '-o',  amerge, '-S', asamfn, '-i', 'ID=@%s' % varlist, afn])
    subprocess.run(['tabix', '-p', 'vcf', amerge])
        
    print('\tCohort B: %s' % bfn)
    if bsamfn == None:
        subprocess.run(['bcftools', 'view', '-Oz', '-o',  bmerge, '-i', 'ID=@%s' % varlist, bfn])
    else:
        subprocess.run(['bcftools', 'view', '-Oz', '-o',  bmerge, '-S', bsamfn, '-i', 'ID=@%s' % varlist, bfn])
    subprocess.run(['tabix', '-p', 'vcf', bmerge])
    
    print('####################')
    print('Merging cohorts')
    subprocess.run(['bcftools', 'merge', '-m' ,'none', '-Oz', '-o', output, amerge, bmerge])
    
    print('####################')
    print('Merging done')
    print('Merged VCF saved to:\n\t%s' % output)
    print('Runtime: ' + str(datetime.now()-start) + '\n')
