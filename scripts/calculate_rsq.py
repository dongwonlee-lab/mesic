from datetime import datetime
import sys
import logging

def calculate_rsq(d, p_hat, rsq_topmed, term1):    
    term2 = ((d-p_hat)**2).sum()
    term3 = p_hat*(1-p_hat)
    
    rsq = (term1*term2)/term3
    
    return (p_hat, rsq, rsq_topmed)


def run(vcf_fn, out_fn, sam_fn, python_lib):
    if (python_lib != ''):
        sys.path.append(python_lib)
   
    import cyvcf2
    
    start = datetime.now()

    logging.info('###########################')
    logging.info('##### Calculating Rsq #####')
    logging.info('###########################')
    logging.info('### Arguments')
    logging.info('VCF\t%s' % vcf_fn)
    logging.info('Output\t%s' % out_fn)
    logging.info('Samples\t%s' % sam_fn)
    
    logging.info('####################')
    logging.info('Getting list of samples for calculation...')
    if sam_fn != None:
        fo = open(sam_fn, 'r')
        samples = fo.read().split('\n')[:-1]
        fo.close()        
        vcf = cyvcf2.VCF(fname=vcf_fn, samples=samples)
    else:
        logging.info('Using all samples in VCF.')
        vcf = cyvcf2.VCF(fname=vcf_fn)
        samples = vcf.samples
        
    n_samples = len(samples)
    logging.info('Number of samples: %s' % n_samples)
    logging.debug("Samples' head: %s" % samples[:5])
    
    
    logging.info('####################')
    logging.info('Calculating rsq...')
    t1 = 1/(2*n_samples) 
    fo = open(out_fn, 'w')
    fo.write('ID\tAAF\tRSQ\tRSQ_TOPMED\tER2\n')
    
    passed = 0
    failed = 0
    n_variants = 0
    count = 0
    step = 500000
    
    for variant in vcf:
        n_variants = n_variants + 1
        var_id = variant.ID  
        
        # estimated haploid alternate allele dosage
        d = variant.format('HDS')
        assert d.shape == (n_samples, 2)

        # calculate alternative allele frequency
        af = d.sum()/(2*n_samples) 
        maf = af > 0 and af < 1 # monomorphic sites not included

        rsq = variant.INFO.get('R2')
        er2 = variant.INFO.get('ER2')
        if er2 == None:
            er2 = '-'

        # calculate rsq
        if maf:
            passed = passed + 1
            rsq_calc = calculate_rsq(d, af, rsq, t1)
            fo.write('%s\t' % var_id)
            fo.write('\t'.join(list(map(str, rsq_calc))) + '\t%s\n' % er2)
        else:
            # (p_hat, rsq, rsq_topmed, er2)
            rsq_calc = (0, '-', rsq)
            failed = failed+1
            fo.write('%s\t' % var_id)
            fo.write('\t'.join(list(map(str, rsq_calc))) + '\t%s\n' % er2)

        if n_variants % step == 0:
            count = count + 1
            
            logging.debug('%s variants processed...' % (step*count))
    
    vcf.close()
    fo.close()
    logging.info('%s variants processed...' % n_variants)

    logging.info('Number of variants non-monomorphic: %s' % passed)
    logging.info('Number of variants monomorphic (rsq could not be calculated): %s' % failed)
    
    logging.info('####################')
    logging.info('Rsq calculations done.')
    logging.info('Rsq calculations saved to:\n\t%s' % out_fn)
    logging.info('Runtime: ' + str(datetime.now()-start) + '\n')
                         