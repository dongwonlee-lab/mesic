from datetime import datetime
import cyvcf2


def calculate_rsq(d, p_hat, rsq_topmed, term1):    
    term2 = sum(sum((d-p_hat)**2))    
    term3 = p_hat*(1-p_hat)
    
    rsq = (term1*term2)/term3
    
    return (p_hat, rsq, rsq_topmed)


def run(vcf_fn, out_fn, sam_fn):
    start = datetime.now()

    print('###########################')
    print('##### Calculating Rsq #####')
    print('###########################')
    print('### Arguments')
    print('VCF\t%s' % vcf_fn)
    print('Output\t%s' % out_fn)
    print('Samples\t%s' % sam_fn)
    
    print('####################')
    print('Getting list of samples for calculation...')
    if sam_fn != None:
        fo = open(sam_fn, 'r')
        samples = fo.read().split('\n')[:-1]
        fo.close()        
        vcf = cyvcf2.VCF(fname=vcf_fn, samples=samples)
    else:
        print('Using all samples in VCF.')
        vcf = cyvcf2.VCF(fname=vcf_fn)
        samples = vcf.samples
        
    n_samples = len(samples)
    print('Number of samples: %s' % n_samples)
    print("Samples' head: %s" % samples[:5])
    
    
    print('####################')
    print('Calculating rsq...')
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
        af = sum(sum(d))/(2*n_samples) 
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
            
            print('%s variants processed...' % (step*count))
    
    vcf.close()
    fo.close()
    print('%s variants processed...' % n_variants)

    print('Number of variants non-monomorphic: %s' % passed)
    print('Number of variants monomorphic (rsq could not be calculated): %s' % failed)
    
    print('####################')
    print('Rsq calculations done.')
    print('Rsq calculations saved to:\n\t%s' % out_fn)
    print('Runtime: ' + str(datetime.now()-start) + '\n')
                         