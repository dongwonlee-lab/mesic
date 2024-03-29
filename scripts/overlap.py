from datetime import datetime
import sys
import logging

from check_chrom import check_chrom

def load_varlist(var_fn):
    logging.info('####################')
    logging.info('Reading variants: %s' % var_fn)
    varlist = set()
    f1 = open(var_fn, 'r')
    for line in f1.readlines():
        var = line.strip()
        varlist.add(var)
    f1.close()
    logging.info('%s variants detected' % len(varlist))
    return varlist


def run(filelist_fn, chrom, output):
    start = datetime.now()
    
    filelist = set()
    f1 = open(filelist_fn, 'r')
    for line in f1.readlines():
        fn = line.strip()
        filelist.add(fn)
    f1.close()
    
    logging.info('###########################')
    logging.info('##### Variant Overlap #####')
    logging.info('###########################')
    logging.info('### Arguments')
    logging.info('Chromosome:\t%s' % chrom)
    logging.info('Variant lists:\n\t%s' % '\n\t'.join(filelist))
    logging.info('Output:\t%s' % output)
    
    min_num = check_chrom(chrom)
    
    varlist = {}
    for fn in filelist:
        varlist[fn] = load_varlist(fn)
        
    skip = True
    for fn in varlist.keys():
        if skip:
            list1 = varlist[fn]
            skip = False
        else:
            list2 = varlist[fn]
            overlap = list1.intersection(list2)
            list1 = overlap
            
    overlap = sorted(overlap)
    f1 = open(output, 'w')
    f1.write('\n'.join(overlap))
    f1.write('\n')
    f1.close()

    logging.info('####################')
    logging.info(' %s variants overlapping' % len(overlap))

    if len(overlap) < min_num:
        logging.warning('WARNING: low number of variants passing QC for chromosome %s (ideally want %s variants or more)' % (chrom, min_num))


    logging.info('####################')
    logging.info('Finding variant overlap done.')
    logging.info('Overlapping variant IDs saved to:\n\t%s' % output)
    logging.info('Runtime: ' + str(datetime.now()-start) + '\n')
