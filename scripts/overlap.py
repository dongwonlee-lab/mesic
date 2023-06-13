import sys
from datetime import datetime
from check_chrom import check_chrom


def load_varlist(var_fn):
    print('####################')
    print('Reading variants: %s' % var_fn)
    varlist = set()
    f1 = open(var_fn, 'r')
    for line in f1.readlines():
        var = line.strip()
        varlist.add(var)
    f1.close()
    print('%s variants detected' % len(varlist))
    return varlist


def run(chrom, afn, bfn, output):
    start = datetime.now()
    
    print('###########################')
    print('##### Variant Overlap #####')
    print('###########################')
    print('### Arguments')
    print('Chromosome\t%s' % chrom)
    print('Variant list A\t%s' % afn)
    print('Variant list B\t%s' % bfn)
    print('Output\t%s' % output)
    
    min_num = check_chrom(chrom)
    
    a_varlist = load_varlist(afn)
    b_varlist = load_varlist(bfn)
    overlap = a_varlist.intersection(b_varlist)
    f1 = open(output, 'w')
    f1.write('\n'.join(overlap))
    f1.close()

    print('####################')
    print(' %s variants overlapping' % len(overlap))

    if len(overlap) < min_num:
        print('WARNING: low number of variants passing QC for chromosome %s (ideally want %s variants or more)' % (chrom, min_num))


    print('####################')
    print('Finding variant overlap done.')
    print('Overlapping variants IDs saved to:\n\t%s' % output)
    print('Runtime: ' + str(datetime.now()-start) + '\n')
