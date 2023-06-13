from datetime import datetime
import gzip
from check_chrom import check_chrom


def load_rsq(rfn, rfilter, rcol, rvarcol, efilter, ecol, r2_drop=0, er2_drop=0, track_bin=500000):  
    print('####################')
    print('Loading Rsqs...')
    count = 0
    rdic = {} # only if pass Rsq filter and ER2 filter
    rsq_zip = True if rfn.endswith('.gz') else False
    with gzip.open(rfn, 'rt') if rsq_zip else open(rfn) as f1:
        for line in f1:
            ln = line.strip().split('\t')
            var = ln[rvarcol]
            rsq = ln[rcol]
            er2 = ln[ecol]
            try:
                rsq = float(rsq)
                if rsq >= rfilter:
                    try:
                        er2 = float(er2)
                        if er2 >= efilter:
                            rdic[var] = rsq
                        else:
                            er2_drop = er2_drop+1
                    except:
                        rdic[var] = rsq
                else:
                    r2_drop = r2_drop+1
                count = count + 1
                if count%track_bin == 0:
                    print('Processed %s variants' % count)
            except:
                if rsq != '-':
                    print('%s is not a number' % rsq)
        print('Processed %s variants' % count)
    return (rdic, r2_drop, er2_drop)


def load_maf(mfn, mfilter, mcol, mvarcol, maf_drop=0, track_bin=500000):
    print('####################')
    print('Loading AFs...')
    count = 0
    mdic = {} # only if pass MAF
    maf_zip = True if mfn.endswith('.gz') else False
    with gzip.open(mfn, 'rt') if maf_zip else open(mfn) as f1:
        for line in f1:
            ln = line.strip().split('\t')
            var = ln[mvarcol]
            maf = ln[mcol]
            try:
                maf = float(maf)
                if maf >= mfilter and maf <= 1-mfilter:
                    mdic[var] = maf
                else:
                    maf_drop = maf_drop+1
                count = count+1
                if count%track_bin == 0:
                    print('Processed %s variants' % count)
            except:
                if maf != '-':
                    print('%s is not a number' % maf)
        print('Processed %s variants' % count)
    return (mdic, maf_drop)


def load_hwe(hfile, hfilter, track_bin=500000):
    print('####################')
    print('Loading HWE p-values...')
    count = 0
    hdic = {} # only if fail HWE
    with open(hfile) as f1:
        for line in f1:
            ln = line.strip().split()
            var = ln[1]
            hwe = ln[8]
            test = [ln[2] if len(ln) == 9 else 'UNAFF'][0]
            try:
                hwe = float(hwe)
                if hwe <= hfilter and test == 'UNAFF':
                    hdic[var] = hwe
                count = count+1
                if count%track_bin == 0:
                    print('Processed %s lines' % count)
            except:
                print('%s is not a number' % hwe)
        print('Processed %s lines (may be multiple of total variants)' % count)
    return hdic


def run(chrom, rfn, mfn, hfn, ofn,
        rfilter, mfilter, efilter, hfilter,
        rcol, mcol, ecol, rvarcol, mvarcol):
    start = datetime.now()
    
    print('######################')
    print('##### Variant QC #####')
    print('######################')
    print('### Arguments')
    print('Chromosome\t%s' % chrom)
    print('Rsq file\t%s' % rfn)
    print('Rsq filter\t%s' % rfilter)
    print('ER2 filter\t%s' % efilter)
    print('Rsq column\t%s' % rcol)
    print('ER2 column\t%s' % ecol)
    print('Rsq variant column\t%s' % rvarcol)
    print('MAF file\t%s' % mfn)
    print('MAF filter\t%s' % mfilter)
    print('MAF column\t%s' % mcol)
    print('MAF variant column\t%s' % mvarcol)
    print('HWE file\t%s' % hfn)
    if hfn != None:
        print('HWE filter\t%s' % hfilter)
    print('Output\t%s' % ofn)
    
    min_num = check_chrom(chrom)   
    
    (rdic, r2_drop, er2_drop) = load_rsq(rfn, rfilter, rcol-1, rvarcol-1, efilter, ecol-1)
    (mdic, maf_drop) = load_maf(mfn, mfilter, mcol-1, mvarcol-1)
    
    if hfn != None:
        hdic = load_hwe(hfn, hfilter)
    
    print('####################')
    print('Finding high quality variants')
    count = 0
    hq_keep = 0
    total = len(rdic)

    if hfn == None:
        f1 = open(ofn, 'w')
        for var in rdic:
            if var in mdic.keys():
                f1.write('%s\n' % var)
                hq_keep = hq_keep+1
            count = count+1
            if count%10000 == 0:
                print('Processed %s/%s variants' % (count, total))
        print('Processed %s/%s variants' % (count, total))
        f1.close()
    else:
        f1 = open(ofn, 'w')
        for var in rdic:
            if var in mdic.keys() and var not in hdic.keys():
                f1.write('%s\n' % var)
                hq_keep = hq_keep+1
            count = count+1
            if count%10000 == 0:
                print('Processed %s/%s variants' % (count, total))
        print('Processed %s/%s variants' % (count, total))
        f1.close()

    print('####################')
    print(' %s variants fail RSQ filter (%s)' % (r2_drop, rfilter))
    print(' %s genotyped variants fail ER2 filter (%s)' % (er2_drop, efilter))
    print(' %s variants fail MAF filter (%s)' % (maf_drop, mfilter))
    if hfn != None:
        print(' %s variants fail HWE (%s)' % (len(hdic.keys()), hfilter))
    print(' %s variants pass all filters' % hq_keep)

    if hq_keep < min_num:
        print('WARNING: low number of variants passing QC for chromosome %s (ideally want %s variants or more)' % (chrom, min_num))

    print('####################') 
    print('Filtering high quality variants done.')
    print('High quality variants IDs saved to:\n\t%s' % ofn)
    print('Runtime: ' + str(datetime.now()-start) + '\n')
