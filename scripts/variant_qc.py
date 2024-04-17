from datetime import datetime
import sys
import logging
import gzip

from check_chrom import check_chrom

def load_rsq(rfn, rfilter, rcol, rvarcol, efilter, ecol, r2_drop=0, er2_drop=0, track_bin=500000, badcol=10): 
    wrong_col1 = 0
    wrong_col2 = 0
    
    logging.info('####################')
    logging.info('Loading Rsqs...')
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
                        if er2 != '-' and er2 != '.':
                            logging.debug('%s is not a number' % rsq)
                            wrong_col2 = wrong_col2+1
                        else:
                            rdic[var] = rsq
                else:
                    r2_drop = r2_drop+1
                count = count + 1
                if count%track_bin == 0:
                    logging.debug('Processed %s variants' % count)
            except:
                if rsq != '-' and rsq != '.':
                    logging.debug('%s is not a number' % rsq)
                    wrong_col1 = wrong_col1+1
            
            if wrong_col1 > badcol:
                logging.error("Check your column number for RSQ (-rc). It doesn't seem right.")
                logging.error('%s is not a number' % rsq)
                sys.exit(1)
                
            if wrong_col2 > badcol:
                logging.error("Check your column number for ER2 (-ec). It doesn't seem right.")
                logging.error('%s is not a number' % er2)
                sys.exit(1)
        logging.info('Processed %s variants' % count)
    return (rdic, r2_drop, er2_drop)


def load_maf(mfn, mfilter, mcol, mvarcol, maf_drop=0, track_bin=500000, badcol=10):
    wrong_col = 0

    logging.info('####################')
    logging.info('Loading AFs...')
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
                    logging.debug('Processed %s variants' % count)
            except:
                if maf != '-' and maf != '.':
                    logging.debug('%s is not a number' % maf)
                    wrong_col = wrong_col+1
                    
            if wrong_col > badcol:
                logging.error("Check your column number for MAF (-mc). It doesn't seem right.")
                logging.error('%s is not a number' % maf)
                sys.exit(1)
                
        logging.info('Processed %s variants' % count)
    return (mdic, maf_drop)


def load_hwe(hfile, hfilter, nocases, track_bin=500000, badcol=10):
    wrong_col = 0
    
    logging.info('####################')
    logging.info('Loading HWE p-values...')
    
    if nocases:
        logging.debug('No cases. Using ALL for HWE filter.')
    
    count = 0
    hdic = {} # only if fail HWE
    
    with open(hfile) as f1:
        for line in f1:
            ln = line.strip().split()
            var = ln[1]
            hwe = ln[8]
            test = [ln[2]][0]
                
            try:
                hwe = float(hwe)
                
                if nocases:
                    if hwe <= hfilter and 'ALL' in test:
                        hdic[var] = hwe
                else:
                    if hwe <= hfilter and test == 'UNAFF':
                        hdic[var] = hwe
                        
                count = count+1
                if count%track_bin == 0:
                    logging.debug('Processed %s lines' % count)
            except:
                logging.debug('%s is not a number' % hwe)
                wrong_col = wrong_col+1
                    
            if wrong_col > badcol:
                logging.error("Check your HWE file. It doesn't seem right.")
                logging.error('%s is not a number' % hwe)
                sys.exit(1)
                
        logging.info('Processed %s lines (may be multiple of total variants)' % count)
    return hdic


def run(chrom, rfn, mfn, hfn, ofn,
        rfilter, mfilter, efilter, hfilter,
        rcol, mcol, ecol, rvarcol, mvarcol, 
        nocases):
    start = datetime.now()
    
    logging.info('######################')
    logging.info('##### Variant QC #####')
    logging.info('######################')
    logging.info('### Arguments')
    logging.info('Chromosome\t%s' % chrom)
    logging.info('Rsq file\t%s' % rfn)
    logging.info('Rsq filter\t%s' % rfilter)
    logging.info('ER2 filter\t%s' % efilter)
    logging.info('Rsq column\t%s' % rcol)
    logging.info('ER2 column\t%s' % ecol)
    logging.info('Rsq variant column\t%s' % rvarcol)
    logging.info('MAF file\t%s' % mfn)
    logging.info('MAF filter\t%s' % mfilter)
    logging.info('MAF column\t%s' % mcol)
    logging.info('MAF variant column\t%s' % mvarcol)
    logging.info('HWE file\t%s' % hfn)
    if hfn != None:
        logging.info('HWE filter\t%s' % hfilter)
    logging.info('Output\t%s' % ofn)
    
    min_num = check_chrom(chrom)   
    
    (rdic, r2_drop, er2_drop) = load_rsq(rfn, rfilter, rcol-1, rvarcol-1, efilter, ecol-1)
    (mdic, maf_drop) = load_maf(mfn, mfilter, mcol-1, mvarcol-1)
    
    if hfn != None:
        hdic = load_hwe(hfn, hfilter, nocases)
    
    logging.info('####################')
    logging.info('Finding high quality variants')
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
                logging.debug('Processed %s/%s variants' % (count, total))
        logging.info('Processed %s/%s variants' % (count, total))
        f1.close()
    else:
        f1 = open(ofn, 'w')
        for var in rdic:
            if var in mdic.keys() and var not in hdic.keys():
                f1.write('%s\n' % var)
                hq_keep = hq_keep+1
            count = count+1
            if count%10000 == 0:
                logging.debug('Processed %s/%s variants' % (count, total))
        logging.info('Processed %s/%s variants' % (count, total))
        f1.close()

    logging.info('####################')
    logging.info(' %s variants fail RSQ filter (%s)' % (r2_drop, rfilter))
    logging.info(' %s genotyped variants fail ER2 filter (%s)' % (er2_drop, efilter))
    logging.info(' %s variants fail MAF filter (%s)' % (maf_drop, mfilter))
    if hfn != None:
        logging.info(' %s variants fail HWE (%s)' % (len(hdic.keys()), hfilter))
    logging.info(' %s variants pass all filters' % hq_keep)

    if hq_keep < min_num:
        logging.warning('WARNING: low number of variants passing QC for chromosome %s (ideally want %s variants or more)' % (chrom, min_num))

    logging.info('####################') 
    logging.info('Filtering high quality variants done.')
    logging.info('High quality variants IDs saved to:\n\t%s' % ofn)
    logging.info('Runtime: ' + str(datetime.now()-start) + '\n')
