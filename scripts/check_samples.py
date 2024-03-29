import os
import sys
import logging


def check_samples(filelists, python_lib):
    if (python_lib != ''):
        sys.path.append(python_lib)
        
    import cyvcf2
    
    logging.info('####################')
    logging.info('Checking if all samples are in VCFs...')
    
    samples_txt = set()
    samples_vcf = set()
    
    for idx in filelists:
        vcf_fn = filelists[idx][0]
        sample_fn = filelists[idx][2]
        
        vcf = cyvcf2.VCF(fname=vcf_fn).samples
        for y in vcf:
            samples_vcf.add(y)
            
        f1 = open(sample_fn, 'r')
        txt = [line.strip() for line in f1.readlines()]
        f1.close()
        for x in txt:
            samples_txt.add(x)
            
    if samples_txt.issubset(samples_vcf):
        logging.info('All samples in sample lists are in at least 1 VCF.')
        logging.info(" _________________________________________")
        logging.info("/   Please ignore bcftools warn messages  \\")
        logging.info("\       about samples not existing.       /")
        logging.info(" -----------------------------------------")
        logging.info("                           \\")
        logging.info("                            \\")       
        logging.info("                             .---.         ,,")
        logging.info("                  ,,        /     \       ;,,'")
        logging.info("                 ;, ;      (  o  o )      ; ;")
        logging.info("                   ;,';,,,  \  \/ /      ,; ;")
        logging.info("                ,,,  ;,,,,;;,`   '-,;'''',,,'")
        logging.info("               ;,, ;,, ,,,,   ,;  ,,,'';;,,;''';")
        logging.info("                  ;,,,;    ~~'  '';,,''',,;''''")
        logging.info("                                     ''''")
        
    else:
        logging.error(samples_txt - samples_vcf)
        logging.error('These samples are not in any VCFs.')
        sys.exit(1)
    
        