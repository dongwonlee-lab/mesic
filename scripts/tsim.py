import os
import sys
import logging
import argparse


global_parser = argparse.ArgumentParser(prog='tsim.py',
                                  description='##########################################\n### Two-Stage Imputation Method ###\n##########################################\n\nConducts per chromosome analysis.',
                                  formatter_class=argparse.RawTextHelpFormatter,
                                  prefix_chars='-',
                                  argument_default=None,
                                  add_help=True,
                                  allow_abbrev=True,
                                  exit_on_error=True)
subparsers = global_parser.add_subparsers(title='subcommands',
                                          dest='command',
                                          description=None, 
                                          required=True)

global_parser.add_argument('--version', action='version', version='%(prog)s 1.0 (2023)')


rsq_parser = subparsers.add_parser('rsq', 
                                   formatter_class=argparse.RawTextHelpFormatter,
                                   description='Calculates Minimac4 rsq from imputed dosages.\nOutputs TSV containing variant ID, AAF, recalculated RSQ, original Rsq, and empirical Rsq (if genotyped).\nParticularly useful if you only want to analyze a subset of the samples that were imputed.',
                                   help='calculate Rsq')
rsq_parser.add_argument('-v', '--vcf',
                        required=True,
                        help='VCF containing dosage info')
rsq_parser.add_argument('-o', '--output',
                        required=True,
                        help='output file (TSV)')
rsq_parser.add_argument('-s', '--samples',
                        help='file containing list of samples to include in calculation')
rsq_parser.add_argument('-p', '--pythonlib', default='',
                        help='specify python site-packages location')
rsq_parser.add_argument('--verbose', action='store_true',
                        help='run with more verbose logging')
                        

qc_parser = subparsers.add_parser('qc', 
                                  formatter_class=argparse.RawTextHelpFormatter,
                                  description='Perform variant QC based on Rsq, empirical Rsq (if genotyped), allele frequency, and HWE.\nOutputs file containing variants passing filters, one per line.',
                                  help='perform variant QC')
qc_parser.add_argument('-r', '--rsq', required=True,
                       help='TSV containing Rsq and empirical rsq calculations (can be gzipped)')
qc_parser.add_argument('-m', '--maf', required=True,
                       help='TSV containing allele frequency (can be gzipped)')
qc_parser.add_argument('-o', '--output', required=True,
                       help='output file (txt)')
qc_parser.add_argument('-c', '--chrom', required=True,
                       help='chromosome of analysis')
qc_parser.add_argument('--hwe', default=None,
                       help='file containing PLINK `--hardy` output')
qc_parser.add_argument('-rf', '--rfilter', default=0.99, type=float,
                       help='Rsq filter [default: 0.99]')
qc_parser.add_argument('-mf', '--mfilter', default=0.01, type=float,
                       help='MAF filter [default: 0.01]')
qc_parser.add_argument('-ef', '--efilter', default=0.9, type=float,
                       help='empirical Rsq filter [default: 0.90]')
qc_parser.add_argument('-hf', '--hfilter', default=1e-6, type=float,
                       help='HWE p-value filter [default: 1e-6]')
qc_parser.add_argument('-rc', '--rcol', default=3, type=int,
                       help='column # containing Rsq (1-based) [default: 3]')
qc_parser.add_argument('-mc', '--mcol', default=2, type=int,
                       help='column # containing allele frequency (1-based) [default: 2]')
qc_parser.add_argument('-ec', '--ecol', default=5, type=int,
                       help='column # containing empirical Rsq (1-based) [default: 5]')
qc_parser.add_argument('-rvc', '--rvarcol', default=1, type=int,
                       help='column # containing variant ID in rsq file (1-based) [default: 1]')
qc_parser.add_argument('-mvc', '--mvarcol', default=1, type=int,
                       help='column # containing variant ID in maf file (1-based) [default: 1]')
qc_parser.add_argument('-nc', '--nocases', action='store_true',
                       help='if used, indicates there are no cases in QC (relevant for HWE filtering)')
qc_parser.add_argument('--verbose', action='store_true',
                       help='run with more verbose logging')


overlap_parser = subparsers.add_parser('overlap', 
                                       formatter_class=argparse.RawTextHelpFormatter,
                                       description='Find variant overlap between multiple cohorts.\nOutputs file containing variants present in both cohorts, one per line.',
                                       help='find variant overlap')
overlap_parser.add_argument('-l', '--varlist', required=True,
                            help='file containing list of variant lists to overlap')
overlap_parser.add_argument('-o', '--output', required=True,
                            help='output file (txt)')
overlap_parser.add_argument('-c', '--chrom', required=True,
                            help='chromosome of analysis')


merge_parser = subparsers.add_parser('merge', 
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='Merge VCF files from multiple cohorts based on variant overlap.\nOutputs merged VCF file.',
                                     help='merge VCF files')
merge_parser.add_argument('-l', '--list', required=True,
                          help='common-separated file contianing list of VCFs to merge, variant lists to merge on, and samples to include for each file. Column 1 = VCF files, Column 2 = variant lists, Column 3 = sample lists.')
merge_parser.add_argument('-o', '--output', required=True,
                          help='output VCF')
merge_parser.add_argument('-c', '--chrom', required=True,
                          help='chromosome of analysis')
merge_parser.add_argument('-t', '--tempdir', default=os.getcwd(),
                          help='directory for temporary files [default is current working directory]')
merge_parser.add_argument('-s', '--snpsonly', action='store_true',
                          help='filter for SNPs only')
merge_parser.add_argument('-p', '--pythonlib', default='',
                        help='specify python site-packages location')


args = global_parser.parse_args()


logfmt_str = '%(levelname)s %(asctime)s: %(message)s'
datefmt_str = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(stream=sys.stdout,
            format=logfmt_str, datefmt=datefmt_str,
            level=logging.INFO)
logger = logging.getLogger()

if args.command == 'rsq':
    import calculate_rsq
    
    if args.verbose:
        logging.info('Verbosity on')
        logger.setLevel('DEBUG')
    calculate_rsq.run(args.vcf, args.output, args.samples, args.pythonlib)
    
if args.command == 'qc':
    import variant_qc
    
    if args.verbose:
        logging.info('Verbosity on')
        logger.setLevel('DEBUG')
    variant_qc.run(args.chrom, args.rsq, args.maf, args.hwe, args.output,
                   args.rfilter, args.mfilter, args.efilter, args.hfilter,
                   args.rcol, args.mcol, args.ecol, args.rvarcol, args.mvarcol, 
                   args.nocases)

if args.command == 'overlap':
    import overlap
    
    # no debug logging
    
    overlap.run(args.varlist, args.chrom, args.output)
    
if args.command == 'merge':
    import merge
    
    # no debug logging
    
    merge.run(args.chrom, args.list, args.output, args.tempdir, args.snpsonly, args.pythonlib)

