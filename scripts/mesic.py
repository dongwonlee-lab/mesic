import os
import argparse


global_parser = argparse.ArgumentParser(prog='mesic.py',
                                  description='##########################################\n### Merging Separately Imputed Cohorts ###\n##########################################\n\nConducts per chromosome analysis.',
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
                        

qc_parser = subparsers.add_parser('qc', 
                                  formatter_class=argparse.RawTextHelpFormatter,
                                  description='Perform variant QC based on Rsq, empirical Rsq (if genotyped), allele frequency, and HWE.\nOutputs file containing variants passing filters, one per line.',
                                  help='perform variant QC')
qc_parser.add_argument('-r', '--rsq', required=True,
                       help='TSV containing Rsq and empirical rsq calculations (can be gzipped)')
qc_parser.add_argument('-m', '--maf', required=True,
                       help='TSV containing allele frequency (can be gzipped)')
qc_parser.add_argument('-o', '--output', required=True,
                       help='output file')
qc_parser.add_argument('-c', '--chrom', required=True,
                       help='chromosome of analysis')
qc_parser.add_argument('--hwe', default=None,
                       help='file containing PLINK `--hardy` output')
qc_parser.add_argument('-rf', '--rfilter', default=0.30, type=float,
                       help='Rsq filter [default: 0.30]')
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


overlap_parser = subparsers.add_parser('overlap', 
                                       formatter_class=argparse.RawTextHelpFormatter,
                                       description='Find variant overlap between 2 cohorts.\nOutputs file containing variants present in both cohorts, one per line.',
                                       help='find variant overlap')
overlap_parser.add_argument('-a', '--avarlist', required=True,
                            help='list of variants in cohort A passing QC')
overlap_parser.add_argument('-b', '--bvarlist', required=True,
                            help='list of variants in cohort B passing QC')
overlap_parser.add_argument('-o', '--output', required=True,
                            help='output file')
overlap_parser.add_argument('-c', '--chrom', required=True,
                            help='chromosome of analysis')


merge_parser = subparsers.add_parser('merge', 
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='Merge VCF files from 2 cohorts based on variant overlap.\nOutputs merged VCF file.',
                                     help='merge VCF files')
merge_parser.add_argument('-a', '--avcf', required=True,
                          help='VCF for cohrot A')
merge_parser.add_argument('-b', '--bvcf', required=True,
                          help='VCF for cohort B')
merge_parser.add_argument('-v', '--varlist', required=True,
                          help='list of variants to keep')
merge_parser.add_argument('-o', '--output', required=True,
                          help='output VCF')
merge_parser.add_argument('-c', '--chrom', required=True,
                          help='chromosome of analysis')
merge_parser.add_argument('-as', '--asamples', default=None,
                          help='list of cohort A samples to keep')
merge_parser.add_argument('-bs', '--bsamples', default=None,
                          help='list of cohort B samples to keep')
merge_parser.add_argument('-t', '--tempdir', default=os.getcwd(),
                          help='directory for temporary files [default is current working directory]')


args = global_parser.parse_args()


if args.command == 'rsq':
    import calculate_rsq
    calculate_rsq.run(args.vcf, args.output, args.samples)

if args.command == 'qc':
    import variant_qc
    variant_qc.run(args.chrom, args.rsq, args.maf, args.hwe, args.output,
                   args.rfilter, args.mfilter, args.efilter, args.hfilter,
                   args.rcol, args.mcol, args.ecol, args.rvarcol, args.mvarcol)

if args.command == 'overlap':
    import overlap
    overlap.run(args.chrom, args.avarlist, args.bvarlist, args.output)
    
if args.command == 'merge':
    import merge
    merge.run(args.chrom, args.avcf, args.bvcf, args.varlist, args.output, 
              args.asamples, args.bsamples, args.tempdir)

