import sys
import logging

def check_chrom(chrom):
    logging.info('####################')
    logging.info('Checking if valid chromosome...')
    min_dic = {1:49800,
               2:48440,
               3:39660,
               4:38050,
               5:36310,
               6:34170,
               7:31870,
               8:29030,
               9:27680,
               10:26760,
               11:27020,
               12:26660,
               13:22880,
               14:21410,
               15:20400,
               16:18070,
               17:16660,
               18:16080,
               19:11730,
               20:12890,
               21:9350,
               22:10170,
               999:0}
    try:
        chrom = int(chrom)
        min_num = min_dic[chrom]
        logging.info('Valid!')
    except:
        logging.error('Invalid chromosome. Accepts 1-22.')
        sys.exit(1)
    return min_num
