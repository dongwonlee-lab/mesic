import sys

def check_chrom(chrom):
    print('####################')
    print('Checking if valid chromosome...')
    min_dic = {1:99600,
               2:96900,
               3:79400,
               4:76100,
               5:72700,
               6:68400,
               7:63800,
               8:58100,
               9:55400,
               10:53600,
               11:54100,
               12:53400,
               13:45800,
               14:42900,
               15:40800,
               16:36200,
               17:33400,
               18:32200,
               19:23500,
               20:25800,
               21:18700,
               22:20400}
    try:
        chrom = int(chrom)
        min_num = min_dic[chrom]
        print('Valid!')
    except:
        print('Invalid chromosome. Accepts 1-22.')
        sys.exit()
    return min_num
