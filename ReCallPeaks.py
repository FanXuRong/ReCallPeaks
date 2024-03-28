import pandas as  pd
import glob2 as glob
import re
import numpy as np
import multiprocessing as mp
import os
import scipy
import math
from optparse import OptionParser

class ReCallPeaks:
    ## 将 deeptools 生成的 coverage 文件按照需要的长度切片
    ## 接受两个参数 (切片长度)，bedgraph文件，为多线程设计
    @staticmethod
    def Fliter_Bedgraph (binsize_length: int,File)->None:
        outfile = File + 'fliter.csv'
        a = os.getpid()
        print (f'PID is {a}\noutfile is {outfile}\n')
        a = pd.read_csv(File,sep='\t',header=None,names=['Chr','start','end','score'],dtype={'Chr':'str','start':'int','end':'int','score':'int'})
        a.loc[:,'index_number'] = a['end'] - a['start']
        insert = a[a['index_number']>binsize_length]
        use_df = a[a['index_number']==binsize_length]
        select_df = ReCallPeaks.Insert(binsize_length,insert)
        result = pd.concat([select_df,use_df])
        result = result.sort_values(by=['Chr', 'start']).reset_index(drop=True)
        print (result.head(3))
        result_df = result.drop('index_number',axis=1)
        result_df.to_csv(outfile,index=False)
    '''
    对于每个 样本 的bedgraph文件,均返回一个 {file}fliter.csv 文件

    Chr,start,end,score
    '''
    ## 实际的切片函数,直接被 Fliter_Bedgraph 使用
    ## 接受被 Fliter_Bedgraph 筛选的大于滑窗的片段
    @staticmethod
    def Insert (length,data):
        insert = []
        for idx, row in data.iterrows():
            for i in range(row['start'], row['end'], length):
                insert.append({'Chr': row['Chr'], 'start': i, 'end': min(i+length, row['end']), 'score': row['score']})
        select_df = pd.DataFrame(insert)
        return select_df
    '''
    返回一个 dataframe
    '''

    ## 返回执行多线程中返回的 error and warning
    @staticmethod
    def call_back(res):
        print(f'successful! warning is {res}')
    @staticmethod
    def err_call_back(err):
        print(f'error!：{str(err)}')
    ############################################################################################################################
    ## 读入处理后的 *fliter.bedgraph文件, 只取Coverage那一行
    @staticmethod
    def Read_Fliter_File(File):
        data = pd.read_csv(File)
        Series = data.loc[:,'score']
        return Series
    '''
    返回每个重复的 Coverage , 数据类型为Series
    '''
    ## 处理空白组（除数），因为是除法运算，空白组<1 ，赋值为 1，不小于 1,返回实际值
    @staticmethod
    def Blank_compute (Dataframe):
        Series = (Dataframe['TB0135r1']+Dataframe['TB0135r2'])/2
        if (Series < 1):
            return 1
        else:
            return Series

    ## 计算处理后的 bedgraph 总矩阵的处理组, 获取每个片段 coverage 对于对照组倍数
    ## 接受样本字典 {gene:[rep_1,rep_2]},与 Coverage 总矩阵 dataframe
    @staticmethod
    def Sample_compute (Samples: dict,DataFrame):
        for i in Samples:
            a = len(Samples[i])
            Sample_threshold = i + '_threshold'     ## 此列代表 处理组与对照组的倍数
            Sample_Mean = i + '_Mean'               ## 此列代表平均数
            if a == 1:
                DataFrame.loc[:,Sample_Mean] = DataFrame[Samples[i][0]]
                DataFrame.loc[:,Sample_threshold] = (DataFrame[Samples[i][0]]/DataFrame['TB0135'])
                #DataFrame.loc[:,i] = (DataFrame[Sample_Mean]/DataFrame[Sample_threshold])     ## 此列代表,样品的差异倍数与Coverage的绝对值,用于抑制假阳性
            else:
                DataFrame.loc[:,Sample_Mean] = (DataFrame[Samples[i][0]] + DataFrame[Samples[i][1]]) / 2
                DataFrame.loc[:,Sample_threshold] = (DataFrame[Sample_Mean]/DataFrame['TB0135'])
                #DataFrame.loc[:,i] = (DataFrame[Sample_Mean]/DataFrame[Sample_threshold])

    '''
    返回每个基因重复的均值以及倍数
    '''

    ## 输出每个 gene 每个重复的 Coverage, 与 Coverage 的平均值
    ## 接受三个参数 包含处理组倍数、处理组均值的结果均值, 样本字典 {gene:[rep_1,rep_2]},以及重复的绝对值阈值
    @staticmethod
    def SelectData (DataFrame,GeneList: list,Mean_Coverage: int):
        for sample in GeneList:
            Sample_Mean = sample + '_Mean'
            Sample_threshold = sample + '_threshold'
            Selected_df = DataFrame.loc[(DataFrame[Sample_threshold] >= 2) & (DataFrame[Sample_Mean] > Mean_Coverage)]
            Reselt_df = Selected_df[['Chr','start','end',Sample_Mean,Sample_threshold,'TB0135']]
            outfile = sample + '_Result.csv'
            Sorted_Reselt_df = Reselt_df.sort_values(by=Sample_Mean,ascending=False)
            Sorted_Reselt_df.to_csv(outfile,index=False)
    '''
    返回 '*_Result.csv文件'
    'Chr','start','end',sample(compare blank),Sample_Mean,'TB0135'
    '''

    @staticmethod
    ## 寻找最高峰，限定Coverage均值，写出文件 '_Peak'
    ## 接受6个参数 总结果矩阵,基因列表,峰高阈值,峰间距,threshold倍数阈值,滑窗长度
    def FindPeaks (Dataframe,sample,Peak_height: list,Peak_space: int,Mean_threshold: int,Windows_size) -> None:
        FindPeaks = []
        Sample_Mean = sample + '_Mean'
        Sample_threshold = sample + '_threshold'
        outfile = sample + '_Peaks.message.csv'
        bedfile = sample + '_Peaks.select.region.bed'
        Top_out_1000_file = sample + 'top1000.select.region'
        a = os.getpid()
        print (f'The PID is {a}, outfile is {outfile} {Top_out_1000_file} {bedfile} \n')
        ## 将对应 gene 的倍数列转成数组
        Testnp = Dataframe[Sample_Mean].to_numpy()
        ## 寻峰函数,接受峰高[min,max],以及峰间距 int
        Peaks_index,_ = scipy.signal.find_peaks(Testnp,height=Peak_height,distance=Peak_space)
        ## 将寻找到的每一个峰添加到新的表格中
        for index in Peaks_index:
            result = Dataframe.iloc[index]
            FindPeaks.append(result)
        FindPeaksResult = pd.DataFrame(FindPeaks)
        ## 按照 Coverage 限定阈值
        result = FindPeaksResult.loc[FindPeaksResult[Sample_threshold] > Mean_threshold]
        ## 提取 threshold 对应列
        FindPeaksOut = result[['Chr','start','end',Sample_threshold,Sample_Mean,'TB0135']]
        ## 排序
        FindPeaksResultOut = FindPeaksOut.sort_values(by=[Sample_Mean],ascending=False).reset_index(drop=True)
        ## 获取峰高前1000个
        Top_One_th = FindPeaksResultOut.head(1000)[['Chr','start','end']]
        Top_One_th.loc[:,'bedstart'] = Top_One_th.apply(ReCallPeaks.Addstart,args=(Windows_size,),axis=1)
        Top_One_th.loc[:,'bedend'] = Top_One_th['bedstart'] + 5*Windows_size
        Top_Result = Top_One_th.sort_values(by=['Chr','start'])[['Chr','bedstart','bedend']]
        Top_Result.to_csv(Top_out_1000_file,sep='\t',index=False,header=None)
        ## 输出每个峰的信息 csv
        FindPeaksResultOut.to_csv(outfile,index=False)
        BedOut = FindPeaksOut[['Chr','start','end']]
        ## 由于寻峰函数本身的限制,寻找到的不一定在峰顶,这里取前后 windows 的区域,重新写表输出
        BedOut.loc[:,'bedstart'] = BedOut.apply(ReCallPeaks.Addstart,args=(Windows_size,),axis=1)
        ## 强制限制区域为五倍滑窗长度,方便后面计算
        BedOut.loc[:,'bedend'] = BedOut['bedstart'] + Windows_size
        BedOut_result = BedOut[['Chr','bedstart','bedend']]
        BedOut_result.to_csv(bedfile,sep='\t',index=False,header=None)
        print (f'PID {a} is done! outfile is {outfile} {Top_out_1000_file} {bedfile} \n')

    '''
    输出包含每个峰区域位置的信息 '_Peaks.message.csv' ,包含 Coverage 倍数,Coverage 均值,对照组 Coverage
    'Chr','start','end',sample,Sample_Mean,'TB0135'
    -----------------------------------------------------------------------------------------------
    输出峰高前1000的区域 bed 文件
    'Chr','bedstart','bedend'
    -----------------------------------------------------------------------------------------------
    输出峰所在区域前后一定区域 '_Peaks.selectregion'
    'Chr','bedstart','bedend'
    '''
    ## 直接被FindPeaks调用,主要防止start出现负数
    @staticmethod
    def Addstart (Dataframe,Windows_size):
        if Dataframe['start'] < Windows_size:
            return Dataframe['start']
        else:
            return Dataframe['start'] - Windows_size

    ## 对于不同的表格列表,返回对应的基因列表与样本重复列表
    ## 接受dataframe 列名列表
    @staticmethod
    def Fliter_gene_base_sample (data_columns: list)-> dict:
        gene = []
        for column in data_columns:
            if column.startswith("LITCHI"):
                pattern = r'(^.*?)r\d'
                a = re.findall(pattern,column)[0]
                gene.append(a)
        set_gene = set(gene)
        gene = list(set_gene)
        samples = {}
        for column in data_columns:
            for i in gene:
                if column.startswith(i):
                    if i in samples:
                        samples[i].append(column)
                    else:
                        samples[i] = []
                        samples[i].append(column)
        return gene,samples

################################################################################################################
#        二次寻峰,在第一次按照滑窗大小找到峰所在区域后,按照 samtools depth 结果进行二次寻峰,寻找蜂最高的位点
#       在执行二次寻峰前,应该对每个基因生成 逐 bp 的 depth 信息文件,文件中包含所有重复
# shell:
############################################################
#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 30
#SBATCH -J depth
#SBATCH -o bamdepth.log
#SBATCH -e bamdepth.err
# WORKDIR=$SLURM_SUBMIT_DIR
# source ~/.bashrc
# export PATH=/home/xurong_fan/bin:$PATH
# for file in `ls *.bam|grep -v 'LITCHI005122'|sed 's/r[0-2].sorted.bam//g'|uniq`;do
# echo "samtools depth -a -H -l 150 -o ${file}.depth -@ 10 -b ${file}_Peaks.selectregion ${file}r1.sorted.bam ${file}r2.sorted.bam"
# done|parallel -j 3
################################################################################################################
    
    ## 读取测序深度文件,并对于每个基因深度信息文件,计算多个重复的均值
    @staticmethod
    def Read_depth_file(Sample_depth_file):
        depth = pd.read_csv(Sample_depth_file,sep='\t')
        columns = depth.columns
        pattern = r'(^.*r\d?).sorted.bam'
        rename_columns = []
        for name in columns:
            if name.startswith('#'):
                rename_columns.append('Chr')
            elif (name == 'POS'):
                rename_columns.append('POS')
            else:
                a = re.findall(pattern,name)[0]
                rename_columns.append(a)
        depth.columns = rename_columns
        gene_pattern = r'(^.*?).depth'
        gene = re.findall(gene_pattern ,Sample_depth_file)[0]
        Samples = [name for name in depth.columns if name.startswith('LITCHI')]
        if len(Samples) == 1:
            depth.loc[:,gene] = depth[Samples[0]]
        else:
            depth.loc[:,gene] = (depth[Samples[0]] + depth[Samples[1]])/2
        return gene,depth
    
    ## 二次寻峰函数,为多线程设计
    @staticmethod
    def Find_peaks_position(Sample_depth_file,Windows_size):
        ## 接受读取函数的返回值
        gene,Peaks_depth = ReCallPeaks.Read_depth_file(Sample_depth_file)
        Position = []
        outfile = gene + '_Peaks.region.bed'
        print (f'Outfile is {outfile}')
        ## 按照制定的Windows_size,读取测序深度信息文件
        for i in range(0,len(Peaks_depth),Windows_size):
            select_df = Peaks_depth.iloc[i:i+Windows_size,:]
            Compu_np = select_df[gene].to_numpy()
            ## 寻峰,为了简化,直接寻找最高峰
            Peaks_index,_ = scipy.signal.find_peaks(Compu_np,height=np.max(Compu_np))
            # ## 不符合正态分布,放弃
            # _,pvalue = scipy.stats.normaltest(Compu_np)
            # if (pvalue < 0.05):
            #     continue
            ## peaks 的情况太复杂，放弃
            ## 无峰情况,抛弃
            if (len(Peaks_index)==0):
                continue
            ## 多个高峰，取中间位点
            average = math.floor(np.mean(Peaks_index))  ## 取最高点平均值,向下取整
            bed = select_df.iloc[average][['Chr','POS']]
            Position.append(bed)
        Position_df = pd.DataFrame(Position)
        ## 前后延长 50 bp , 即为最终 peaks 区域
        Position_df.loc[:,'start'] = Position_df['POS'] - 50
        Position_df.loc[:,'end'] = Position_df['POS'] + 50
        Result = Position_df[['Chr','start','end']]
        Result.to_csv(outfile,sep='\t',index=False,header=None)
    '''
    输出包含每个峰顶前后100bp的信息 '_Peaks.region'
    'Chr','bedstart','bedend'
    '''
def Call_Peaks():
    parser = OptionParser()
    parser.add_option("--bedgraphDirectory",help="The directory where the Bedgraph file resides")
    parser.add_option("--binsize_length",default="100",help="Cut the sliding window length of the Bedgraph file")
    parser.add_option("--threads",default="1",help="The number of threads used by the program")
    parser.add_option("--control_prefix",help="Input a control file")
    parser.add_option("--outRessultcsv",default="ReCallPeaks_Results.csv",help="Output file name")
    parser.add_option("--peaks_height_min",default="20",help="Peaks height min")
    parser.add_option("--peaks_height_max",default="1000",help="Peaks height max")
    parser.add_option("--peaks_spacing",default="100",help="Peaks spacing [int] ")
    parser.add_option("--threshold",default="2",help="The multiples difference between the control group and the experimental group")
    parser.add_option("--C_peaks_length",default="1000",help="Candidate peaks length")
    options, args = parser.parse_args()
    threads = int(options.threads)
    control = options.control_prefix
    control_file = options.control_prefix + 'r1' +'.bedgraphfliter.csv'
    outRessultcsv = options.outRessultcsv
    height_min = float(options.peaks_height_min)
    height_max = float(options.peaks_height_max)
    Peaks_height = [height_min,height_max]
    spacing = float(options.peaks_spacing)
    threshold = float(options.threshold)
    binsize = int(options.binsize_length)
    Windows_size = int(options.C_peaks_length)

    os.chdir(options.bedgraphDirectory)
    BedgraphFile = glob.glob('*.bedgraph')


    if __name__ == '__main__':
        p = mp.Pool(processes=threads)
        for file in BedgraphFile:
            a = ReCallPeaks()
            p.apply_async(a.Fliter_Bedgraph,args=(binsize,file),callback=a.call_back,error_callback=a.err_call_back)
        p.close()
        p.join()

    Flitered_File = glob.glob('*bedgraphfliter.csv')
    Flitered_File.sort()
    data = pd.read_csv(control_file,dtype={'Chr':'str','start':'int','end':'int','score':'int'}).iloc[:,0:3]
    for file in Flitered_File:
        pattern = r'(^.*r\d?).bedgraphfliter'
        column_name = re.findall(pattern,file)[0]
        data.loc[:,column_name] = ReCallPeaks.Read_Fliter_File(file).astype('int')
    print ('step 3 done!')

    gene,samples = ReCallPeaks.Fliter_gene_base_sample(data.columns)
    print ('step 4 done!')

    Site = data.iloc[:,:3]
    Sitedata = data.iloc[:,3:]
    print ('step 5 done!')

    Sitedata.loc[:,control] = Sitedata.apply(ReCallPeaks.Blank_compute,axis=1)
    print ('step 6 done!')

    ReCallPeaks.Sample_compute(samples,Sitedata)
    print ('step 7 done!')

    Result = pd.concat([Site,Sitedata],axis=1)
    print ('step 8 done!')

    Result.to_csv(outRessultcsv,index=False)
    print ('step 9 done!')

    if __name__ == '__main__':
        p = mp.Pool(threads)
        for sample in samples:
            a = ReCallPeaks()
            p.apply_async(a.FindPeaks,(Result,sample,Peaks_height,spacing,threshold,Windows_size),callback=a.call_back,error_callback=a.err_call_back)
        p.close()
        p.join()
    
    print ('Samtools start...')
    order=f"for file in `ls *.bam|grep -v {control}|sed 's/r[0-2].sorted.bam//g'|uniq`;do " +\
    'echo "samtools depth -a -H -l 150 -o ${file}.depth -@ 1 -b ${file}_Peaks.select.region.bed ${file}r1.sorted.bam ${file}r2.sorted.bam";done'+\
    f"|parallel -j {threads}"
    os.system(order)

    if __name__ == '__main__':
        p = mp.Pool(threads)
        DepthFile = glob.glob('*depth')
        for file in DepthFile:
            a = ReCallPeaks()
            p.apply_async(a.Find_peaks_position,args=(file,Windows_size),callback=a.call_back,error_callback=a.err_call_back)
        p.close()
        p.join()

Call_Peaks()