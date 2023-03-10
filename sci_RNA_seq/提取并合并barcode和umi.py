import pyfastx
import pandas as pd
import gzip
import os
import time

# 读取barcode文件
df = pd.read_csv('/root/wangje/Project/Yin/Script/p5_barcode.txt', sep='\t')
barcode1_values = set(df.loc[:, 'barcode1'].values)
barcode2_values = set(df.loc[:, 'barcode2'].values)

def extract_umi_barcode(input_folder, sample_id, umi_seq, out_path):
    time_start = time.time()
    # 存放匹配barcode和umi的标志物
    id_markers = []
    barcode_umi = []
    
    # 文件的不同是可以通过后缀的序号去识别的
    read1 = os.path.join(input_folder, f"{sample_id}1_1.fq.gz")
    read2 = os.path.join(input_folder, f"{sample_id}1_2.fq.gz")
    # 创建输出文件的路径
    output_file = os.path.join(out_path, f"{sample_id}1_2.barcode.fastq.gz")
    f3 = gzip.open(output_file, 'wb')
    
    fa = pyfastx.Fastx(read1)
    for name, seq, comment in fa:
        if umi_seq in seq:
            umi_index = seq.find(umi_seq)
            UMI_seq = seq[umi_index: umi_index + 6 + 9]
            barcode1 = seq[umi_index - 10: umi_index]
            barcode2 = seq[umi_index + 6 + 8: umi_index + 6 + 18]
            if barcode1 in barcode1_values and barcode2 in barcode2_values:
                id_markers.append(name)
                barcode_umi.append(f"{barcode1}_{UMI_seq}_{barcode2}")
    
    # 根据Read1中得到的tag信息截取Read1中的序列
    id_markers_new = [x.replace('/1', '/2') for x in id_markers]
    
    # 生成name和barcoe_umi_barcode的字典
    barcode_dict = {id_markers_new[i]: barcode_umi[i] for i in range(len(id_markers_new))}
    
    for name, seq, comment in pyfastx.Fastx(read2):
        if name in barcode_dict:
            # print(name + barcode_dict[name])
            # print(seq)
            # print('+')
            # print(comment)
            f3.write(('@' + name + barcode_dict[name] + "\n").encode())
            f3.write((seq+'\n').encode())
            f3.write(('+' + '\n').encode())
            f3.write((comment+'\n').encode())
    time_end = time.time()
    run_time = time_end - time_start
    print("""
          ------------------------summary---------------------------
          umi和两个barcode均一致的序列: %s 
          耗时：%s s
          """ %(len(id_markers_new), run_time))
    
    # print(len(id_markers_new))

input_folder = "/root/wangje/Project/Yin/E100067938-中山大学孙逸仙纪念医院"
sample_id = "I230303_I230303-"
umi_seq = "CAGAGC"
out_path = "/root/wangje/Project/Yin/Data/extract_umi_barcode"
extract_umi_barcode(input_folder, sample_id, umi_seq, out_path)
