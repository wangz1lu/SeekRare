import pandas as pd
import gzip
import os
import glob

def annotate_vcf_safe(vcf_path, tissue_dir, output_path):
    tissue_files = glob.glob(os.path.join(tissue_dir, "*.eQTLs.signif_pairs.parquet"))
    
    # 1. 获取后缀
    sample_id = pd.read_parquet(tissue_files[0], columns=['variant_id']).iloc[0]['variant_id']
    suffix = "_" + sample_id.split('_')[-1]
    
    # 2. 读取 VCF 位点并转为 Set (利用哈希查找，极快且省内存)
    print(f"正在读取 VCF 并构建位点集...")
    vcf_set = set()
    _open = gzip.open if vcf_path.endswith('.gz') else open
    with _open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            vid = f"chr{cols[0].replace('chr','')}_{cols[1]}_{cols[3]}_{cols[4]}{suffix}"
            vcf_set.add(vid)
    print(f"VCF 中共有 {len(vcf_set)} 个唯一变异位点。")

    # 3. 逐个组织处理并立即写出
    print("开始逐个组织筛选...")
    first_write = True
    
    for tissue_file in tissue_files:
        t_name = os.path.basename(tissue_file).split('.v11')[0]
        
        # 仅读取 GTEx 中在我们 VCF 里的行
        # 我们分块读取或直接过滤，避免一次性加载过多
        print(f"--> 处理组织: {t_name}")
        gtex_chunk = pd.read_parquet(tissue_file, columns=['variant_id', 'phenotype_id', 'slope', 'pval_nominal'])
        
        # 核心过滤：只保留存在于 VCF 集合中的行
        matched = gtex_chunk[gtex_chunk['variant_id'].isin(vcf_set)].copy()
        matched['tissue'] = t_name
        
        if not matched.empty:
            # 写入模式：第一次覆盖，后续追加
            mode = 'w' if first_write else 'a'
            header = True if first_write else False
            matched.to_csv(output_path, sep='\t', index=False, mode=mode, header=header)
            first_write = False
            print(f"   匹配到 {len(matched)} 条记录并已写入。")
        
        del gtex_chunk, matched

    print(f"\n分析完成！结果保存至: {output_path}")

# 执行
annotate_vcf_safe("/mnt/wzl/3.recessive.nocommon.vcf.gz", "Optic_System_Related/", "4.recessive_nocommon_GTEx_Long_Format.tsv")
