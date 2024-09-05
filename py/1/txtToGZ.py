import gzip
import shutil

# 定义输入和输出文件名
input_file = 'cells_center.txt'
output_file = 'barcodes_pos.tsv.gz'

# 使用gzip模块进行压缩
with open(input_file, 'rb') as f_in:
    with gzip.open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

print(f"{input_file} 已成功压缩为 {output_file}")
