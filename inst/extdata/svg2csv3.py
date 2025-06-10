import xml.etree.ElementTree as ET
import csv
import argparse
import re

# 设置命令行参数解析器
parser = argparse.ArgumentParser(description='Extract <circle> elements from an SVG file and export to CSV.')
parser.add_argument('--input', required=True, help='Path to input SVG file')
parser.add_argument('--output', required=True, help='Path to output CSV file')

args = parser.parse_args()

# 解析 SVG 文件
tree = ET.parse(args.input)
root = tree.getroot()

# SVG 命名空间
ns = {'svg': 'http://www.w3.org/2000/svg'}

# 查找所有 <circle> 元素
circles = root.findall('.//svg:circle', ns)

# 提取并排序
def extract_id_number(cls):
    match = re.search(r'id_(\d+)', cls)
    return int(match.group(1)) if match else float('inf')

circle_data = []
for circle in circles:
    cx = circle.get('cx')
    cy = circle.get('cy')
    cls = circle.get('class')
    circle_data.append((cls, cx, cy))

# 按 id_ 后的数字排序
circle_data.sort(key=lambda x: extract_id_number(x[0]))

# 写入 CSV 文件，只保留 cx 和 cy 列
with open(args.output, mode='w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['cx', 'cy'])  # 只写入两列标题
    for _, cx, cy in circle_data:
        writer.writerow([cx, cy])

print(f"提取完成，结果已保存到 {args.output}")
