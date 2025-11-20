#!/usr/bin/env python3

import sys
import os

if len(sys.argv) != 3:
    print("Usage: python 21.bionano-trans.py file1.txt file2.txt")
    sys.exit(1)

file1_path = sys.argv[1]
file2_path = sys.argv[2]

# derive output file names
file1_out = '21.'+ file1_path + ".processed"
file2_out = '21.'+ file2_path + ".marked"

# build mapping: file2 col1 -> col2
mapping = {}
file2_lines = []  # cache raw or parsed lines from file2

with open(file2_path, 'r') as f2:
    for line in f2:
        line = line.rstrip('\n')
        if line.startswith("#") or not line.strip():
            file2_lines.append((line, False))
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]
            file2_lines.append((parts, False))

# process file1 and replace column 2 using mapping
file1_col2_values = set()
new_file1_lines = []

with open(file1_path, 'r') as f1:
    for line in f1:
        line = line.rstrip('\n')
        if line.startswith("#") or not line.strip():
            new_file1_lines.append(line)
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            original_col2 = parts[1]
            file1_col2_values.add(original_col2)
            if original_col2 in mapping:
                parts[1] = mapping[original_col2]
        new_file1_lines.append('\t'.join(parts))

# mark file2 rows whose first column appears in file1 column 2
updated_file2_lines = []
for line_info, _ in file2_lines:
    if isinstance(line_info, str):
        updated_file2_lines.append(line_info)
    else:
        if line_info[0] in file1_col2_values:
            line_info.append("yes")
        updated_file2_lines.append('\t'.join(line_info))

# write outputs
with open(file1_out, 'w') as out1:
    for line in new_file1_lines:
        out1.write(line + '\n')

with open(file2_out, 'w') as out2:
    for line in updated_file2_lines:
        out2.write(line + '\n')

print(f"Done. Output files:\n - {file1_out}\n - {file2_out}")