import os
infile = r"c:\Users\sarachchige\OneDrive - UNSW\AIM3\Sample INPUTs\Patristic_distance_calculator\aln_std_nt.fasta"
outfile = infile + ".tmp"
keep_prefixes = (">13_", ">38_", ">42_")

kept_headers = []
count_headers = 0
count_kept = 0

with open(infile, 'r', encoding='utf-8') as inf, open(outfile, 'w', encoding='utf-8') as outf:
    write = False
    for line in inf:
        if line.startswith('>'):
            count_headers += 1
            write = any(line.startswith(p) for p in keep_prefixes)
            if write:
                kept_headers.append(line.strip())
                count_kept += 1
        if write:
            outf.write(line)

# Replace original file with filtered file
os.replace(outfile, infile)

print(f"Total headers in original file: {count_headers}")
print(f"Kept headers: {count_kept}")
for h in kept_headers:
    print(h)
