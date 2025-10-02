import gzip
import sys

gtf_file = sys.argv[1]
output_bed = sys.argv[2]

with gzip.open(gtf_file, 'rt') as f, open(output_bed, 'w') as out:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        chrom = fields[0]
        feature_type = fields[2]
        start = fields[3]
        end = fields[4]
        strand = fields[6]
        attributes = fields[8]
        
        if feature_type != 'gene':
            continue
        
        attr_dict = {}
        for attr in attributes.split(';'):
            attr = attr.strip()
            if not attr:
                continue
            key_value = attr.split(' ', 1)
            if len(key_value) == 2:
                key = key_value[0]
                value = key_value[1].strip('"')
                attr_dict[key] = value
        
        gene_id = attr_dict.get('gene_id', '')
        gene_name = attr_dict.get('gene_name', '')
        gene_type = attr_dict.get('gene_type', '')
        
        out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t{gene_name}\t{strand}\t{gene_type}\n")
