
import re

def proc_attr(a):
    
    m = reattr.match(a)
    
    return m.groups() if m is not None else None

reattr = re.compile(r'([A-z0-9_]+) "(.+)"')

fname = 'ensembl/Mus_musculus.GRCm38.92.gtf'
outf  = 'mouse_biotypes.tsv'

result = []

with open(fname, 'r') as fp:
    
    for l in fp:
        
        l = l.split('\t')
        
        if len(l) < 7:
            continue
        
        attr = dict(
            filter(
                lambda a: a is not None,
                map(
                    proc_attr,
                    (s.strip() for s in l[8].split(';'))
                )
            )
        )
        
        result.append((attr['gene_name'], attr['gene_biotype']))

with open(outf, 'w') as fp:
    
    _ = fp.write('name\tbiotype\n')
    
    _ = fp.write(
        '\n'.join(
            '\t'.join(l) for l in result
        )
    )
