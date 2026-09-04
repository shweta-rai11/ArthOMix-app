import re,gzip,sys,csv
def gsm_map(g):
    txt=open(f"{g}.gsm.txt").read()
    out={}
    for r in re.split(r'^\^SAMPLE = ',txt,flags=re.M)[1:]:
        gsm=r.split('\n',1)[0].strip()
        t=re.search(r'!Sample_title = (.*)',r).group(1).strip()
        ch=dict(x.split(': ',1) for x in re.findall(r'!Sample_characteristics_ch1 = (.*)',r))
        out[gsm]=(t,ch)
    return out
def run(g,infile,outfile,sep,keyfun):
    m=gsm_map(g); key2gsm={}
    for gsm,(t,ch) in m.items():
        k=keyfun(t,ch); key2gsm[k]=gsm
    with gzip.open(infile,'rt') as f, gzip.open(outfile,'wt') as o:
        hdr=next(f).rstrip('\n').split(sep)
        keep=[0]; names=['cpg']; miss=[]
        for i,h in enumerate(hdr[1:],1):
            if 'etection' in h: continue
            h2=h.strip('"')
            if h2 in key2gsm: keep.append(i); names.append(key2gsm[h2])
            else: miss.append(h2)
        print(g,"mapped",len(keep)-1,"of",len(m),"GSMs; unmapped columns:",miss[:10])
        o.write(','.join(names)+'\n'); n=0
        for line in f:
            p=line.rstrip('\n').split(sep)
            o.write(','.join(p[i].strip('"') for i in keep)+'\n'); n+=1
        print(g,"rows",n)
#run("GSE193879","GSE193879_Matrix_processed.csv.gz","GSE193879_betas_gsm.csv.gz",",",lambda t,ch:t)
run("GSE175758","GSE175758_GEO_processed.txt.gz","GSE175758_betas_gsm.csv.gz","\t",
    lambda t,ch: "S%s.c%s.%s.d%s"%(re.search(r'sample (\d+)',t).group(1),ch['course'],ch['cell.type'],ch['day']))
