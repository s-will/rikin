#!/usr/bin/env python

import argparse

def dict_map(f,d):
    return dict(zip(d.keys(),map(f,d.values())))

def parse_row(names,row):
    xs=dict()
    for idx,x in enumerate(row.split()):
        xs[names[idx]]=x
    time = float(xs["Time"])
    del xs["Time"]
    xs = dict_map(int, xs)
    return (time,xs)

def unparse_row(row,names):
    (time,xs)=row
    time = "%e" % time
    return '\t'.join(map(str,[time]+[xs[k] for k in names]));

def read_table(fname):
    tab=list()
    names=list()
    with open(fname,'r') as f:
        names = f.readline().split()
        for row in f:
            xs = parse_row(names,row)
            tab.append(xs)
    names.remove("Time")
    return (tab,names)


def add_dict(x,y):
    res=dict(x)
    for k in y:
        if not k in x:
            res[k]=0
        res[k] = res[k] + y[k]
    return res

def sub_dict(x,y):
    res=dict(x)
    for k in y:
        if not k in x:
            res[k]=0
        res[k] = res[k] - y[k]
    return res

        
def diff_tab(tab,names):
    newtab=list()
    last=dict(zip(names,[0]*len(names)))
    for (time,row) in tab:
        diffrow=sub_dict(row,last)
        newtab.append((time,diffrow))
        last=row
    return newtab


def fst(xy):
    (x,y)=xy
    return x

def snd(xy):
    (x,y)=xy
    return y

def undiff_tab(tab,names):
    newtab=list()
    last=dict(zip(names,[0]*len(names)))
    lasttime=-1
    for (time,row) in tab:
        if time==lasttime:
            last=snd(newtab[len(newtab)-1])
            newtab=newtab[0:len(newtab)-2]
        newrow=add_dict(last,row)
        newtab.append((time,newrow))
        last=newrow
        lasttime=time
    return newtab


def merge_tab(left,right):
    
    result = []
    i, j = 0, 0
    while i < len(left) and j < len(right):
        if fst(left[i]) <= fst(right[j]):
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1

    result += left[i:]
    result += right[j:]
    
    return result

def normalize_tab(tab,n):
    newtab=list()
    for (time,row) in tab:
        row = dict_map(lambda x: x/float(n),row)
        newtab.append((time,row))
    return newtab

def main():
    description = 'Average over a set of Stoachastirator output files'
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument('infiles', nargs='+')
    
    args = parser.parse_args()

    n=len(args.infiles)
    
    total=None
    
    for f in args.infiles:
        (tab,names) = read_table(f)
        tab = diff_tab(tab,names)      
        if total==None:
            total=tab
        else:
            total = merge_tab(total,tab)

    total = undiff_tab(total,names)
    
    total = normalize_tab(total,n)
    
    print "Time","\t".join(names)
    for line in total:
        print unparse_row(line,names)

if __name__ == "__main__":
    main()

