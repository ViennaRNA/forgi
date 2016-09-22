from __future__ import print_function, absolute_import, division, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 



import argparse
import random
import math
import sys
import warnings
from collections import defaultdict
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.dssr as ftud
import pandas as pd
import numpy as np
import scipy.stats
import copy

try:
    import readline
except:
    pass #readline not available

import matplotlib
matplotlib.use("TkAgg")

def av_angle(angles):
    # http://stackoverflow.com/a/491907/5069869
    try:
        return math.atan(sum(math.sin(x) for x in angles)/sum(math.cos(x) for x in angles))
    except ZeroDivisionError:
        return float("nan")

def bin_angular(data, num_bins = 100):
    bins = np.linspace(0, 2*math.pi, num_bins)
    groups = data.angle.groupby(pd.cut(data.angle, bins))
    return groups.count().values, bins

def angular_chisquare(data1, data2):
    binned1,_ = bin_angular(data1)
    binned2,_ = bin_angular(data2)
    mask=(binned2>0)
    binned1=binned1[mask]
    binned2=binned2[mask]
    all1 = sum(binned1)
    binned2 = binned2/sum(binned2)*all1
    #print (binned1, binned2)
    return scipy.stats.chisquare(binned1, binned2)
    
def show_circlehist(data, title):
    """
    http://stackoverflow.com/a/22568292/5069869
    """
    try:
        import matplotlib.pyplot as plt
    except:
        print ("Could not import matplotlib. Histogram not shown.", file=sys.stderr)
        return
    fig, ax = plt.subplots(2,2, subplot_kw=dict(projection='polar'))
    fig.suptitle(title)
    fig.text(0.05,0.9, "{} datapoints".format(len(data)) )
    values, bins = bin_angular(data)
    bars = ax[0,0].bar((bins[:-1]+bins[1:])/2, values, width=2*math.pi/100, linewidth=0.25)    
    maxc = max(values)

    for r, bar in zip(values, bars):
        bar.set_facecolor(plt.cm.jet(r / maxc))

    ax[0,0].set_title("All")
    ax[0,0].set_theta_direction(-1)
    ax[0,0].set_theta_zero_location("W")
    for ang_type, ax in [(2, ax[0,1]),(3, ax[1,0]),(4, ax[1,1])]:
        data_f = data[(data.angle_type == ang_type) | (data.angle_type == -ang_type)]
        values, bins = bin_angular(data_f)
        bars = ax.bar((bins[:-1]+bins[1:])/2, values, width=2*math.pi/100, linewidth=0.25)
        for r, bar in zip(values, bars):
            bar.set_facecolor(plt.cm.jet(r / maxc))
        ax.set_title("Angle_type {}".format(ang_type))
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location("W")

    plt.show(block=False)
def show(data):
    for a_type in [2,3,4]:
        st = pd.value_counts(data[(data.angle_type == a_type) | (data.angle_type == -a_type)]["is_stacking_dssr"])
        try: t = st[True]
        except: t = 0
        try: f = st[False]
        except: 
            f=0
        try:
            r=t/(t+f)
        except:
            r=float("nan")
        ang = av_angle(data[(data.angle_type == a_type) | (data.angle_type == -a_type)]["angle"])
        print("Angle type {}: {}/{} stack ({:.2f} %). average angle {}".format(a_type, t, t+f, r, ang))
    
def eq(data, key, value):
    return data[data[key]==value]
def gt(data, key, value):
    return data[data[key]>value]
def lt(data, key, value):
    return data[data[key]<value]
def ge(data, key, value):
    return data[data[key]>=value]
def le(data, key, value):
    return data[data[key]<=value]

def interactive_analysis(data):
    try:
        from gi.repository import Notify
    except: 
        notification=None
    else:
        Notify.init("coaxial_stacking.py")
        notification = Notify.Notification.new("Loading of Data complete.", "Interactive coaxial stacking analysis is ready")
        notification.show()
    print(data.columns.values)
    filtered_data = data
    history=[]
    stored={}
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    OUTPUT = '\033[92m' #expected output of commands
    OKBLUE = '\033[94m' 
    BOLD = '\033[1m'
    ops = {
        "==": eq,
        "=": eq,
        ">": gt,
        "<": lt,
        "<=": le,
        ">=": ge
    }
    try:
        while True:            
            if history:
                title = "<"+";".join(history)+">"
            else:
                title = "<no filters>"                
            print (BOLD+title+OKBLUE)
            show(filtered_data)
            print(ENDC)
            f = input("Please add a filter (HELP to show help): ") #imported from future
            try:
                notification.close()
            except:
                pass
            if f=="HELP":
                print(  OUTPUT+"* Use HELP to show this help\n"
                        "* Input filers like 'ml_length == 3'/ 'segment_length < 4'\n"
                        "  Note that the space is important!\n"
                        "  Valid keys are:\n"
                        "    ml_length       \t# of multiloop segemnts in the multiloop\n"
                        "    segment_length  \t# of nucleotides in the multiloop segment\n"
                        "    pseudoknot      \t0 or 1 (1 if multiloop segment is part of a pseudoknot)\n"
                        "    all_same        \t0 or 1 (1 if all segments of the ultiloop have the same number of nucleotides)\n"
                        "    min_nt_pos      \tThe lowest nucleotide number in any multiloop segment of the multiloop\n"
                        "    broken          \tIf the multiloop segment is broken in the standard forgi graph (=longest segment of multiloop)\n"
                        "    angle           \tThe angle (in rad) between the two helices connected by the ml-segment\n"
                        "    is_stacking_dssr\tWhether the two helices connected by this ml-segment stack according to DSSR\n"
                        "* Use R to reset all filters.\n"
                        "* Use S to show a plot.\n"
                        "* Use 'PRINT key' (e.g. 'PRINT ml_length') to print all values of this key in the current dataset.\n"
                        "* Use 'SAVE name' to save the current filters in memory under the given name\n"
                        "* Use 'LOAD name' to load filters that were previousely saved with 'SAVE name'.\n"
                        "* USE 'SHOW_SAVED' to show all saved sub-datasets\n"
                        "* Use 'COMPARE name1 name2' to compare to subsets of the data, previousely stored with 'SAVE name1' and 'SAVE name2'\n"+ENDC)
            elif f=="R":
                filtered_data = data
                history=[]
            elif f=="S":
                show_circlehist(filtered_data, title)
            elif f.startswith("SAVE"):
                try:
                    name = f.split()[1]
                    stored[name]=(copy.copy(history), filtered_data)
                except:
                    print (FAIL+"Could not save. (Type HELP for more info)"+ENDC)
            elif f.startswith("LOAD"):
                try:
                    name = f.split()[1]
                    history, filtered_data = stored[name]
                except:
                    print (FAIL+"Could not load. (Type HELP for more info)"+ENDC)
            elif f == "SHOW_SAVED":
                print(OUTPUT, end="")
                for k, v in stored.items():
                    if v[0]:
                        title = "<"+";".join(v[0])+">"
                    else:
                        title = "<no filters>"                
                    print ( k, "\t", title)
                print(ENDC, end="")
            elif f.startswith("COMPARE"):
                try:
                    _, name1, name2 = f.split()
                    data1 = stored[name1][1]
                    data2 = stored[name2][1]
                except KeyError as e:
                    print (FAIL+"Dataset {} was not saved.".format(e)+ENDC)
                except:
                    print(FAIL+"Could not compare. Type 'HELP' for more info."+ENDC)
                else:
                    print(OUTPUT, end="")
                    v, p = angular_chisquare(data1, data2) #Null hypothesis: The same distribution
                    if p<0.01:
                        print ("Datasets are not the same (p={})".format(p))
                    else:
                        print ("Datasets might be correlated (p={})".format(p))
                    print(ENDC, end="")

            elif f.startswith("PRINT"):
                try:
                    key = f.split()[1]
                    print (filtered_data[key])
                except:
                    print (FAIL+"Invalid key for command PRINT"+ENDC)
            else:
                try:
                    key, op, val = f.split()
                    if op not in ops.keys():
                        print(FAIL+"Invalid filter. Expecting Operator ==, =, >, <, <= or >="+ENDC)
                        continue
                        
                    fd=ops[op](filtered_data, key, float(val))
                except:
                    print(FAIL+"Filter not understood"+ENDC)
                else:
                    filtered_data = fd
                    history.append(f)
 
    except (KeyboardInterrupt, EOFError):
        print()
        
    
def update_ml_data(annot, data):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annot_stacks = annot.coaxial_stacks()
    cg = annot._cg
    multiloops, nucleotides = cg.find_multiloop_loops()
    for ml in multiloops:
        if cg.is_loop_pseudoknot(ml):
            pn = True
        else:
            pn = False
        elem_lengths = [cg.element_length(m1) for m1 in ml if m1[0]=="m" ]
        all_same = (len(set(elem_lengths))==1)
        try:
            min_nt_pos = min( x for m1 in ml for x in cg.defines[m1] if m1[0]=="m")
        except ValueError:
            min_nt_pos = float("nan")
        for m in ml:            
            if m[0]!="m": continue
            #loop = cg.shortest_bg_loop(m)
            #print(loop)abs(cg.get_angle_type(m)            
            
            data['pseudoknot'].append(pn)            
            data['all_same'].append(all_same)
            data['min_nt_pos'].append(min_nt_pos)

            try:
                data['angle_type'].append(abs(cg.get_angle_type(m)))
            except TypeError: #angle type is None
                data['broken'].append(True)
                conn = cg.connections(m)
                data['angle_type'].append(abs(cg.connection_type(m, conn)))
            else:
                data['broken'].append(False)
            angle = cg.get_bulge_angle_stats(m)[0]
            data['angle'].append(angle.get_angle())
            s1, s2 = cg.edges[m]
            if [s1,s2] in annot_stacks or [s2,s1] in annot_stacks:
                data['is_stacking_dssr'].append(True)
            else:
                data['is_stacking_dssr'].append(False)
            data['ml_length'].append(len([x for x in ml if x[0]=="m"]))
            data['segment_length'].append(cg.element_length(m))


def generateParser():
    parser=argparse.ArgumentParser( description="Report coaxial stacking.")
    parser.add_argument("files", type=str, nargs="+", help="One or more cg files that all have the same bulge graph!")
    parser.add_argument("--dssr-json", type=str, nargs="*", help="One or more json files generated by x3dna-dssr. They have to be in the same order as the cg files.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Do not be so verbose!!!")
    parser.add_argument("-l", "--per-loop", action="store_true", help="Print statistics per multiloop.")
    parser.add_argument("-i", "--interactive", action="store_true", help="In combination with -l: Enter interactive mode for data analysis.", default="Tyagi")
    parser.add_argument("-m", "--method", type=str, help="'CG' or 'Tyagi'. Method used for stacking detection in forgi.", default="Tyagi")
    return parser

parser = generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    if args.dssr_json and len(args.dssr_json)!=len(args.files):
        parser.error( '--dssr-json must have the same number of arguments as files.' )
        
    if args.per_loop:
        data = defaultdict(list)
        for i, filename in enumerate(args.files):
            cg = ftmc.CoarseGrainRNA(filename)
            try:
                annot = ftud.DSSRAnnotation(args.dssr_json[i], cg)
            except LookupError:
                continue;
            update_ml_data(annot, data)
        if not args.interactive:
            print(data)
        interactive_analysis(pd.DataFrame(data))
        sys.exit(0)
    forgi_count = 0
    forgi_not_stacking = 0
    both_count = 0
    dssr_count = 0
    for i, filename in enumerate(args.files):
        if not args.quiet: print("=== FILE ", filename, args.dssr_json[i], " ===")
        cg = ftmc.CoarseGrainRNA(filename)
        try:
            annot = ftud.DSSRAnnotation(args.dssr_json[i], cg)
            #assert "coaxStacks" in annot._dssr, "{}".format(annot._dssr)
        except LookupError:
            for d in cg.defines:
                if d[0] in "mi" and cg.is_stacking(d):
                    print (cg.connections(d), "stack along", d)
        else:
            #annot.compare_dotbracket()
            annot.basepair_stacking(args.method)
            continue
            forgi, dssr = annot.compare_coaxial_stack_annotation(args.method)
            both = forgi & dssr
            forgi = forgi - both
            dssr = dssr - both
            dssr_count+=len(dssr)
            forgi_count+=len(forgi)
            forgi_not_stacking += len(list(f for f in dssr if f.forgi=="not stacking"))
            both_count+=len(both)
            if not args.quiet: 
                print ("{} found by forgi, {} by dssr, {} by both".format(len(forgi), len(dssr), len(both)))
                for f in forgi:
                    print ("{} and {} stacking in forgi".format(f.stems[0], f.stems[1]))
                for d in dssr:
                    print ("{} and {} stacking in dssr (in forgi {})".format(d.stems[0], d.stems[1], d.forgi))
                for b in both:
                    print ("{} and {} stacking in both".format(b.stems[0], b.stems[1]))
    print("======= SUMMARY ======")
    if args.dssr_json:
        total_count = forgi_count + dssr_count + both_count
        if total_count:
            print ("forgi found {}, dssr {}, thereof both {} of all stacks (Forgi not stacking: {})".format(
                                                  (forgi_count+both_count),
                                                  (dssr_count+both_count),
                                                  (both_count),
                                                  (forgi_not_stacking)))
            print ("forgi found {}%, dssr {}%, thereof both {}% of all stacks (Forgi not stacking: {}%)".format(
                                                  int((forgi_count+both_count)/total_count*100),
                                                  int((dssr_count+both_count)/total_count*100),
                                                  int((both_count)/total_count*100),
                                                  int((forgi_not_stacking)/total_count*100)))
        else:
            print("No coaxial stacks found")
