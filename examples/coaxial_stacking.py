# coding: utf-8
from __future__ import print_function, absolute_import, division, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 



import argparse
import random
import math
import sys
import warnings
from collections import defaultdict, namedtuple
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.dssr as ftud
import pandas as pd
import numpy as np
import scipy.stats
import copy
from pprint import pprint
import logging
try:
    import readline
except:
    pass #readline not available
import itertools as it
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MaxNLocator
matplotlib.use("TkAgg")
from sklearn.neighbors.kde import KernelDensity
from sklearn.grid_search import GridSearchCV

def polar_twin(ax):
    #http://stackoverflow.com/a/19620861/5069869
    ax2 = ax.figure.add_axes(ax.get_position(), projection='polar', 
                             label='twin', frameon=False,
                             theta_direction=ax.get_theta_direction(),
                             theta_offset=ax.get_theta_offset())
    ax2.xaxis.set_visible(False)
    # There should be a method for this, but there isn't... Pull request?
    ax2._r_label_position._t = (22.5 + 180, 0.0)
    ax2._r_label_position.invalidate()
    # Ensure that original axes tick labels are on top of plots in twinned axes
    for label in ax.get_yticklabels():
        ax.figure.texts.append(label)
    plt.setp(ax2.get_yticklabels(), color='red')
    ax2.get_xaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])

    return ax2

def av_angle(angles):
    # http://stackoverflow.com/a/491907/5069869
    try:
        return math.atan(sum(math.sin(x) for x in angles)/sum(math.cos(x) for x in angles))
    except ZeroDivisionError:
        return float("nan")

def bin_angular(data, num_bins = 100):
    bins = np.linspace(0, 2*math.pi, num_bins+1)
    groups = pd.value_counts(pd.cut(data.angle, bins), sort=False)
    return groups.values, bins

def angular_chisquare(data1, data2):
    binned1,_ = bin_angular(data1)
    binned2,_ = bin_angular(data2)
    mask=(binned2>0)
    binned1=binned1[mask]
    binned2=binned2[mask]
    binned2 = binned2/sum(binned2)*sum(binned1)
    #print (binned1, binned2)
    return scipy.stats.chisquare(binned1, binned2)

def plot_kde(data, ax, settings):
    x = np.linspace(0, np.pi, 200)
    if "bandwidth" in settings and settings["bandwidth"] is not None:
        bw = settings["bandwidth"]        
        kde = KernelDensity(kernel=settings.get("kernel", "tophat"), 
                            bandwidth=bw).fit(data.reshape(-1, 1))

    else:
        grid = GridSearchCV(KernelDensity(kernel=settings.get("kernel", "tophat")),
                  {'bandwidth': np.linspace(math.radians(2), math.radians(30), 40)},
                cv=min(20, len(data))) # 20-fold cross-validation
        grid.fit(data.reshape(-1, 1))
        #print("Bandwidth = {}".format(grid.best_params_))
        kde = grid.best_estimator_
    ax2 = polar_twin(ax)
    ax2.plot(x, np.exp(kde.score_samples(x.reshape(-1,1))), label="kde", linewidth = settings.get("kde_linewidth",2), color = settings.get("kde_color", "red"))

    
def show_circlehist(data, title, num_bins = 100, subplots=True, max_val=None, settings = {}):
    """
    http://stackoverflow.com/a/22568292/5069869
    """
    if subplots:
        fig, ax = plt.subplots(2,2, subplot_kw=dict(projection='polar'))
        fig.suptitle(title)
        mainAx = ax[0,0]
    else:
        fig, mainAx = plt.subplots(subplot_kw=dict(projection='polar'))
        fig.suptitle(title)

    fig.text(0.05,0.9, "{} datapoints".format(len(data)) )
    values, bins = bin_angular(data, num_bins)
    bars = mainAx.bar(bins[:-1], values, width=2*math.pi/num_bins, linewidth=0.25)    
    maxc = max(values)
    if max_val:
        mainAx.set_ylim(0,max_val)
    for r, bar in zip(values, bars):
        bar.set_facecolor(plt.cm.jet(r / maxc))

    if subplots:
        mainAx.set_title("All")
    mainAx.set_theta_direction(-1)
    mainAx.set_theta_zero_location("W")    
    
    if settings["show_kde"]:
        plot_kde(data.angle, mainAx, settings)
        
    if subplots:
        for ang_type, ax in [(2, ax[0,1]),(3, ax[1,0]),(4, ax[1,1])]:
            data_f = data[(data.angle_type == ang_type) | (data.angle_type == -ang_type)]
            values, bins = bin_angular(data_f, num_bins)
            bars = ax.bar(bins[:-1], values, width=2*math.pi/num_bins, linewidth=0.25)
            for r, bar in zip(values, bars):
                bar.set_facecolor(plt.cm.jet(r / maxc))
            ax.set_title("Angle_type {}".format(ang_type))
            ax.set_theta_direction(-1)
            ax.set_theta_zero_location("W")
            if max_val:
                ax.set_ylim(0,max_val)
            if settings["show_kde"]:
                plot_kde(data_f.angle, ax, settings)


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

def distribution_change(data, key, target_range):
    vals = defaultdict(list)
    averages = []
    significancy = []
    keys = []
    first = True
    for i in target_range:
        data_f = data[data[key]==i]
        averages.append(data_f.angle.mean())
        print(data_f.angle.mean())
        val, bins = bin_angular(data_f, 8)
        for j in range(4):
            k = "{:d}° - {:d}°".format(int(math.degrees(bins[j])), int(math.degrees(bins[j+1])))
            if first:
                keys.append(k)
            vals[k].append(val[j]/len(data_f))
        significancy.append(len(data_f))
        first = False
    fig, ax = plt.subplots(2)
    ax2 = ax[1].twinx()
    ax2.plot(target_range, significancy, "--", color="black", label = "num samples")
    for k in keys:
        ax[0].plot(target_range, vals[k], "o-", label = "{}".format(k))
    ax[0].legend(loc="upper center")
    ax2.legend(loc="upper right")
    ax2.set_ylabel("number of samples")
    ax[0].set_ylabel("Fraction")
    ax[0].set_xlabel(key.replace("_", " "))
    print(list(map(math.degrees, averages)))
    ax[1].plot(target_range, averages, label="Average angle")
    #ax[0].set_title("Fraction of junctions with certain angles")
    ax[1].legend(loc="lower left")
    ax[1].set_ylabel("Radians")
    ax[1].set_xlabel(key.replace("_", " "))
    
    # y-ticks in pi multiples from nye7 via http://stackoverflow.com/a/10731637/5069869
    unit   = 0.25
    y_tick = np.arange(0, 1+unit, unit)
    y_label = ["$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$", r"$\frac{3\pi}{4}$", r"$\pi$"]
    ax[1].set_yticks(y_tick*np.pi)
    ax[1].set_yticklabels(y_label, fontsize=15)
    plt.show(block = False)

def eq(data, key, value):
    try:
        if np.isnan(value):
            return data[np.isnan(data[key])]
    except Exception as e:
        print(e)
        raise
    return data[data[key]==value]
def ne(data, key, value):
    try:
        if np.isnan(value):
            return data[not np.isnan(data[key])]
    except Exception as e:
        print(e)
        raise
    return data[data[key]!=value]
def gt(data, key, value):
    return data[data[key]>value]
def lt(data, key, value):
    return data[data[key]<value]
def ge(data, key, value):
    return data[data[key]>=value]
def le(data, key, value):
    return data[data[key]<=value]

def hist_to_title(history):
    if history:
        title = "<"+";".join(history)+">"
    else:
        title = "<no filters>"
    return title

def interactive_analysis(data):
    def save(name):
        stored[name]=(copy.copy(history), filtered_data)
        
    try:
        from gi.repository import Notify
    except: 
        notification=None
    else:
        try:
            Notify.init("coaxial_stacking.py")
            notification = Notify.Notification.new("Loading of Data complete.", "Interactive coaxial stacking analysis is ready")
            notification.show()
        except: pass
    print(data.columns.values)
    filtered_data = data
    history=[]
    stored={}
    settings = {"num_bins":100, "ylim":None, "show_kde": True, "kernel": "epanechnikov"}
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
        ">=": ge,
        "!=": ne
    }
    while True:
        try:
            while True: 
                title = hist_to_title(history)
                print (BOLD+title+OKBLUE)
                show(filtered_data)
                print(ENDC)
                fs = input("Please add a filter (HELP to show help): ") #imported from future
                commands = fs.split(";")
                for f in commands:
                    f=f.strip()
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
                                "* Use 'COMPARE name1 name2' to compare to subsets of the data, previousely stored with 'SAVE name1' and 'SAVE name2'\n"
                                "* Use 'E pythoncode' to do an exec(pythoncode)\n"
                                "* Use 'SET num_bins 10' to set the number of bins for the histograms shown with 'S'\n"+ENDC)
                    elif f=="R":
                        filtered_data = data
                        history=[]
                    elif f=="S":
                        try:
                            show_circlehist(filtered_data, title, settings["num_bins"], max_val=settings["ylim"], settings = settings)
                        except Exception as e:
                            print(FAIL,e,ENDC)
                    elif f=="S1":
                        try:
                            show_circlehist(filtered_data, title, settings["num_bins"], subplots = False, max_val=settings["ylim"], settings = settings)
                        except Exception as e:
                            print(FAIL,e,ENDC)
                    elif f.startswith("S1 for"):
                        try:
                            to_load = f.split()[2:]
                            for dataname in to_load:
                                 show_circlehist(stored[dataname][1], hist_to_title(stored[dataname][0]), 
                                                 settings["num_bins"], subplots = False, 
                                                 max_val=settings["ylim"], settings = settings)
                        except Exception as e:
                            print(FAIL,e,ENDC)                        
                    elif f.startswith("PLOT_DC"):
                        try:
                            _, key, from_, to_ = f.split()
                            distribution_change(filtered_data, key, range(int(from_), int(to_)))
                        except Exception as e:
                            print (FAIL+"Error: {}".format(e)+ENDC)

                    elif f.startswith("SET"):
                        try:
                            _,key,to_set = f.split()
                            if key in ["num_bins", "bandwidth"]:
                                settings[key] = int(to_set)
                            elif key == "kernel":
                                if to_set not in ['gaussian', 'tophat', 'epanechnikov', 'exponential', 
                                                  'linear', 'cosine']:
                                    print (FAIL+"Kernel {} not supported".format(to_set)+ENDC)
                                else:
                                    settings[key] = to_set
                            elif key == "ylim":
                                if to_set == "None":
                                    settings["ylim"] = None
                                else:
                                    settings["ylim"] = int(to_set)
                            else:
                                print (FAIL+"Key {} not understood".format(key)+ENDC)
                        except Exception as e:
                            print (FAIL+"Could not change settings: {}".format(e)+ENDC)
                    elif f.startswith("SAVE"):
                        try:
                            fields = f.split()[1:]
                            if len(fields)==1: 
                                name = fields[0]
                                save(name)
                            elif fields[1] != "for":
                                print (FAIL+"Could not save. (Type HELP for more info)"+ENDC)
                            else:
                                name = fields[0]
                                key = fields[2]
                                from_ = int(fields[3])
                                to_ = int(fields[4])
                                f_data2 = filtered_data
                                for i in range(from_, to_):
                                    fd=ops["=="](f_data2, key, i)
                                    h_new = copy.copy(history)
                                    h_new.append("{} == {}".format(key, i))
                                    print("SAVING {}".format(name+"_"+str(i)))
                                    stored[name+"_"+str(i)]=(h_new, fd)
                        except:
                            print (FAIL+"Could not save. (Type HELP for more info)"+ENDC)
                    elif f.startswith("LOAD"):
                        try:
                            name = f.split()[1]
                            history, filtered_data = stored[name]
                            history = copy.copy(history)
                        except:
                            print (FAIL+"Could not load. (Type HELP for more info)"+ENDC)
                    elif f == "SHOW_SAVED":
                        print(OUTPUT, end="")
                        for k, v in sorted(stored.items()):
                            if v[0]:
                                title = "<"+";".join(v[0])+">"
                            else:
                                title = "<no filters>"                
                            print ( k, "\t", title)
                        print(ENDC, end="")
                    elif f=="COMPARE ALL":
                        for s1, s2 in it.combinations(stored.items, 2):
                            print ("Kolmogorov-Smirnov")
                            s, p = scipy.stats.ks_2samp(s1[1][1].angle, s2[1][1].angle)
                            if p<0.01:
                                print (s1[0], " IS DIFFERENT TO ", s2[0])
                            else:
                                print (s1[0], " ~ ", s2[0])

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
                                print ("Chi Square: Datasets are not the same (p={})!!!".format(p))
                            else:
                                print ("Chi Square: Datasets might be correlated (p={})".format(p))

                            s, p = scipy.stats.ks_2samp(data1.angle, data2.angle)
                            if p<0.01:
                                print ("Kolmogorov-Smirnov: Datasets are not the same (p={})!!!".format(p))
                            else:
                                print ("Kolmogorov-Smirnov: Datasets might be correlated (p={})".format(p))

                            s, c, p = scipy.stats.anderson_ksamp([data1.angle, data2.angle])
                            #print(s,c,p)
                            if s>c[-1]:
                                print ("Anderson-Darling: Datasets are not the same (p={})!!!".format(p))
                            else:
                                print ("Anderson-Darling: Datasets might be correlated (p={})".format(p))

                            print(ENDC, end="")

                    elif f.startswith("PRINT"):
                        try:
                            key = f.split()[1]
                            print (filtered_data[key])
                        except:
                            print (FAIL+"Invalid key for command PRINT"+ENDC)
                    elif f.startswith("E"):
                        try:
                            exec f[2:] in globals(), locals()
                        except Exception as e:
                            print (FAIL,type(e), e, "Cannot eval command '{}'".format(f[2:])+ENDC)

                    else:
                        try:
                            key, op, val = f.split()
                            if op not in ops.keys():
                                print(FAIL+"Invalid filter. Expecting Operator ==, =, >, <, <= or >="+ENDC)
                                continue

                            fd=ops[op](filtered_data, key, float(val))
                        except:
                            print(FAIL+"Filter '{}' not understood".format(f)+ENDC)
                        else:
                            filtered_data = fd
                            history.append(f)

        except (KeyboardInterrupt, EOFError):
            print()
            while True:
                if input("Exit? Y/N")=="Y":
                    exit = 1
                    break
                else:
                    exit = 0
                    break
            if exit:
                break
        
    
def update_ml_data(annot, data):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annot_stacks = annot.coaxial_stacks()
    cg = annot._cg
    print(cg.name)
    multiloops = cg.find_mlonly_multiloops()
    for ml in multiloops:
            descriptors = cg.describe_multiloop(ml)
            pk = "pseudoknot" in descriptors
            regular_ml  = "regular_multiloop" in descriptors
            elem_lengths = [cg.element_length(m1) for m1 in ml ]
            all_same = (len(set(elem_lengths))==1)
            nuc_len = sum(cg.element_length(m1) for m1 in ml)
            try:
                min_nt_pos = min( x for m1 in ml for x in cg.defines[m1] )
            except ValueError as e:
                min_nt_pos = float("nan") #All 0-length segments
            for m in ml:      
                #if m in known_segments:
                #    continue
                #known_segments.add(m)
                if m[0]!="m": continue
                #loop = cg.shortest_bg_loop(m)
                #print(loop)abs(cg.get_angle_type(m)            
                data['regular_ml'].append(regular_ml) 
                data['pseudoknot'].append(pk) 
                data['open'].append("open" in descriptors) 
                data['all_same'].append(all_same)
                data['min_nt_pos'].append(min_nt_pos)
                data['pdb_name'].append(cg.name)
                data['nucleotide_length'].append(nuc_len)
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
                left_m = ml[ml.index(m)-1]
                if left_m[0]!="m":
                   data["neighbor_left"].append(0)
                else:
                  data["neighbor_left"].append(abs(cg.connection_type(left_m,cg.connections(left_m))))
                right_m = ml[(ml.index(m)+1)%len(ml)]
                if right_m[0]!="m":
                   data["neighbor_right"].append(0)
                else:
                  data["neighbor_right"].append(abs(cg.connection_type(right_m,cg.connections(right_m))))

    

def generateParser():
    parser=argparse.ArgumentParser( description="Report coaxial stacking.")
    parser.add_argument("files", type=str, nargs="+", help="One or more cg files that all have the same bulge graph!")
    parser.add_argument("--dssr-json", type=str, nargs="*", help="One or more json files generated by x3dna-dssr. They have to be in the same order as the cg files.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Do not be so verbose!!!")
    parser.add_argument("-l", "--per-loop", action="store_true", help="Print statistics per multiloop.")
    parser.add_argument("-i", "--interactive", action="store_true", help="In combination with -l: Enter interactive mode for data analysis.", default="Tyagi")
    parser.add_argument("-m", "--method", type=str, help="'CG' or 'Tyagi'. Method used for stacking detection in forgi.", default="Tyagi")
    parser.add_argument("-v", "--verbose", action="store_true", help="Be verbose")
    parser.add_argument("--csv", type=str, help="Store data in csv with this filename", default="")
    parser.add_argument("-c", "--continue-from", type=str, help="Load Data from csv with this filename", default="")

    return parser

class DummyAnnotation(object):
    def __init__(self, cg):
        self._cg = cg
    def coaxial_stacks(self):
        return []
parser = generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
        #logging.getLogger('forgi.threedee').setLevel(logging.ERROR)

    if args.dssr_json and len(args.dssr_json)!=len(args.files):
        parser.error( '--dssr-json must have the same number of arguments as files.' )
        
    if args.per_loop:
        if args.continue_from:
            data = pd.read_csv(args.continue_from)
        else:
            data = defaultdict(list)
            for i, filename in enumerate(args.files):
                cg = ftmc.CoarseGrainRNA(filename)
                try:
                    annot = ftud.DSSRAnnotation(args.dssr_json[i], cg)
                except (LookupError, TypeError):
                    annot = DummyAnnotation(cg)
                update_ml_data(annot, data)
            if args.csv:
                pd.DataFrame(data).to_csv(args.csv)
        if not args.interactive:
            print(data)
        try:
            interactive_analysis(pd.DataFrame(data))
        finally:
            print('\033[0m')
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
