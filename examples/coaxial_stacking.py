# coding: utf-8
from __future__ import print_function, absolute_import, division, unicode_literals
from builtins import str
from builtins import zip
from builtins import next
from builtins import input
from builtins import map
from builtins import range
from builtins import object
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 



import argparse
import random
import math
import sys
import warnings
from collections import defaultdict, namedtuple, OrderedDict, Mapping
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.dssr as ftud
import pandas as pd
import numpy as np
import scipy.stats
import copy
import traceback
from pprint import pprint
import logging
logging.basicConfig()
log = logging.getLogger(__name__)
try:
    import readline
except:
    pass #readline not available
import itertools as it
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.pyplot import cm 
from matplotlib.ticker import MaxNLocator
matplotlib.use("TkAgg")
from sklearn.neighbors.kde import KernelDensity
from sklearn.grid_search import GridSearchCV
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import PolynomialFeatures
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering, DBSCAN
from sklearn.naive_bayes import GaussianNB
import functools
pd.set_option('display.max_colwidth',1000)
import os.path
try:
    import pymbar
    if pymbar.version.version.split(".")[0] != '3':
        log.warning("Wrong version of optional dependency pymbar: Found {}, need 3+.".format(pymbar.version.version))
        pymbar = None

except ImportError as e:
    log.info("Optional dependency pymbar could not be imported. Reason: '{}'".format(e))
    pymbar = None
    
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
    if settings["bandwidth"] is not None:
        bw = settings["bandwidth"]        
        kde = KernelDensity(kernel=settings.get("kernel", "tophat"), 
                            bandwidth=bw).fit(data.reshape(-1, 1))

    else:
        grid = GridSearchCV(KernelDensity(kernel=settings.get("kernel", "tophat")),
                  {'bandwidth': np.linspace(math.radians(2), math.radians(30), 40)},
                cv=min(20, len(data))) # 20-fold cross-validation
        try:
            grid.fit(data.reshape(-1, 1))
        except ValueError:
            return #Do not plot kde, if we do not have enough datapoints
        #print("Bandwidth = {}".format(grid.best_params_))
        kde = grid.best_estimator_
    ax2 = polar_twin(ax)
    ax2.plot(x, np.exp(kde.score_samples(x.reshape(-1,1))), label="kde", linewidth = settings.get("kde_linewidth",2), color = settings.get("kde_color", "red"))

    
def show_circlehist(data, title, subplots, settings):
    num_bins = settings["num_bins"]
    max_val = settings["ylim"]
    if subplots:
        fig, ax = plt.subplots(2,2, subplot_kw=dict(projection='polar'))
        mainAx = ax[0,0]
    else:
        fig, mainAx = plt.subplots(subplot_kw=dict(projection='polar'))
    fig.suptitle(title)
    fig.text(0.05,0.9, "{} datapoints".format(len(data)) )
    if len(data):
        if num_bins is not None:
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
        
        if settings["kernel"] is not None:
            plot_kde(data.angle, mainAx, settings)
            
        if subplots:
            for ang_type, ax in [(2, ax[0,1]),(3, ax[1,0]),(4, ax[1,1])]:
                data_f = data[(data.angle_type == ang_type) | (data.angle_type == -ang_type)]    
                if num_bins is not None:
                    values, bins = bin_angular(data_f, num_bins)
                    bars = ax.bar(bins[:-1], values, width=2*math.pi/num_bins, linewidth=0.25)
                for r, bar in zip(values, bars):
                    bar.set_facecolor(plt.cm.jet(r / maxc))
                ax.set_title("Angle_type {}".format(ang_type))
                ax.set_theta_direction(-1)
                ax.set_theta_zero_location("W")
                if max_val:
                    ax.set_ylim(0,max_val)
                if settings["kernel"] is not None:
                    plot_kde(data_f.angle, ax, settings)
    plt.show(block=False)

def show(data):
    for a_type in [2,3,4]:
        try:
            st = pd.value_counts(data[(data.angle_type == a_type) | (data.angle_type == -a_type)]["is_stacking_dssr"])
        except KeyError:
            t = 0
            f = len(data[(data.angle_type == a_type) | (data.angle_type == -a_type)])
        else:
            try: 
                t = st[True]
            except: 
                t = 0
            try: 
                f = st[False]
            except: 
                f=0
        try:
            r=t/(t+f)
        except:
            r=float("nan")
        ang = av_angle(data[(data.angle_type == a_type) | (data.angle_type == -a_type)]["angle"])
        print("Angle type {}: {}/{} stack ({:.2f} %). average angle {}".format(a_type, t, t+f, r, ang))

def decomposition(data):    
    colors = ['navy', 'turquoise', 'darkorange', 'red']
    # Preprocessing
    data = data.select_dtypes(include=[np.number,np.bool_]).dropna(axis = 1, how="all").dropna()    
    o1 = np.array(data.angle<np.pi/4)
    o2 = np.array((data.angle<np.pi/2) & (data.angle>=np.pi/4))
    o3 = np.array((data.angle<3*np.pi/4) & (data.angle>=np.pi/2))
    o4 = np.array((data.angle>=3*np.pi/4))
    labels = ["o1", "o2", "o3", "o4"]
    

    """
    dat1 = (data -data.min()) / (data.max() - data.min())
    data = data.dropna(axis = 1, how="all").dropna()
    data.angle*=np.pi
    # PCA
    pca = PCA(n_components=2)
    data_r = pca.fit(data).transform(data)
    print("FIRST COMPONENT")
    for c, w in zip(data.columns, pca.components_[0]):
        print(c, w)
    print("SECOND COMPONENT")
    for c, w in zip(data.columns, pca.components_[0]):
        print(c, w)
    for o, c, l in zip([o1,o2,o3,o4], colors, labels):
      plt.scatter(data_r[o,0], data_r[o,1], color = c, label=l)
    plt.legend()
    plt.show(block=False)
    plt.figure()
    k_means = DBSCAN(eps=0.25)#n_clusters = 3)
    k_means.fit(data_r)
    labels = k_means.labels_
    plt.scatter(data_r[:, 0], data_r[:, 1], c=labels.astype(np.float))

    plt.show(block=False)
    
    return"""
    # LDA    

    plt.figure()
    data_noA = data.drop('angle', 1)
    lda = LDA(n_components=2)
    bins = np.linspace(0, 1, 5)
    print(data.angle.describe())
    categories =  pd.cut(data.angle, 5)
    print(categories)
    X_r2 = lda.fit(data_noA, categories).transform(data_noA)
    print("F")
    for o, c, l in zip([o1,o2,o3,o4], colors, labels):
        plt.scatter(X_r2[o,0], X_r2[o,1], color = c, label=l)
    plt.legend()
    plt.show(block=False)
    
    gnb = GaussianNB()
    pred = gnb.fit(data_noA, categories).predict(data_noA)
    print(sum(pred==categories), "ALL", len(data_noA))
    
    return
    # ================== Regression =====================
    # Lasso
    plt.figure()
    lasso = Lasso()
    lasso.fit(data_noA, data.angle)
    prediction = lasso.predict(data_noA)
    plt.scatter(data.angle, prediction)
    plt.xlabel("observed angle")
    plt.ylabel("predicted angle")
    plt.title("LASSO")
    plt.show(block=False)
    # Elsatic Net
    plt.figure()
    en = ElasticNetCV()
    en.fit(data_noA, data.angle)
    prediction = en.predict(data_noA)
    plt.scatter(data.angle, prediction)
    plt.xlabel("observed angle")
    plt.ylabel("predicted angle")
    plt.title("Elastic Net")
    plt.show(block=False)
    # Polynomial Features
    poly = PolynomialFeatures(degree=3)
    data_noAPoly = poly.fit_transform(data_noA)
    plt.figure()
    en = ElasticNetCV()
    en.fit(data_noAPoly, data.angle)
    prediction = en.predict(data_noAPoly)
    plt.scatter(data.angle, prediction)
    plt.xlabel("observed angle")
    plt.ylabel("predicted angle")
    plt.title("Elastic Net Poly")
    plt.show(block=False)

 
def distribution_change(data, key, target_range):
    print(target_range)
    vals = defaultdict(list)
    averages = []
    significancy = []
    keys = []
    first = True
    x_values = []
    for i in target_range:
        try:
            len(i)
        except:
            data_f = data[data[key]==i]
            x_values.append(i)
        else:
            data_f = data[(data[key]>=i[0])&(data[key]<i[1])]
            x_values.append((i[0]+i[1])/2)

        averages.append(data_f.angle.mean())

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
    ax2.plot(x_values, significancy, "--", color="black", label = "num samples")
    for k in keys:
        ax[0].plot(x_values, vals[k], "o-", label = "{}".format(k))
    ax[0].legend(loc="upper center")
    ax2.legend(loc="upper right")
    ax2.set_ylabel("number of samples")
    ax[0].set_ylabel("Fraction")
    ax[0].set_xlabel(key.replace("_", " "))
    print(list(map(math.degrees, averages)))
    ax[1].plot(x_values, averages, label="Average angle")
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

def step_from_str(stri):
    if not "step" in stri:
        return float('nan')
    last = stri.split("step")[-1]
    steps = last.split(".coord")[0]
    return int(steps)

def datafilter(f):
    @functools.wraps(f)
    def wrapped(data, key, value):
        newdata = f(data, key, value)
        newdata.history = data.history + ["{} {} {}".format(key, f.__doc__, value)]
        newdata.full_dataset_name = data.full_dataset_name
        return newdata
    return wrapped
@datafilter
def eq(data, key, value):
    """=="""
    try:
        if np.isnan(value):
            return data[np.isnan(data[key])]
    except TypeError: #Datatype does not support nan (int/ string)
        pass 
    return data[data[key]==value]
@datafilter
def ne(data, key, value):
    """!="""
    try:
        if np.isnan(value):
            return data[not np.isnan(data[key])]
    except TypeError: #Datatype does not support nan (int/ string)
        pass
    return data[data[key]!=value]
@datafilter
def gt(data, key, value):
    """>"""
    return data[data[key]>value]
@datafilter
def lt(data, key, value):
    """<"""
    return data[data[key]<value]
@datafilter
def ge(data, key, value):
    """>="""
    return data[data[key]>=value]
@datafilter
def le(data, key, value):
    """<="""
    return data[data[key]<=value]
@datafilter
def in_(data, key, value):
    """contains"""
    return data[data[key].str.contains(value)]

@datafilter
def notin(data, key, value):
    """doesnot_contain"""
    return data[np.invert(data[key].str.contains(value))]

def subsample_pymbar(data):                            
    if pymbar:
        #We do not sort in-place, because ds might or might not be a view to the original data
        sorted_data = data.sort_values("step", inplace = False)
        sorted_data.full_dataset_name = data.full_dataset_name
        sorted_data.history = data.history
        t, g, n_eff = pymbar.timeseries.detectEquilibration(sorted_data.angle)
        sorted_data = ge(sorted_data, "step", t)
        indices = pymbar.timeseries.subsampleCorrelatedData(sorted_data.angle, g, conservative = True)
        subsampled =  sorted_data.iloc[indices]
        subsampled.history = sorted_data.history + ["subsampled ({})".format(math.ceil(g))]
        subsampled.full_dataset_name = sorted_data.full_dataset_name
        return subsampled
    else:
        raise InvalidInput("Requires pymbar")
def hist_to_title(history):
    if history:
        title = "<"+";".join(history)+">"
    else:
        title = "<no filters>"
    return title

FAIL = '\033[91m'
ENDC = '\033[0m'
OUTPUT = '\033[92m' #expected output of commands
OKBLUE = '\033[94m' 
BOLD = '\033[1m'

class UnknownCommand(LookupError):
    pass

class InvalidInput(ValueError):
    pass

class ConfigValues(Mapping):
    def __init__(self):
        self.values = OrderedDict()
        self.type = {}
        self.doc = {}
        self.none_allowed = {}
        self.restricted = {}
    def add_item(self, name, default, type_, doc, none_allowed = False, restricted = None):
        if none_allowed and default is None:
            self.values[name] = None
        else:
            self.values[name] = type_(default)
        self.type[name] = type_
        self.doc[name] = doc
        self.none_allowed[name] = none_allowed
        self.restricted[name] = restricted
    def __getitem__(self, key):
        return self.values[key]
    def __setitem__(self, key, val):
        if key not in self:
            raise ValueError("Can only update values")
        if (val == "None" or val is None) and self.none_allowed[key]:
            self.values[key] = None
        elif val is None and not self.none_allowed[key]:
            raise TypeError("None is not allowed for key {}".format(key))
        else:
            val = self.type[key](val)
            if self.restricted[key] and val not in self.restricted[key]:
                raise ValueError("Value for configuration variable {} must be one of {}".format(name, self.restricted[name]))
            self.values[key] = val
    def __len__(self):
        return len(self.values)
    def __iter__(self):
        return iter(self.values)
    def get_doc(self, key):
        typestr = "{}".format(self.type[key])
        if self.restricted[key]:
            typestr += ", one of {}".format(self.restricted[key])
        elif self.none_allowed[key]:
            typestr += " or None"
        return  ("{}\n"
                "\t\t[{}]\n"
                "\t\t{}".format(self.values[key], typestr, self.doc[key]))
    
class interactive_analysis(object):
    def __init__(self, data):
        #We add a helpful column to the dataframe.
        try:
            data["step"]=data.filename.apply(step_from_str)
        except AttributeError:
            pass
        data.history = []
        data.full_dataset_name = "main"
        self.data={"main":data}
        self.filtered_data = data

        #self.history = []
        self.stored = {}
        self.settings = ConfigValues()
        
        self.settings.add_item("kernel", "epanechnikov", str, 
                               "What kernel schoud be used for the kernel density estimates used"
                               "in plots? (None to disable kde plots)",
                               none_allowed = True ,
                               restricted = ['gaussian', 'tophat', 'epanechnikov', 'exponential', 
                                             'linear', 'cosine', None])
        self.settings.add_item("num_bins", 36, int, 
                               "How many bins should be used for histograms of angles? (None to disable histograms)", 
                               none_allowed = True)
        self.settings.add_item("ylim", None, int, 
                               "The maximal y value for histograms of angles. None for auto-detect.", 
                               none_allowed = True)
        self.settings.add_item("bandwidth", None, float, 
                               "The bandwidth used for the kdes (radians). None for auto-detect via cross-validation.", 
                               none_allowed = True)

        self.ops = {
            "==": eq,
            "=": eq,
            ">": gt,
            "<": lt,
            "<=": le,
            ">=": ge,
            "!=": ne,
            "contains": in_, 
            "doesnot_contain": notin}
        showS1 = functools.partial(self.show_circularKDE, one=True)
        showS1.__doc__ = "Like S, but show only one plot for all angle stats combined."
        self.allowed_commands = OrderedDict([
            ("HELP", self.show_help),
            ("R", self.reset),
            ("S", self.show_circularKDE),
            ("S1", showS1),
            ("PLOT_DC", self.plot_dc),
            ("SET", self.set_config),
            ("SAVE", self.save),
            ("LOAD", self.load),
            ("DEL", self.del_saved),
            ("SHOW_SAVED", self.show_saved),
            ("COMPARE", self.compare),
            ("PRINT", self.print_data),
            ("E", self.exec_),
            ("W",self.workflow),
            ("IMPORT_DATASET", self.import_dataset),
            ("SELECT_DATASET", self.select_dataset)
        ])
        
        
    def start(self):
        try:
            import gi
            gi.require_version('Notify', '0.7')
            from gi.repository import Notify
        except: 
            notification=None
        else:
            try:
                Notify.init("coaxial_stacking.py")
                notification = Notify.Notification.new("Loading of Data complete.", "Interactive coaxial stacking analysis is ready")
                notification.show()
            except: pass
        print("Fields are:", self.filtered_data.columns.values)
        while True:
            try:
                # Show the filtered data
                title = hist_to_title(self.filtered_data.history)
                print (BOLD+title+OKBLUE)
                show(self.filtered_data)
                print(ENDC)
                command = input("Please enter a command (HELP to show help): ") #input was imported from future
                commands = command.split(";")
                for comm in commands:
                    try:
                        print(OUTPUT, end="")
                        self._perform(comm)
                        print(ENDC, end="")
                    except UnknownCommand as e:
                        print(FAIL,"Command '{}' not understood".format(e),ENDC)
                    except InvalidInput as e:
                        print(FAIL,"Invalid input: {}".format(e),ENDC)
                    except Exception as e:
                        print(FAIL)
                        traceback.print_exc()
                        print("An exception of type {} occurred: {}".format(type(e), e),ENDC)
                        #raise
            except (KeyboardInterrupt, EOFError):
                try:
                    if input("Exit? y/N") in ["Y","y"]:
                        return
                except (KeyboardInterrupt, EOFError): #Second time
                    return
    def _perform(self, command):
        parts = command.strip().split()
        if not parts:
            return #Empty command
        if parts[0] in self.allowed_commands:
            self.allowed_commands[parts[0]](*parts[1:])
        elif parts[0] in self.filtered_data.columns.values:
            self._apply_filter(*parts)
        else:
            raise UnknownCommand(parts[0])

    def _apply_filter(self, *args):
        try:
            key, operator, value = args
        except ValueError:
            raise InvalidInput("Filters need to be specified as space-seperated triples: "
                               "key, operator, value. E.g. 'ml_length == 3'")
       
        d_type = self.filtered_data[key].dtype
        try:
            if (d_type == np.bool_ and value in ["0", "False", "false"] ):
                value = False
            else:
                value = d_type.type(value)
        except ValueError:
            raise InvalidInput("Cannot convert value {} to type {} of column {}.".format(value, d_type, key))
        if operator not in self.ops:
            raise InvalidInput("Operator {} not understood. Must be one of "
                               "{}".format(operator, ",".join(map(repr, self.ops.keys()))))
        #print(repr(value))
        self.filtered_data = self.ops[operator](self.filtered_data, key, value)

    def _get_range(self, *args):
        if len(args)==3:
            key, from_, to_ = args
            step = None
        elif len(args)==4:
            key, from_, to_, step = args
        else:
            raise InvalidInput("A range has to be specified by 3-4 parameters "
                               "(key, from, to, [step]), found {} of length {}".format(" ".join(map(str,args)), len(args)))
        use_ranges = False
        if "." in from_ or "." in to_ or (step and "." in step):
            use_ranges = True
        if self.filtered_data[key].dtype == np.float_:
            use_ranges = True
        if use_ranges:
            from_, to_ = sorted([float(to_),float(from_)])
            if not step:
                step =(from_-to_)/10
                if step<0:
                    step*=-1
            else:
                step=abs(float(step))
            target_range = []
            for i in it.count():
                target_range.append((step*i, step*(i+1)))
                if step*(i+1)>max(to_, from_):
                    break
        else:
            from_, to_ = sorted([int(to_),int(from_)])
            if step:
                step = abs(int(step))
                target_range=list(range(int(from_), int(to_), int(step)))
            else:
                target_range=list(range(int(from_), int(to_)))
        return target_range

    def _filter_from_r(self, key, r):
        """param r: an integer (for '==') or a tuple (start, stop)"""
        if isinstance(r, int):
            data = self.ops["=="](self.filtered_data, key, r)
        else:
            data = self.ops["<"](self.ops[">="](self.filtered_data, key, r[0]), key, r[1])
        return data

    def show_help(self, *args):
        """
        Show this help message
        """
        if args:
            warnings.warn(FAIL+"Arguments ignored: {}".format(args)+ENDC)
        print("Usage of the interactive prompt:\n"
                     "###  Input filers like 'ml_length == 3' or 'segment_length < 4'\n"
                     "     Note that the space is important!\n"
                     "     Valid keys depend on the loaded dataset. Currently they are:\n")
        for header in self.filtered_data.columns.values:
            print(   "          {}".format(header))
        print     ("###  In addition, the following commands are available:")
        max_commandlength = max(len(c) for c in self.allowed_commands.keys())
        template = "     * {: <"+str(max_commandlength)+"}\t{}"
        for command, function in self.allowed_commands.items():
            print (template.format(command, function.__doc__))
        print     ("###  The following configuration values are available:")
        template = "     * {: <"+str(max(len(c) for c in self.settings.keys()))+"}\t{})"
        for name in self.settings:
            print (template.format(name, self.settings.get_doc(name)))


    def reset(self, *args):
        """
        Reset all filters (restore original dataset).
        """
        self.filtered_data = self.data[self.filtered_data.full_dataset_name]

    def show_circularKDE(self, *args, **kwargs): #In python3: (self, *args, one=False)
        """
        Show a polar histogram and kernel density estimate of the angles found in the dataset.
        
        Use 'S' without arguments for a plot of the current dataset
        Use 'S for saved NAME1, NAME2' to show the plots for the saved datasets NAME1, NAME2
        Use 'S for KEY FROM TO [STEP]' to show the plots for subsets of the current dataset where
            the key KEY falls into different regions of the range.
        """
        #http://stackoverflow.com/a/22568292/5069869
        if "one" in kwargs and kwargs["one"]:
            subplots = False
        else:
            subplots = True
        if args and args[0] == "for":
            if len(args)>1 and args[1]=="saved":
                for dataname in args[2:]:
                    if dataname not in self.stored:
                        raise InvalidInput("The commands 'S for saved' and 'S1 for saved' expect a list of saved "
                                           "datasets, but nothing was saved under the name "
                                           "'{}'.".format(dataname))
                    show_circlehist(self.stored[dataname], hist_to_title(self.stored[dataname].history), 
                                    subplots, self.settings)
            elif len(args)>1 and args[1] in self.filtered_data.columns.values:
                range_ = self._get_range(*args[1:])
                for r in range_:
                    data = self._filter_from_r(args[1], r)
                    show_circlehist(data, hist_to_title(data.history), 
                                    subplots, self.settings)
            else:
                raise InvalidInput("The commands 'S for' and 'S1 for' expect either the keyword "
                                   "'saved' or a vaild key that identifies one column of the data,"
                                   " not {}".format(args[1]))
        elif args:
            raise InvalidInput("Use the commands S and S1 with no arguments or with 'for' as the first argument. Found {}".format(repr(args)))
        else:
            title = hist_to_title(self.filtered_data.history)
            show_circlehist(self.filtered_data, title, subplots, self.settings)

    def plot_dc(self, *args):
        """
        Plot the distribution of angles on different intervals as a function of the specified range.
        
        Use 'PLOT_DC [for] KEY FROM TO [STEP]'.
        """
        if args[:1] == ["for", "saved"]:
            raise InvalidInput("The PLOT_DC command does not support the 'for saved' syntax.")
        if args[0]=="for": #for is optional
            args=args[1:]
        key = args[0]
        if key not in self.filtered_data.columns.values:
            raise InvalidInput("Invalid key '{}' used in 'for key' syntax".format(key))
        range_ = self._get_range(*args)
        distribution_change(self.filtered_data, key, range_)
    def set_config(self, *args):
        """
        Set configuration values (see below)

        Use 'SET NAME VALUE' to set the configuration variable with name NAME to the value VALUE.
        """
        try:
            name, val = args
        except ValueError:
            raise InvalidInput("SET expects to be followed by the name of a config variable"
                               " and its desired value.")
        try:
            self.settings[name] = val
        except Exception as e:
            raise InvalidInput(e)
    def save(self, *args):
        """
        Save current dataset under the given name.

        Use 'SAVE NAME' to save the current dataset under the name NAME.
        Use 'SAVE NAME for KEY FROM TO [STEP]' to save subsets of the current dataset according
            to the range specified with FROM and TO. 
            They will have '_NUMBER' appended to their name.        
        """
        if len(args)>1 and args[1]=="for":
            range_ = self._get_range(*args[2:])
            for r in range_:
                data = self._filter_from_r(args[2], r)
                if isinstance(r, int):
                    name = "{}_{}".format(args[0], r) 
                else:
                    name = "{}_[{},{})".format(args[0], r[0], r[1]) 
                print("SAVING {}".format(name))
                self.stored[name]=data
        elif len(args)==0:
            raise InvalidInput("Need a name to save data")
        elif len(args)>1:
            raise InvalidInput("Name for data save must not contain any spaces! "
                               " (Use the 'SAVE NAME for'-syntax for saving several sub-datasets")
        else:
            if args[0]=="for":
                raise InvalidInput("Name 'for' is not allowed as a name for saving data, because it is a keyword")
            self.stored[args[0]]=self.filtered_data
    def load(self, *args):
        """ 
        Load a previousely saved dataset.
        """
        if len(args)==0:
            raise InvalidInput("Name of dataset to load was not specified.")
        elif len(args)>1:
            raise InvalidInput("Can load only one dataset at a time.")
        self.filtered_data = self.stored[args[0]]
    def del_saved(self, *args):
        """
        Delete a dataset that was previousely saved under this name(s)

        Use 'DEL NAME1 [NAME2...] to delete saves NAME1,... 
        """
        if len(args)==0:
            raise InvalidInput("Need at least one name of a save which should be deleted.")
        for name in args:
            del self.stored[name]
    def show_saved(self,*args):
        """
        Show all saved datasets
        """
        for k, v in sorted(self.stored.items()):
            title = hist_to_title(v.history)
            print ( k, "\t", title)
    
    @staticmethod
    def _compare_data(datasets):
        for data1, data2 in it.combinations(datasets, 2):
            s, p = scipy.stats.ks_2samp(data1.angle, data2.angle)
            if p<0.01:
                print ("'{:>24}' != '{:<24}' (p={:1.2e}, s={:1.2f})".format(
                                                        hist_to_title(data1.history)[-24:],
                                                        hist_to_title(data2.history)[-24:], p, s))
            else:
                print ("'{:>24}' ~~ '{:<24}' (p={:1.1e}, s={:1.1f})".format(
                                                        hist_to_title(data1.history)[-24:],
                                                        hist_to_title(data2.history)[-24:], p, s))
            
    def compare(self, *args):
        """
        Compare two saved datasets using the Kolmogorov-Smirnov test.

        Use 'COMPARE name1 name2' to compare two subsets of the data, previousely stored with 'SAVE name1' and 'SAVE name2'
        Use 'COMPARE name1 name2 name3 [name4...] to perform pairwise comparison of multiple datasets.
        Use 'COMPARE for 
        """            
        print ("Kolmogorov-Smirnov-Test:")
        datasets = []
        for name in args:
            datasets.append(self.stored[name])
            self._compare_data(datasets)
            
    def print_data(self, *args):
        """
        Print some values and a summary for a column of data.

        Use 'PRINT key' to print information about the column with name key
        """
        if len(args)!=1 or args[0] not in  self.filtered_data.columns.values:
            raise InvalidInput("Expecting one argument containing a valid key!")
        else:
            print (self.filtered_data[args[0]])
            print(self.filtered_data[args[0]].describe())
    def exec_(self, *args):
        """
        DANGEROUSE! Execute arbitrary python code
        """
        exec(" ".join(args))
    
    def import_dataset(self, *args):
        """
        Import another dataset into the application.
        
        Use 'IMPORT_DATASET path/to/csv/file/csv name' to import the file and assign the name 'name' to it.
        """
        if len(args)!=2:
            raise InvalidInput("Expecting 2 arguments: path and name. Found {}".format(args))
        data = pd.read_csv(os.path.expanduser(args[0]))
        if args[1] in self.data:
            raise InvalidInput("Name {} exists already. Please choose another name.")
        data.history = []
        data.full_dataset_name = args[1]
        self.data[args[1]] = data
        print ("Dataset has been imported to {}. Use SELECT_DATASET to select it.".format(args[1]))
    
    def select_dataset(self, *args):
        """
        Switch to another dataset.
        """
        self.filtered_data = self.data[args[0]]
        print("Dataset {} has been selected".format(args[0]))
    def _save_by_seed(self, data, prefix = ""):            
        names = []
        for seed in range(1, 10):
            ds = self.ops["contains"](data, "parent_dir", "seed{}".format(seed))
            if len(ds):
                fname = "_{}seed_{}".format(prefix, seed)
                self.stored[fname] = ds
                names.append(fname)
        return names
    
    def _bin_by_run(self, data):
        parent_d = self.filtered_data.parent_dir.unique()
        names = []
        for d in parent_d:            
            data_d = self.ops["=="](self.filtered_data, "parent_dir", d)
            fname = "_{}".format(d)
            self.stored[fname] = data_d
            names.append(fname)
        return names
    
    def workflow(self, *args):
        """
        High-level workflow commands that make certain assumptions about the data.

        """ 
        if not args:
            return
        if args[0].upper()=="Bin_By_Seed".upper():
            saved_names = self._save_by_seed(self.filtered_data)
            for fname in saved_names:
                print ("File {} saved".format(fname))
        elif args[0].upper() == "Bin_By_Run".upper():
            saved_names = self._bin_by_run(self.filtered_data)
            for fname in saved_names:
                print ("File {} saved".format(fname))
        elif args[0].upper()=="Sampling_Efficiency".upper():
            STEPS = [62,93,140,210,315,473, 710, 1065, 1598, 2397,3596,5394,8091,12137,18206,27309,40964]
            #Get names for different directories without the seed part
            parent_d = self.filtered_data.parent_dir.unique()
            runs = set()
            for d in parent_d:
                if "seed" in d:
                    runs.add(d.partition("seed")[0])
                else:
                    raise InvalidInput("Command Sampling_Efficiency does not work for current dataset."
                                       "Expecting 'seed' in directory name")
            color=iter(cm.rainbow(np.linspace(0,1,len(runs))))
            fig, ax = plt.subplots()                           
            for run in runs:                   
                #print("RUN {}".format(run))
                run_data = self.ops["contains"](self.filtered_data, "parent_dir", run)            
                plot_p = []
                plot_step = []
                for max_steps in STEPS:
                    #print("============= {} STEPS ==============".format(max_steps))
                    step_data = self.ops["<="](run_data,"step", max_steps)
                    seed_data = defaultdict(list)
                    for seed in range(1, 10):
                        ds = self.ops["contains"](step_data, "parent_dir", "seed{}".format(seed))
                        ds = self.ops["contains"](ds, "filename", "step")
                        elem_names = ds.element_name.unique()
                        for m in elem_names:
                            elem_data = self.ops["=="](ds, "element_name", m)
                            if len(elem_data):
                                try:
                                    seed_data[m].append(subsample_pymbar(elem_data))
                                except ValueError as e:
                                    log.exception(e)
                    p_values = []
                    for key, datasets in seed_data.items():
                        #print("+++ {} ({} datasets) +++".format(key, len(datasets)))
                        for data1, data2 in it.combinations(datasets, 2):
                            s, p = scipy.stats.ks_2samp(data1.angle, data2.angle)
                            #print(p)
                            p_values.append(p)
                    if p_values:
                        plot_p.append(p_values)
                        plot_step.append(max(step_data.step))
                plot_p = np.array(plot_p)
                print(plot_p, plot_p.shape, len(STEPS))
                c=next(color)
                if plot_step:
                    ax.plot(plot_step, plot_p, "o", c=c)
                    ax.plot(plot_step, np.mean(plot_p, axis = 1), "-", c=c, label = run)
            box = ax.get_position()       
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show(block = False)
        else:
            raise InvalidInput("Workflow step {} not defined".format(args[0]))

def update_ml_data(annot, data, pathname):    
    cg = annot._cg
    if not cg.defines:
        warnings.warn("Ignoring empty structure {}".format(cg.name))
        return
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annot_stacks = annot.coaxial_stacks()
    print(cg.name)
    multiloops = cg.find_mlonly_multiloops()
    path, filename = os.path.split(pathname)
    ppath, directory = os.path.split(path)
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
            """
            try:
                steric_r = cg.steric_value(ml, method = "r**-1")
                steric_r2 = cg.steric_value(ml, method = "r**-2")
                steric_r3 = cg.steric_value(ml, method = "r**-3")
                steric_c10 = cg.steric_value(ml, method = "cutoff 10")
                steric_c20 = cg.steric_value(ml, method = "cutoff 20")
                steric_c50 = cg.steric_value(ml, method = "cutoff 50")
                steric_c100 = cg.steric_value(ml, method = "cutoff 100")
                steric_r1e = cg.steric_value(ml, method = "r**-1e")
                steric_r2e = cg.steric_value(ml, method = "r**-2e")
                steric_r3e = cg.steric_value(ml, method = "r**-3e")
            except:
                print(pathname)
                return
            """
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
                """data['steric_r'].append(steric_r)
                data['steric_r2'].append(steric_r2)
                data['steric_r3'].append(steric_r3)
                data['steric_r1e'].append(steric_r1e)
                data['steric_r2e'].append(steric_r2e)
                data['steric_r3e'].append(steric_r3e)
                data['steric_c10'].append(steric_c10)
                data['steric_c20'].append(steric_c20)
                data['steric_c50'].append(steric_c50)
                data['steric_c100'].append(steric_c100)"""
                data['filename'].append(filename)
                data['parent_dir'].append(directory)
                data['complete_path'].append(pathname)
                data['element_name'].append(m)
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
                    data["segmentlength_left"].append(float("nan"))
                else:
                    data["neighbor_left"].append(abs(cg.connection_type(left_m,cg.connections(left_m))))
                    data["segmentlength_left"].append(cg.element_length(left_m))
                right_m = ml[(ml.index(m)+1)%len(ml)]
                if right_m[0]!="m":
                    data["neighbor_right"].append(0)
                    data["segmentlength_right"].append(float("nan"))
                else:
                    data["neighbor_right"].append(abs(cg.connection_type(right_m,cg.connections(right_m))))
                    data["segmentlength_right"].append(cg.element_length(right_m))


def generateParser():
    parser=argparse.ArgumentParser( description="Report coaxial stacking.")
    parser.add_argument("files", type=str, nargs="+", help="One or more cg files that all have the same bulge graph!")
    parser.add_argument("--dssr-json", type=str, nargs="*", help="One or more json files generated by x3dna-dssr. They have to be in the same order as the cg files.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Do not be so verbose!!!")
    parser.add_argument("-l", "--per-loop", action="store_true", help="Print statistics per multiloop.")
    parser.add_argument("-i", "--interactive", action="store_true", help="In combination with -l: Enter interactive mode for data analysis.")
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
                update_ml_data(annot, data, filename)
            if args.csv and len(data):
                pd.DataFrame(data).to_csv(args.csv, mode='a')
        if not args.interactive:
            print(data)
        else:
            try:
                ia = interactive_analysis(pd.DataFrame(data))
                print(type(ia.data))
                ia.start()
            finally:
                print('\033[0m')
            sys.exit(0)
    if args.dssr_json:
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
