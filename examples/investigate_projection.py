#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

__metaclass__ = type #New style classes in Python 2.x

import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as fpp
import forgi.projection.hausdorff as fph
import forgi.threedee.utilities.vector as ftuv
import Tkinter as tk
import ttk
from tkFileDialog import askopenfilename, asksaveasfilename
import tkMessageBox
import sys, argparse, math
from PIL import Image, ImageTk
import numpy as np
import scipy

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('cgfile', nargs="*", help='One *.cg/*.coord files holding the RNA.')
    parser.add_argument('-w', '--width', type=int ,help='Width in Angstrom', default = 200)
    parser.add_argument('-d', '--dpi', type=int, help='Number of pixels in the image', default = 30)
    return parser

class RnaDisplay(object):
    def __init__(self, args):
        self.root = tk.Tk()
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.root.title("Manual Projection Rasterization")
        mainframe = ttk.Frame(self.root, padding="3 3 12 12")
        self.mainframe = mainframe
        mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))

        self.width = tk.IntVar()
        self.width.set(args.width)
        self.dpi = tk.IntVar()
        self.dpi.set(args.dpi)

        if args.cgfile:
            filename,  = args.cgfile
        else:
            filename = askopenfilename(parent=mainframe, title='Please choose a *.cg file.',
                                        filetypes = [("cg files", ("*.cg", "*.coord")), 
                                                     ("All", "*")])
        self.cg = ftmc.CoarseGrainRNA(filename)
        #Image Display
        self.imgDisplay = ttk.Label(mainframe)
        mainframe.columnconfigure(0, minsize = min(int(args.width*3.75+1), 100))
        mainframe.columnconfigure(1, minsize = min(int(args.width*3.75+1), 100))
        mainframe.columnconfigure(2, minsize = min(int(args.width*3.75+1), 100))
        mainframe.columnconfigure(3, minsize = min(int(args.width*3.75+1), 100))
        mainframe.rowconfigure(6, minsize = min(int(args.width*15+2), 400))
        self.imgDisplay.grid(row=6, column = 0, columnspan = 4, sticky = "N")

        #Zoom Selection
        self.zoom = tk.IntVar()
        self.zoom.set(10)
        #self.zoom.set(10)
        ttk.Label(mainframe, text="Zoom:").grid(row = 0, column = 0, sticky = "E")
        tk.Scale(mainframe, from_=1, to_=15, orient=tk.HORIZONTAL, command = self.update, 
                 variable = self.zoom, length=180).grid(row = 0, column = 1, sticky = "W")
        #Nucleotide selection
        self.nucleotidePosition = tk.IntVar()
        self.nucleotidePosition.set(1)
        ttk.Label(mainframe, text="Nucleotide:").grid(row = 1, column = 0, sticky = "E")
        tk.Scale(mainframe, from_=1, to_=len(self.cg.seq), orient=tk.HORIZONTAL, 
                 command = self.update, length=len(self.cg.seq), tickinterval = 100,
                 variable = self.nucleotidePosition).grid(row = 1, column = 1, columnspan=3, sticky = "W")
        #Show virtual residues and stems
        self.showVres = tk.IntVar()
        self.showStems = tk.IntVar()
        ttk.Checkbutton(mainframe,text = "Show virtual residues", command = self.update, 
                        variable = self.showVres).grid(row = 0, column = 2)
        ttk.Checkbutton(mainframe,text = "Show stems", command = self.update, 
                        variable = self.showStems).grid(row = 0, column = 3)
        # Offset
        self.offset_x = tk.IntVar()
        self.offset_y = tk.IntVar()
        ttk.Label(mainframe, text="offset x (Angstrom):").grid(row = 2, column = 0, sticky = "E")        
        tk.Scale(mainframe, from_=-100, to_=100, orient=tk.HORIZONTAL, command = self.update, 
                 variable = self.offset_x, length=180).grid(row = 2, column = 1, sticky = "W")
        ttk.Label(mainframe, text="offset y (Angstrom):").grid(row = 3, column = 0, sticky = "E")        
        tk.Scale(mainframe, from_=-100, to_=100, orient=tk.HORIZONTAL, command = self.update,
                 variable = self.offset_y, length=180).grid(row = 3, column = 1, sticky = "W")
        # Rotation
        self.inplane_rot = tk.IntVar()
        ttk.Label(mainframe, text="rotate (degrees)").grid(row = 2, column = 2, sticky = "E")
        tk.Scale(mainframe, from_=-180, to_=180, orient=tk.HORIZONTAL, command = self.update,
                 variable = self.inplane_rot, length=180).grid(row = 2, column = 3, sticky = "W")
        #Projection Angles
        self.theta = tk.IntVar()
        self.phi = tk.IntVar()
        if self.cg.project_from is not None:
            r,t,p = ftuv.spherical_cartesian_to_polar(self.cg.project_from)
            self.theta.set(int(t))
            self.phi.set(int(p))
        ttk.Label(mainframe, text="Projection angle: THETA").grid(row = 4, column = 0, sticky = "E")
        tk.Scale(mainframe, from_=-180, to_=180, orient=tk.HORIZONTAL, command = self.updateProjection,
                 variable = self.theta, length=180).grid(row = 4, column = 1, sticky = "W")
        ttk.Label(mainframe, text="Projection angle: PHI").grid(row = 4, column = 2, sticky = "E")
        tk.Scale(mainframe, from_=-90, to_=90, orient=tk.HORIZONTAL, command = self.updateProjection,
                 variable = self.phi, length=180).grid(row = 4, column = 3, sticky = "W")
        # Action Buttons
        ttk.Button(mainframe, text="Print Coordinates", 
                   command=self.createLandmark).grid(row = 3, column = 2)
        ttk.Button(mainframe, text="Save Image as...", 
                   command=self.saveImage).grid(row = 3, column = 3)

        # Select width and dpi
        ttk.Label(mainframe, text="Image width in Angstrom").grid(row = 5, column = 0, sticky = "E")
        tk.Scale(mainframe, from_=50, to_=1000, orient=tk.HORIZONTAL, command = self.update,
                 variable = self.width, length=180).grid(row = 5, column = 1, sticky = "W")
        ttk.Label(mainframe, text="Image width in Pixels").grid(row = 5, column = 2, sticky = "E")
        tk.Scale(mainframe, from_=5, to_=80, orient=tk.HORIZONTAL, command = self.update,
                 variable = self.dpi, length=180).grid(row = 5, column = 3, sticky = "W")

        self.proj = None
        self.updateProjection()
    def show(self):      
        ttk.Style().theme_use('clam')      
        self.root.mainloop()

    def updateProjection(self, _=None):
        direction = ftuv.spherical_polar_to_cartesian(np.array([1, math.radians(self.theta.get()), math.radians(self.phi.get())]))
        self.proj=fpp.Projection2D(self.cg, project_virtual_atoms=True, proj_direction=direction,
                          project_virtual_residues = list(range(1, len(self.cg.seq)+1)))
        self.update()
    def update(self, _=None):
        width = self.width.get()
        rest = width % 10
        if rest<5:
            self.width.set(width-rest)
        else:
            self.width.set(width-rest+10)
        proj = self.proj
        #Offset
        box = self.getBox() 
        #Rasterize the image
        if self.showVres.get():
            img, _= proj.rasterize(self.dpi.get(), box, virtual_atoms = True, 
                                   rotate = self.inplane_rot.get(), warn = False)
        else:
            img, _= proj.rasterize(self.dpi.get(), box, virtual_atoms = True, 
                                   rotate = self.inplane_rot.get(), 
                                   virtual_residues = False, warn = False)            
            img = fpp.to_rgb(img)
        pilImg = Image.fromarray(img)
        pilImg = pilImg.convert("RGB")
        #Show stems
        if self.showStems.get():
            stemres=[]
            for s in self.cg.defines.keys():
                if s[0]!="s": continue
                for pos in self.cg.define_residue_num_iterator(s):
                    stemres.append(proj.get_vres_by_position(pos))
            stemres = np.array(stemres)
            rast = fpp.rasterized_2d_coordinates(stemres, self.width.get()/self.dpi.get(), 
                                                 origin = np.array([box[0],box[2]]), 
                                                 rotate = self.inplane_rot.get())                
            for x,y in rast:
                if 0<=x<self.dpi.get() and 0<=y<self.dpi.get():
                    pilImg.putpixel((int(y), int(x)),(100,255,0))
        #Selected Nucleotide
        x,y = self.getSelected()
        if 0<=x<self.dpi.get() and 0<=y<self.dpi.get():
            pilImg.putpixel((int(y), int(x)),(255,0,0))

        #Zoom
        zoom = int(self.zoom.get())     
        newsize = (self.dpi.get()*int(zoom), self.dpi.get()*int(zoom))        
        pilImg = pilImg.resize(newsize, Image.NEAREST)
        #Show
        tkImg = ImageTk.PhotoImage(pilImg)
        self.imgDisplay.configure(image = tkImg)
        self.image = tkImg # Keep a reference!
        self.pilImage = pilImg # For saving
    def getBox(self):
        offset = np.array([self.offset_y.get(), self.offset_x.get()])
        box=fph.get_box(self.proj, self.width.get(), offset )
        return box
    def getSelected(self):
        box = self.getBox()
        selected = int(self.nucleotidePosition.get())
        nucPos = self.proj.get_vres_by_position(selected)
        x, y = fpp.rasterized_2d_coordinates(np.array([[nucPos[0],nucPos[1]]]), 
                                             self.width.get()/self.dpi.get(), 
                                             origin = np.array([box[0],box[2]]), 
                                             rotate = self.inplane_rot.get())[0]
        return x,y
    def createLandmark(self):
        box=self.getBox()
        selected = self.nucleotidePosition.get()
        nucPos = self.proj.get_vres_by_position(selected)
        x, y = self.getSelected()
        print ("{},{},{}".format(selected, x, y))

    def on_closing(self):
        if tkMessageBox.askokcancel("Quit", "This will exit the complete script. Proceed?"):
            self.root.destroy()
            sys.exit(2)

    def saveImage(self):
        filename = asksaveasfilename(parent=self.mainframe, title='Save png image', defaultextension='.png')
        box=self.getBox()
        imgBW, _= self.proj.rasterize(self.dpi.get(), box, virtual_atoms = True, virtual_residues = False,
                                   rotate = self.inplane_rot.get(), warn = True)
        imgColor = self.pilImage
        scipy.misc.imsave(filename, imgBW, "png")
        imgColor.save(filename+"color.png", format="png")

parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    r = RnaDisplay(args)
    r.show()













