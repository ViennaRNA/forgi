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
import sys, argparse
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


parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    #ttk.Style().theme_use('clam')

    def on_closing():
        if tkMessageBox.askokcancel("Quit", "This will exit the complete script. Proceed?"):
            root.destroy()
            sys.exit(2)
    root = tk.Tk()
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.title("Manual Projection Rasterization")
    mainframe = ttk.Frame(root, padding="3 3 12 12")
    mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))

    if args.cgfile:
        filename,  = args.cgfile
    else:
        filename = askopenfilename(parent=mainframe, title='Please choose a *.cg file.',
                                    filetypes = [("cg files", ("*.cg", "*.coord")), 
                                                 ("All", "*")])
    cg=ftmc.CoarseGrainRNA(filename)
    try:
        proj=fpp.Projection2D(cg, project_virtual_atoms=True, 
                              project_virtual_residues = list(range(1,len(cg.seq)+1)))
    except ValueError:
        proj=fpp.Projection2D(cg, project_virtual_atoms=True, proj_direction=(0,1,1), #ftuv.get_random_vector(), 
                              project_virtual_residues = list(range(1, len(cg.seq)+1)))


    imgDisplay = ttk.Label(mainframe)
    mainframe.columnconfigure(0, minsize = min(int(args.width*3.75+1), 100))
    mainframe.columnconfigure(1, minsize = min(int(args.width*3.75+1), 100))
    mainframe.columnconfigure(2, minsize = min(int(args.width*3.75+1), 100))
    mainframe.columnconfigure(3, minsize = min(int(args.width*3.75+1), 100))
    mainframe.rowconfigure(5, minsize = min(int(args.width*15+2), 400))
    imgDisplay.grid(row=5, column = 0, columnspan = 4, sticky = "N")
    w = ttk.Label(mainframe, text="Zoom:")
    w.grid(row = 0, column = 0, sticky = "E")        
    def show(z=None):
        showImage()
    sc = tk.Scale(mainframe, from_=1, to_=15, orient=tk.HORIZONTAL, command = show)
    sc.grid(row = 0, column = 1, sticky = "W")
    ttk.Label(mainframe, text="Nucleotide:").grid(row = 1, column = 0, sticky = "E")

    showVresVar = tk.IntVar()
    showStemsVar = tk.IntVar()    
    ttk.Checkbutton(mainframe,text = "Show virtual residues", command = show, variable = showVresVar).grid(row = 0, column = 2)
    ttk.Checkbutton(mainframe,text = "Show stems", command = show, variable = showStemsVar).grid(row = 1, column = 2)

    w = ttk.Label(mainframe, text="offset x (Angstrom):").grid(row = 2, column = 0, sticky = "E")        
    xOffset = tk.Scale(mainframe, from_=-100, to_=100, orient=tk.HORIZONTAL, command = show)
    xOffset.grid(row = 2, column = 1, sticky = "W")

    w = ttk.Label(mainframe, text="offset y (Angstrom):").grid(row = 3, column = 0, sticky = "E")        
    yOffset = tk.Scale(mainframe, from_=-100, to_=100, orient=tk.HORIZONTAL, command = show)
    yOffset.grid(row = 3, column = 1, sticky = "W")

    w = ttk.Label(mainframe, text="rotate (degrees)").grid(row = 2, column = 2, sticky = "E")        
    inplaneRotate = tk.Scale(mainframe, from_=-180, to_=180, orient=tk.HORIZONTAL, command = show)
    inplaneRotate.grid(row = 2, column = 3, sticky = "W")

    nucSelection = tk.Scale(mainframe, from_=1, to_=len(cg.seq), orient=tk.HORIZONTAL, command = show)
    nucSelection.grid(row = 1, column = 1, sticky = "W")
    #selectedNucInfo = ttk.Label(mainframe)
    #selectedNucInfo.grid(row = 2, column = 0)      

    def showImage():
        #Offset
        offset = np.array([yOffset.get(), xOffset.get()])
        box=fph.get_box(proj, args.width, offset )
        #Image with or without virtual residues
        if showVresVar.get():
            img, _= proj.rasterize(args.dpi, box, virtual_atoms = True, 
                                   rotate = inplaneRotate.get(), warn = False)
        else:
            img, _= proj.rasterize(args.dpi, box, virtual_atoms = True, 
                                   rotate = inplaneRotate.get(), 
                                   virtual_residues = False, warn = False)
            img = fpp.to_rgb(img)
        newImg = Image.fromarray(img)
        #Show stems:
        if showStemsVar.get():
            stemres=[]
            for s in cg.defines.keys():
                if s[0]!="s": continue
                for pos in cg.define_residue_num_iterator(s):
                    stemres.append(proj.get_vres_by_position(pos))
            stemres = np.array(stemres)
            rast = fpp.rasterized_2d_coordinates(stemres, args.width/args.dpi, 
                                                 origin = np.array([box[0],box[2]]), 
                                                 rotate = inplaneRotate.get())                
            newImg = newImg.convert("RGB")        
            for x,y in rast:
                if 0<=x<args.dpi and 0<=y<=args.dpi:
                    newImg.putpixel((int(y), int(x)),(100,255,0))
        #Selected Nucleotide
        selected = int(nucSelection.get())
        nucPos = proj.get_vres_by_position(selected)
        x, y = fpp.rasterized_2d_coordinates(np.array([[nucPos[0],nucPos[1]]]), 
                                             args.width/args.dpi, 
                                             origin = np.array([box[0],box[2]]), 
                                             rotate = inplaneRotate.get())[0]
        newImg = newImg.convert("RGB")
        if 0<=x<args.dpi and 0<=y<args.dpi:
            newImg.putpixel((int(y), int(x)),(255,0,0))
        #Zoom
        zoom = int(sc.get())     
        newsize = (args.dpi*int(zoom), args.dpi*int(zoom))        
        newImg = newImg.resize(newsize, Image.NEAREST)
        #Show
        showimg = ImageTk.PhotoImage(newImg)
        imgDisplay.configure(image = showimg)
        imgDisplay.image = showimg # Keep a reference!
        imgDisplay.pilImage = newImg # For saving
    def createLandmark():
        offset = np.array([yOffset.get(), xOffset.get()])
        box=fph.get_box(proj, args.width, offset )
        selected = int(nucSelection.get())
        nucPos = proj.get_vres_by_position(selected)
        x, y = fpp.rasterized_2d_coordinates(np.array([[nucPos[0],nucPos[1]]]), 
                                             args.width/args.dpi, 
                                             origin = np.array([box[0],box[2]]), 
                                             rotate = inplaneRotate.get())[0]
        print ("{},{},{}".format(selected, x, y))

    def saveImage():
        filename = asksaveasfilename(parent=mainframe, title='Save png image', defaultextension='.png')
        offset = np.array([yOffset.get(), xOffset.get()])
        box=fph.get_box(proj, args.width, offset )
        imgBW, _= proj.rasterize(args.dpi, box, virtual_atoms = True, virtual_residues = False,
                                   rotate = inplaneRotate.get(), warn = False)
        imgColor = imgDisplay.pilImage
        scipy.misc.imsave(filename, imgBW, "png")
        imgColor.save(filename+"color.png", format="png")
    ttk.Button(mainframe, text="Print Coordinates", command=createLandmark).grid(row = 3, column = 2)
    ttk.Button(mainframe, text="Save Image as...", command=saveImage).grid(row = 3, column = 3)
    ttk.Style().theme_use('clam')
    root.mainloop()

