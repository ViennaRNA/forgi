#!/usr/bin/python

import string
import sys
from sys import argv,stderr,exit,stdout
from os import popen,system
from os.path import exists,dirname,basename,abspath
from optparse import OptionParser
import os

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u'
              }

def convert_pdb_to_rosetta_pdb(lines, chainids, ignore_chain, removechain):
    """
    Convert a regular pdb file to a rosetta ready pdb file. The contents of the original file
    are passed in as an array of lines and the contents of the fasta file and the new
    pdb file are returned as a pair of strings.

    @param chainids Some parameter
    @param ignore_chain Some other parameter
    """
    out_text = ''
    out_fasta = ''

    oldresnum = '   '
    count = 0;

    num_model = 0
    max_model = 0 # for virus

    for i in range( len( chainids ) ) :
        if chainids[i] == '_':
            chainids[i] = ' '

    goodnames = [' rA',' rC',' rG',' rU']
    for line in lines:
        if len(line)>5 and line[:6]=='ENDMDL':
            num_model += 1
            if num_model > max_model:  break #Its an NMR model.
        if len(line) <= 21:  continue
        if (line[21] in chainids or ignore_chain):
            line_edit = line

            if line[0:3] == 'TER':
                continue
            elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit)>75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5BU'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='OMC'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5MC'):
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]==' DC'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='CBR'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='CB2'):
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='2MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='H2U'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='PSU'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5MU'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='OMG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='7MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='1MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]==' YG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='1MA'):
                line_edit = 'ATOM  '+line[6:17]+'  A'+line[20:]

            #Don't save alternative conformations.
            if line[16] == 'A':
                continue;

            if line_edit[0:4] == 'ATOM':
                resnum = line_edit[23:26]
                if line_edit[26] != ' ':
                    # Ignore alternate atoms
                    continue

                if not resnum == oldresnum: #  or line_edit[12:16] == ' P  ':
                    longname = line_edit[17:20]
                    if longname == '  G':
                        longname = ' rG'
                    elif longname == '  A':
                        longname =   ' rA'
                    elif longname == '  C':
                        longname =   ' rC'
                    elif longname == '  U':
                        longname =   ' rU'
                    elif longname == 'G  ':
                        longname =   ' rG'
                    elif longname == 'A  ':
                        longname =   ' rA'
                    elif longname == 'C  ':
                        longname =   ' rC'
                    elif longname == 'U  ':
                        longname =   ' rU'
                    #elif longname == ' DG':
                    #    longname = ' rG'
                    #elif longname == ' DA':
                    #    longname = ' rA'
                    #elif longname == ' DC':
                    #    longname = ' rC'
                    #elif longname == ' DT':
                    #    longname = ' rU'
                    elif longname == 'GUA':
                        longname = ' rG'
                    elif longname == 'ADE':
                        longname = ' rA'
                    elif longname == 'CYT':
                        longname = ' rC'
                    elif longname == 'URA':
                        longname = ' rU'
                    elif longname == 'URI':
                        longname = ' rU'
                    else:
                        if longname not in goodnames:    continue

                    if longer_names.has_key(longname):
                        out_fasta += longer_names[longname]
                    else:
                        out_fasta += 'X'

                    #print "AAH ==> " ,  resnum, oldresnum, line_edit
                    count = count + 1

                oldresnum = resnum

                if not longname in goodnames:
                    continue

                newnum = '%4d' % count
                line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + \
                            newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+'  '+line_edit[23:]

                line_edit = line_edit.replace( 'HO2\'', '2HO*' )
                line_edit = line_edit.replace( 'HO5\'', '5HO*' )
                line_edit = line_edit.replace( 'H5\'\'', '2H5*' )

                line_edit = line_edit.replace('\'','*')
                line_edit = line_edit.replace('OP1','O1P')
                line_edit = line_edit.replace('OP2','O2P')

                out_text += line_edit

    return (out_fasta, out_text)

def get_option_parser():
    parser = OptionParser()
    parser.add_option("-o", "--out-pdb", dest="out_pdb", default=None)
    parser.add_option("-n", "--nochain", dest="nochain", action="store_true", default=False)
    parser.add_option("-i", "--ignorechain", dest="ignorechain", action="store_true", default=False)
    parser.add_option("-f", "--out-fasta", dest="out_fasta", default=None)
    parser.add_option("-b", "--batch", dest="batch", action="store_true", default=False)
    parser.add_option("-d", "--output-dir", dest="out_dir", default='.')

    return parser

def main():
    if len(argv) < 2:
        print "Usage: ./make_rna_rosetta_ready.py pdb_file1 [chain_id1] [chain_id2]"
        print 
        print "Convert a file to the rosetta pdb format"

        exit(1)

    (options, args) = get_option_parser().parse_args()
    #print "options:", options
    #print "args:", args

    if options.out_pdb and options.batch:
        print >>sys.stderr, "Incompatible options --out-pdb and --batch"
        exit(1)

    chainids = []

    if not options.batch and len( args ) > 1:
        chainids = argv[1:]
    else:
        options.ignorechain = True

    if options.batch:
        numstructs = len(args)
    else: 
        numstructs = 1

    for i in range(numstructs):
        pdbname = args[i]

        if options.batch:
            only_pdb = os.path.basename(os.path.splitext(pdbname)[0])
            
            out_filename = only_pdb + ".pdb"
            outid = open(os.path.join(options.out_dir, out_filename), 'w')
        else:
            if options.out_pdb == None:
                outid = stdout
            else:
                outid = open(os.path.join(options.out_dir, options.out_pdb), 'w')

        if options.out_fasta == None:
            fastaid = stderr
        else:
            outid = open(options.out_fasta, 'w')

        print >>sys.stderr, "pdbname", pdbname
        print >>sys.stderr, exists(pdbname)

        if not exists( pdbname ) and pdbname != '-':
            print >>stderr, 'DOES NOT EXIST: ', pdbname
            exit(1)

        if pdbname[-3:] == '.gz':
            lines = popen(' gzcat ' + pdbname).readlines()
        elif pdbname == '-':
            lines = sys.stdin.readlines()
        else:
            lines = open(pdbname,'r').readlines()

        (out_fasta, out_text) = convert_pdb_to_rosetta_pdb(lines, chainids, options.ignorechain, options.nochain)
        outid.write(out_text)

        #print 'Writing ... '+fastafile
        fastaid.write('>'+pdbname+'\n');
        fastaid.write(out_fasta);
        fastaid.write('\n')
        if outid != stdout:
            outid.close()
        if fastaid != stderr:
            fastaid.close()

if __name__ == '__main__':
    main()
