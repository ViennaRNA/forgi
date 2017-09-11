from builtins import object
import os
import forgi

class Configuration(object):
    mids_method="template"
    #mids_method="basenormals"
    base_dir = os.path.expanduser('.')
    #data_base_dir = os.path.expanduser('~/data/ernwin/processed')
    #pdb_base_dir = os.path.expanduser('~/data/ernwin/pdb')
    stem_fragment_dir = os.path.join(base_dir, 'forgi/data')
    #lric_stats_fn = os.path.join(base_dir, 'fess/stats/temp.energy')
    #template_residue_fn = os.path.join(base_dir, 'fess/stats/residue_template.pdb')
    #longrange_contact_stats_fn = os.path.join(base_dir, 'fess/stats/temp.longrange.contact')

    #test_input_dir = os.path.expanduser('~/data/ernwin/processed/')
    #test_output_dir = os.path.join(base_dir, "test_output")
    #sampling_output_dir = 'best'
    #barnacle_dir = '/scr/plastilin/pkerp/apps/Barnacle'
    #stem_library = dict()
