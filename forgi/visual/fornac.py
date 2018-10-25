from __future__ import division

fornac_template = """
<!DOCTYPE html>
<meta charset="utf-8">

This is an RNA container.
<div id='rna_ss'> </div>
This after the RNA container.

    <link rel='stylesheet' type='text/css' href='https://raw.githubusercontent.com/pkerpedjiev/fornac/master/css/fornac.css' />
    <script type='text/javascript' src='https://code.jquery.com/jquery-2.1.4.min.js'></script>
    <script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.6/d3.min.js'></script>
    <script type='text/javascript' src='https://rawgit.com/pkerpedjiev/fornac/master/js/fornac.js'></script>

    <script type='text/javascript'>
        var container = new fornac.FornaContainer("#rna_ss",
            {{'applyForce': true, 'allowPanningAndZooming': true, "initialSize": [500,800], 'cssFileLocation': "https://raw.githubusercontent.com/pkerpedjiev/fornac/master/css/fornac.css"}});

        var options = {{'structure': '{}',
                        'sequence': '{}'}};

        colorStr = "{}" 

        container.addRNA(options.structure, options);
        cs = ColorScheme(colorStr);

        container.addCustomColors(cs.colorsJson);
        container.changeColorScheme('custom'); 
    </script>
"""


def create_fornac_page_for_structure(bg, color_string):
    """
    Create a fornac page dispalying this structure. The dotbracket string
    and sequence will automatically be extracted from the BulgeGraph.

    Colors can be specified as a dictionary containing floating 
    point values. These will be uniformly scaled according to the color
    scale passed in.

    :param color_string: The color string to be passed to fornac.
                          e.g. "11-12: red 14-17: rgb(14,15,120)"
    :return: The html text of the resulting web page.
    """
    return fornac_template.format(bg.to_dotbracket_string(),
                                  bg.seq, color_string)


def scale_colors(colors_dict, cmap=None, reverse=False):
    '''
    A dictionary with values containing scalars which need to 
    be scaled according to some color map. The keys are irrelevant.

    The color map will be normalized to the range of values within 
    color_dict.

    :param colors_dict: The dictionary containing the values to be color scaled.
    :param cmap: A color map to be used to scale the colors.
    :param reverse: Reverse the color map
    :return: Another dictionary containing rgb triples as values.
    '''

    if cmap is None:
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('Blues')

    values = colors_dict.values()
    min_value = min(values)
    max_value = max(values)

    new_dict = {}
    for key in colors_dict:
        if reverse:
            color = cmap(
                1 - ((colors_dict[key] - min_value) / (max_value - min_value)))
        else:
            color = cmap(
                (colors_dict[key] - min_value) / (max_value - min_value))
        new_dict[key] = (int(255 * color[0]),
                         int(255 * color[1]), int(255 * color[2]))

    return new_dict


def element_to_nucleotide_colors(bg, element_colors):
    '''
    Convert a dictionary of per-element colors to a dictionary of per-nucleotide colors

    :param element_colors: A dictionary of element colors (e.g. {'i0': (255,0,0), 'm1': {255,255,255)}
    :return: A dictionary of nucleotide numbers: (e.g {1: (14,15,120), 2: (255,255,255)})
    '''
    new_dict = {}
    for key in element_colors:
        for res in bg.define_residue_num_iterator(key):
            new_dict[res] = element_colors[key]

    return new_dict


def nucleotide_colors_to_fornac_color_string(nucleotide_colors):
    '''
    Convert a dictionary of per nucleotide colors to a fornac
    color string.

    :param nucleotide_colors: A dictionary with nucleotide numbers as keys and colors as values.
                              (e.g. {1: (255,0,0), 2: (255,255,0)})
    :return: A color string (e.g "1:rgb(255,0,0) 2:rgb(255,0,0)")
    '''
    color_string = ""
    for key in nucleotide_colors:
        color_string += "{}:rgb({},{},{}) ".format(key, nucleotide_colors[key][0],
                                                   nucleotide_colors[key][1],
                                                   nucleotide_colors[key][2])

    return color_string
