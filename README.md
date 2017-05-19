#### forgi Library for RNA Secondary Structure Analysis ####

Full documentation: [https://pkerpedjiev.github.io/forgi/]

##### Build the documentation #####

sphinx-build -b html doc ~/public_html/forgi

##### Installation #####

python setup.py install

##### Examples #####

* Color elements by name in a gradient and visualize a coarse-grain RNA

  >>> python examples/visualize_cg.py examples/1y26.cg --color-gradual h0,s2,m1,s0,m0,m2,s3,h1 -x --virtual-atoms --sidechain-atoms

* Create HTML code for a fornaContainer with colors according to the order of coarse grain elements provided on the commandline

  >>> python examples/ordered_elements_to_forna_colors.py examples/1y26.cg h0,s2,m1,s0,m0,m2,s3,h1

* Display an rna structure, save it to file, trim it and view it using eog

  >>> pdb_id="ideal_1_3_4_6.pdb"; venv/bin/python examples/visualize_pdb.py forgi/threedee/data/${pdb_id} --output /tmp/${pdb_id}.png; mogrify -trim /tmp/${pdb_id}.png; eog /tmp/${pdb_id}.png

