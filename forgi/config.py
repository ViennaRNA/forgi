"""
Code for reading and writing the forgi configuration files.

To configure forgi, run forgi_config.py
"""


import json
import os
import os.path
import logging
import appdirs

log = logging.getLogger(__name__)

dirs = appdirs.AppDirs("forgi", "TBI")  # pylint: disable=C0103
ALLOWED_KEY_VALUES = {"PDB_ANNOTATION_TOOL": ["MC-Annotate", "DSSR", "forgi"]}


def iter_configfiles():
    """
    Iterate over config file locations, from lowest priority to highest priority.
    The file does not have to exist.
    """
    for directory in [dirs.site_config_dir, dirs.user_config_dir]:
        filename = os.path.join(directory, "config.json")
        yield filename


def read_config():
    """
    Read all configuration files and return a configuration dictionary

    :returns: A dictionary.
    """
    config = {}
    for filename in iter_configfiles():
        try:
            with open(filename) as f:
                conf = json.load(f)
        except (OSError, IOError):
            log.debug("No configuration file present at %s", filename)
        else:
            log.debug("Reading configuration from %s", filename)
            config.update(conf)
    return config
