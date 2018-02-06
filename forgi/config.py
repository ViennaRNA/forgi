import json
import os
import os.path
import logging
import appdirs

log = logging.getLogger(__name__)

dirs = appdirs.AppDirs("forgi", "TBI")

def read_config():
    config = {}
    for directory in [dirs.site_config_dir, dirs.user_config_dir]:
        filename = os.path.join(directory, "config.json")
        try:
            with open(filename) as f:
                conf = json.load(f)
        except (OSError, IOError):
            log.debug("No configuration file present at %s", filename)
            pass
        else:
            log.debug("Reading configuration from %s", filename)
            config.update(conf)

    return config

def set_config(key, value):
    filename = os.path.join(dirs.user_config_dir, "config.json")
    try:
        with open(filename) as f:
            config = json.load(f)
    except (OSError, IOError):
        config = {}
    config[key]=value
    try:
        os.makedirs(dirs.user_config_dir)
    except OSError:
        pass
    with open(filename, "w") as f:
        json.dump(config, f)
        log.info("Configuration file %s updated", filename)
