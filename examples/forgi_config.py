#!/usr/bin/python
from builtins import input

import argparse
import json
import sys


import forgi.config

parser = argparse.ArgumentParser(
    "Create/ update configuration files for forgi.")

if __name__ == "__main__":
    config = {}
    for filename in forgi.config.iter_configfiles():
        try:
            with open(filename) as f:
                conf = json.load(f)
        except (OSError, IOError):
            print("No config file exists yet at {}".format(filename))
        else:
            if conf:
                print("The following values are set in {}:".format(filename))
                for k, v in conf.items():
                    print("   {:15s}: {}".format(k, repr(v)))
                config.update(conf)
            else:
                print("This config file is empty")
    for filename in forgi.config.iter_configfiles():
        update = input("Update the file {}? (y/N)".format(filename))
        if update in ["y", "Y"]:
            with open(filename) as f:
                conf = json.load(f)
            for key in forgi.config.ALLOWED_KEY_VALUES:
                if key in config:
                    print("After parsing all config files, {} is set to {}".format(
                        key, config[key]))
                    if key in conf:
                        print("In the file you are editing, {} is set to {}".format(
                            key, conf[key]))
                    else:
                        print(
                            "In the file you are editing, {} is not set".format(key))
                else:
                    print("Currently, {} is not set in any config file".format(key))
                while True:
                    new_value = input("Choose value for {} ({}, empty input leave it unchanged, 'X' to delete the entry from the config file)".format(
                        key, "/".join(forgi.config.ALLOWED_KEY_VALUES[key])))
                    if not new_value:
                        break
                    if new_value in forgi.config.ALLOWED_KEY_VALUES[key]:
                        conf[key] = new_value
                        break
                    if new_value == "X":
                        del conf[key]
                        break
                    # Else next loop execution
            with open(filename, "w") as f:
                json.dump(conf, f)
                print("Configuration file {} updated".format(filename))
