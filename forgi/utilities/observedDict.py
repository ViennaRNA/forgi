from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


class observedDict(dict):
    def __init__(self, value=[], on_change=lambda key: None):
        """
        A dictionary superclass.

        A simplified (and less powerful) version of 
        http://code.activestate.com/recipes/306864-list-and-dictionary-observer/,
        which is written by Bernhard Mulder and licensed under the PSF license.

        Whenever a value changes, is added to or deleted from the dictionary, 
        the function on_change will be called with the key as sole argument.
        """
        super(observedDict, self).__init__(value)
        self.on_change = on_change

    def __setitem__(self, key, value):
        super(observedDict, self).__setitem__(key, value)
        self.on_change(key)

    def __delitem__(self, key):
        super(observedDict, self).__delitem__(key)

    def clear(self):
        for key in self.keys():
            self.on_change(key)
        super(observedDict, self).clear()

    def update(self, other_dict):
        for key in other_dict.keys():
            self.on_change(key)
        super(observedDict, self).update(other_dict)

    def setdefault(self, key, value=None):
        if key not in self:
            self.on_change(key)
        super(observedDict, self).setdefault(key, value)

    def pop(self, k, x=None):
        if k in self:
            self.on_change(k)
        return super(observedDict, self).pop(self, k, x)

    def popitem(self):
        key, value = super(observedDict, self).popitem(self)
        self.on_change(key)
        return key, value
