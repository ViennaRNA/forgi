from __future__ import absolute_import
from __future__ import print_function

from builtins import zip
from builtins import next
from builtins import range

import unittest
import itertools as it
from pprint import pprint
import collections as col
import sys
import os
import logging

from forgi.graph._domain_dict import _DefineDictionary, _DomainSupportingDictionary, _EdgeDictionary
import forgi.graph.bulge_graph as fgb

log=logging.getLogger(__name__)

class BaseClassTest(unittest.TestCase):
    def test_dictionary_interface(self):
        """
        Behaves like a dict.
        """
        a=_DomainSupportingDictionary({"s0":"bla", "s1":"blo"}, None)
        self.assertEqual(a["s0"], "bla")
        log.error("Done")
        log.error("Now assert")
        assert "s0" in a
        log.error("Now check IN")
        self.assertIn("s0", a)
        a.update({"s3":"bli", "s1":"blo"})
        self.assertIn("s3", a)


class TestDefineDIctionary(unittest.TestCase):
    db = "(.(.((...))..))"
    bg = fgb.BulgeGraph.from_dotbracket(db)
