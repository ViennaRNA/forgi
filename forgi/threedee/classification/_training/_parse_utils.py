import sys
from collections import namedtuple, defaultdict
"""
Code by me (B.Thiel), borrowed from another, unreleased project.
"""


class LineBasedParser(object):
    comment = "#"
    skip_empty = True
    strip = True
    """
    A baseclass for line-based parsers.

    It provides the member functions `parse` taking a filename,
    and `parse_string`, taking a multiline string.

    If you subclass this class, you have to implement `_parse_line`
    and can implement `_before_parsing` and `_after_parsing`

    The class-level variables `comment`, `skip_empty` and `strip`
    controll how lines without content are handled.

    * If strip is True, all lines are stripped before further processing.
    * If skip_empty is True, empty lines (after stripping) are skipped
      and not passed to _parse_line
    * If comment is non-empty, all lines starting with the comment character(s)
      (after stripping) are skipped as well.
    """

    def parse(self, filename):
        """
        :param filename: A filename. This file will be opened for reading.
        :returns: self.result (has to be filled by self._parse_line)
        """
        with open(filename) as f:
            return self._parse(f)

    def parse_string(self, string):
        """
        :param string: A multiline string
        :returns: see LineBasedParser.parse
        """
        return self._parse(string.splitlines())

    def _before_parsing(self):
        pass

    def _parse_line(self, line):
        pass

    def _after_parsing(self):
        pass


    def _parse(self, lines):
        """
        This function should not be overridden.

        This function handles reading of lines, comments, empty-lines.
        """
        self.result = None
        self._before_parsing()
        for i, line in enumerate(lines):
            if self.strip:
                line = line.strip()
            if self.skip_empty and not line:
                continue
            if self.comment and line.startswith(self.comment):
                continue
            self._parse_line(line)
        self._after_parsing()
        return self.result


class ChainIdMappingParser(LineBasedParser):
    comment = ""  # No comment characters allowed in these files.

    def _before_parsing(self):
        self.result = namedtuple("ChainChainMapping", [
                                 "bundle2mmcif", "mmcif2bundle"])(defaultdict(dict), {})
        self.bundle = None

    def _parse_line(self, line):
        if self.bundle is None:
            if "New chain ID" not in line or "Original chain ID" not in line:
                raise ValueError("Expecting first line to contain "
                                 "'New chain ID' and 'Original chain ID'")
            else:
                self.bundle = "HEADER"
                return
        elif self.bundle == "HEADER" or "-pdb-bundle" in line:
            if not line[4:].startswith("-pdb-bundle"):
                raise ValueError("Expecting '-pdb-bundle'", charno=4)
            else:
                if not line[-1] == ":":
                    raise ValueError(
                        "Expecting pdb-bundle name to be followed by ':'")
                self.bundle = line[:-1]
                return
        else:
            new, old = line.split()
            new = new.strip()
            old = old.strip()
            if len(new) != 1:
                raise ValueError(
                    "Expecting chain in pbd-bundle to be only 1 character.")
            self.result.bundle2mmcif[self.bundle][new] = old
            self.result.mmcif2bundle[old] = (self.bundle, new)
            return

    def _after_parsing(self):
        # Convert the defaultdict to a normal dict.
        self.result = namedtuple("ChainChainMapping",
                                 ["bundle2mmcif", "mmcif2bundle"])(dict(self.result.bundle2mmcif), self.result.mmcif2bundle)
