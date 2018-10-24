from __future__ import print_function
import inspect
import sys
import traceback


def pv(name):
    '''
    Snatched from:

    http://stackoverflow.com/questions/2813227/printing-variable-names-and-contents-as-debugging-tool-looking-for-emacs-python
    '''
    record = inspect.getouterframes(inspect.currentframe())[1]
    frame = record[0]
    val = eval(name, frame.f_globals, frame.f_locals)
    print('{0}: {1}'.format(name, val), file=sys.stderr)


def bt():
    '''
    Print the calls stack.
    '''
    for line in traceback.format_stack()[:-1]:
        print(line.strip())
