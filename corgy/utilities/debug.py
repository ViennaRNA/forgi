import inspect    
import sys

def pv(name):
    '''
    Snatched from:

    http://stackoverflow.com/questions/2813227/printing-variable-names-and-contents-as-debugging-tool-looking-for-emacs-python
    '''
    record=inspect.getouterframes(inspect.currentframe())[1]
    frame=record[0]
    val=eval(name,frame.f_globals,frame.f_locals)
    print >>sys.stderr, '{0}: {1}'.format(name, val)
