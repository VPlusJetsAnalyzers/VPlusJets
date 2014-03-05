import atexit
import readline
import rlcompleter
import os

historyPath = os.path.expanduser("~/.pyhistory")

def save_history(historyPath=historyPath):
    import readline
    try:
        readline.set_history_length(1000)
        readline.write_history_file(historyPath)
    except IOError:
        print 'skipping the history writing'

if os.path.exists(historyPath):
    readline.read_history_file(historyPath)

atexit.register(save_history)
del atexit, readline, rlcompleter, save_history, historyPath
