import pandas as pd
from docopt import docopt
from diffsel_script_utils import param, data
import diffsel_script_utils as utils
from subprocess import Popen, PIPE, call
from os.path import isfile
from os import remove

def strip(s):
    if s[0] == '-' or s[0] == '<':
        return strip(s[1:])
    elif s[-1:] == '>':
        return s[:-1]
    else:
        return s

def to_str(thing):
    if isinstance(thing, int):
        return str(thing)
    elif not isinstance(thing, str):
        return str(thing[0])
    else:
        return str(thing)

def STEP(s, mute = False):
    if not mute:
        print(utils.step(s))

def INFO(s, mute = False):
    if not mute:
        print("-- "+s)

def ERROR(s, mute = False):
    if not mute:
        print("\n\t["+utils.boldred("Error")+"] "+s)

def BAD(s, mute = False):
    if not mute:
        print(utils.bad(s))

def GOOD(s, mute = False):
    if not mute:
        print(utils.good(s))

def getargs(s):
    args = docopt(s)
    for arg in args:
        if arg != "--help":
            INFO(strip(arg)+" set to "+param(to_str(args[arg])))
    return args
