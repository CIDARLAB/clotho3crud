import importlib
import json
from bunch import *
import sys


class ClothoError(Exception):
    pass

class ClothoMethod(object):
    def __init__(self, name, context):
        self._name = name
        self._context = context

    def __call__(self, *args):
        self._context.send_value(
            {"type": "api", "name": self._name, "args": args})
        reply_value = read_value()
        reply_type = reply_value["type"]
        if reply_type == "api":
            return reply_value["return"]
        elif reply_type == "api_error":
            raise ClothoError(reply_value["message"])
        raise RuntimeError("bad API reply")

class Clotho(object):
    def __init__(self, context):
        self._context = context

    def __getattr__(self, name):
        return ClothoMethod(name, self._context)

class Context(object):
    '''Execution context object

    Setup by initialization object
    '''
    def __init__(self, init_obj, swapper):
        if init_obj["type"] != "func":
            raise ValueError
        self.code = init_obj["code"]
        self.args = init_obj["args"]
        self._swapper = swapper

    @staticmethod
    def get_serializable(obj):
        if "_json" in dir(obj):
            return obj._json()
        elif hasattr(obj, "__dict__"):
            return obj.__dict__
        else:
            try:
                return list(iter(obj))
            except TypeError:
                pass
        raise TypeError(repr(obj) + " is not JSON serializable")

    def send_value(self, value):
        '''Send one message to host'''
        with self._swapper:
            sys.stdout.write(
                    json.dumps(value, check_circular=True, 
                        default=Context.get_serializable)
                    .encode("UTF-8"))
            sys.stdout.write(b'\0')
            sys.stdout.flush()

class StdoutSwapper(object):
    def __init__(self, orig, dest):
        self._orig = orig
        self._dest = dest

    def __enter__(self):
        sys.stdout = self._dest

    def __exit__(self, exc_type, exc_value, traceback):
        # need to restore stdout regardless of whether an exception occured
        sys.stdout = self._orig

    
def read_value():
    '''Read one message from host through standard input and return it'''
    buf = bytearray()
    while 1:
        c = sys.stdin.read(1)
        if c == "":
            raise ValueError("received incomplete JSON value")
        if c == "\0":
            break
        buf.append(c)
    return json.loads(bytes(buf).decode("UTF-8"), 
            object_hook=bunchify)

def main():
    # redirect prints to standard error
    orig_stdout = sys.stdout
    sys.stdout = sys.stderr

    # fetch initialization object
    context = Context(read_value(), StdoutSwapper(sys.stderr, orig_stdout))

    # run user code; get return value
    user_return_value = body(context)

    # send return value back to host
    context.send_value({"type": "func", "return": user_return_value})

def body(context):
    # execute user code body with given variable bindings
    scope_dict = {"clotho": Clotho(context), "ClothoError": ClothoError}
    exec(context.code, scope_dict)

    # execute user function, which is always called "run"
    return scope_dict["run"](*context.args)

main()
