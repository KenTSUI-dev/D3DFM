from functools import wraps
from time import time, sleep

def timing(f):
    @wraps(f)
    def wrapper(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print(f'func:${f.__name__} args:[{args}, {kw}] took: {te-ts:.4f} sec' )
        return result
    wrapper.__name__ = f.__name__
    return wrapper