from functools import wraps
from time import time, sleep
from datetime import datetime

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

def retry(exceptions=Exception, tries=3, delay=0):
    def retry_decorator(my_func):
        def retry_wrapper(*args, **kwargs):
            n = tries
            print(f"{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')} Retry of the function: `{my_func.__name__}` started")
            while n != 0:
                try:
                    result = my_func(*args, **kwargs)
                    print(
                        f"{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')} Retry of the function: `{my_func.__name__}` finished.")
                    return result
                except exceptions as e:
                    n -= 1
                    print(f"{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')} {e}")
                    print(f"{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')} Remaining retry: {n}")
                    sleep(delay)
            print(f"{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')} Retry of the function: `{my_func.__name__}` failed.")

        return retry_wrapper

    return retry_decorator