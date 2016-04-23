import os
import simplejson as json  # import json
import cPickle as pickle  # import pickle
from functools import wraps
import numpy
import datetime
import pandas

"""


Create a cache for a function using

```
@money(key_func, cache_directory)
def f(x):
    return y
```

Similar to https://docs.python.org/3/library/functools.html#functools.lru_cache

Use numpy_key where x is a numpy array.

Using json instead of pickle appears to be faster, but requires custom encoder and decoder.

michael.james.asher@gmail.com
"""

"""
Custom json encoder allows sevearal data types to be saved as json
"""
class DateNumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime.datetime):
            return obj.date().isoformat()
        elif isinstance(obj, datetime.date):
            return obj.isoformat()
        elif isinstance(obj, datetime.timedelta):
            return (datetime.datetime.min + obj).date().isoformat()
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        elif isinstance(obj, pandas.DataFrame):
            try:
                return obj.to_json()
            except:
                return obj.as_matrix().tolist()
        else:
            return super(DateNumpyEncoder, self).default(obj)

"""
Some useful key functions
"""

def numpy_key(original_x, *args, **kwargs):
    x = original_x.copy()
    x.flags.writeable = False
    key = str(hash(x.data))
    return key

def join_key(*args, **kwargs):
    return str(hash(''.join([json.dumps(a, cls=DateNumpyEncoder, ignore_nan=True) for a in args + tuple(kwargs.values())])))

def concat_key(*args, **kwargs):
    return '_'.join(str(a) for a in args + tuple(kwargs.values()))

def date_place_key(geom, years):
    return str(hash(''.join([str(y) for y in years]) + json.dumps(geom)))

def one_key(*args, **kwargs):
    return "one_key"

"""
A decorator to cache a some_function(args) by saving each run in <cache_dir>/some_function/<key_func(args)>

"""


def money(key_func=lambda x: str(hash(json.dumps(x))),
                   cache_directory=None):
    
    def decorator(func):
        if cache_directory is None:
            local_cache_directory = os.path.join(os.environ.get('CACHE_DIR', '/tmp'), 'farmlogs-light-' +func.__name__)
        else:
            local_cache_directory = cache_directory



        @wraps(func)
        def f_with_cache(*args, **kwargs):
            # if 'cache_force_refresh' in kwargs:
            #     force_refresh = kwargs['cache_force_refresh']
            #     del kwargs['cache_force_refresh']
            # else:
            #     force_refresh = False

            key = key_func(*args, **kwargs)
            
            if not os.path.exists(local_cache_directory):
                os.mkdir(local_cache_directory)

            cache_file = os.path.join(local_cache_directory, key)

            # if os.path.exists(cache_file) and not force_refresh:
            if os.path.exists(cache_file):
                with open(cache_file, "rb") as file:
                    y = pickle.load(file)
                    # y = json.load(file)
                return y
            else:
                y = func(*args, **kwargs) 
                with open(cache_file, "wb") as file:
                    # json.dump(y, file, cls=DateNumpyEncoder, ignore_nan=True)
                    pickle.dump(y, file, protocol=pickle.HIGHEST_PROTOCOL)
                return y

        return f_with_cache

    return decorator

"""
demo
"""

if __name__ == '__main__':
    
    import time
    import shutil
    
    def slow_function(x, a, b):
        time.sleep(0.003)
        return 2*x

    rng = numpy.random.RandomState( 0 )

    xs = rng.normal(loc=0.0, scale=1.0, size=(1000,1))

    t0 = time.time()
    ys = [slow_function(x, 5, 10) for x in xs] 
    print time.time()-t0, " sec. for ", len(xs), " samples *without* cache" 

    cache_directory = os.path.join(os.environ.get('CACHE_DIR', '/tmp'), 'test_func') 
    @money(key_func=numpy_key, cache_directory=cache_directory)
    def cached_slow_function(x, a, b):
        time.sleep(0.003)
        return 2 * x

    t0 = time.time()
    ys_cached = [cached_slow_function(x, 5, 10) for x in xs] 
    print time.time()-t0, " sec. for ", len(xs), " samples with cache" 

    assert numpy.allclose(ys, ys_cached)

    t0 = time.time()
    ys = [cached_slow_function(x, 5, 10) for x in xs] 
    print time.time()-t0, " sec. for ", len(xs), " samples with decorator" 

    # clean up demo
    shutil.rmtree(cache_directory)