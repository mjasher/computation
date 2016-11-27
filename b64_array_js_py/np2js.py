import base64
import numpy
import datetime

from flask.json import JSONEncoder, JSONDecoder


def js2np(dct):
    a = numpy.frombuffer(base64.b64decode(dct['data']), dtype=dct['data_type']).reshape(dct['data_shape'])
    a.setflags(write=1)
    return a


def np2js(obj):
	return {'type': 'numpy.ndarray', 'data_type': str(obj.dtype), 'data_shape': obj.shape, 'data': base64.b64encode(obj.tobytes())}


class DateNumpyEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return np2js(obj)
        elif isinstance(obj, datetime.datetime):
            return {'type': 'datetime.datetime', 'data': obj.strftime("%Y-%m-%dT%H:%M:%S.%fZ")}
        elif isinstance(obj, datetime.date):
            return {'type': 'datetime.date', 'data': obj.strftime("%Y-%m-%d")}
        else:
            return super(DateNumpyEncoder, self).default(obj)


class DateNumpyDecoder(JSONDecoder):
    def __init__(self, *args, **kwargs):
        super(DateNumpyDecoder, self).__init__(*args,
            object_hook=self.custom_obj_hook, **kwargs)

    def custom_obj_hook(self, dct):
        if isinstance(dct, dict) and 'type' in dct and dct['type'] == 'numpy.ndarray':
        	return js2np(dct)
        elif isinstance(dct, dict) and 'type' in dct and dct['type'] == 'datetime.datetime':
            return datetime.datetime.strptime(dct['data'], "%Y-%m-%dT%H:%M:%S.%fZ")
        elif isinstance(dct, dict) and 'type' in dct and dct['type'] == 'datetime.date':
            return datetime.datetime.strptime(dct['data'], "%Y-%m-%d").date()
        else:
            return dct  
