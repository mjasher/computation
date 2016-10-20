import base64
import numpy


def js2np(json_encoded):
    byte_string = base64.decodestring(json_encoded['__array__'])
    np_array = numpy.frombuffer(byte_string, dtype=json_encoded['dtype']).reshape(json_encoded['shape'])
    return np_array


def np2js(np_array):
	return {
		"shape": np_array.shape,
		"dtype": str(np_array.dtype),
		"__array__": base64.b64encode(np_array.tobytes())
	}
