/*

Taken from
  - https://github.com/hughsk/tab64
  - https://developer.mozilla.org/en-US/docs/Web/JavaScript/Base64_encoding_and_decoding#Appendix.3A_Decode_a_Base64_string_to_Uint8Array_or_ArrayBuffer

*/


// output can be read with eg. numpy.frombuffer(base64.decodestring(request.json['array']), dtype=numpy.float32)
function encode(input) {
     var buffer = new Uint8Array(input.buffer);
     var binary = '';
     var len = buffer.byteLength;
     for (var i = 0; i < len; i++) {
         binary += String.fromCharCode(buffer[i]);
     }
     return window.btoa(binary);
}

function UserException(message) {
   this.message = message;
   this.name = "UserException";
}

var dtype = function(dtype_str) {
  switch (dtype_str) {
    case 'int8':
      return Int8Array
    case 'int16':
      return Int16Array
    case 'int32':
      return Int32Array
    case 'uint8':
      return Uint8Array
    case 'uint16':
      return Uint16Array
    case 'uint32':
      return Uint32Array
    case 'float32':
      return Float32Array
    case 'float64':
      return Float64Array
    case 'array':
      return Array
    case 'uint8_clamped':
      return Uint8ClampedArray
    default:
      throw new UserException("invalid dtype")
  }
}

var ceil = Math.ceil

function b64int(n) {
  return n < 26 ? n + 65
    : n < 52 ? n + 71
    : n < 62 ? n - 4
    : n === 62 ? 43
    : n === 63 ? 47
    : 65
}

function intb64(chr) {
  return chr > 64 && chr < 91 ? chr - 65
    : chr > 96 && chr < 123 ? chr - 71
    : chr > 47 && chr < 58 ? chr + 4
    : chr === 43 ? 62
    : chr === 47 ? 63
    : 0
}

// input can be made eg. base64.b64encode(some_aray.astype(dtype=numpy.float32))
function decode(input, output) {
  input = input.replace(/[^A-Za-z0-9\+\/]/g, "")

  var inputLength = input.length
  var outputLength = inputLength * 3 + 1 >> 2
  var outidx = 0
  var inidx = 0
  var rvalue

  if (!output) output = new Uint8Array(outputLength)
  if (typeof output === 'string') {
    var type = output
    var bytes = parseInt(type.match(/[0-9]+/g), 10) / 8
    var offset = ceil(outputLength / bytes) * bytes - outputLength
    if (bytes) outputLength += offset
    output = new Uint8Array(outputLength)
    rvalue = new (dtype(type))(output.buffer)
  } else {
    rvalue = output
  }

  for (var value = 0; inidx < inputLength; inidx++) {
    var bit = inidx & 3

    value |= intb64(
      input.charCodeAt(inidx)
    ) << (18 - 6 * bit)

    if (bit === 3 || inputLength - inidx === 1) {
      for (var sbit = 0; sbit < 3 && outidx < outputLength; sbit++) {
        output[outidx++] = value >>> (16 >>> sbit & 24) & 255
      }
      value = 0
    }
  }
  return rvalue
}