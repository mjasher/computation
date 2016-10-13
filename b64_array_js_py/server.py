from flask import Flask, request
import os
import json
import numpy
import base64

app = Flask(__name__, static_path='/static', static_folder=os.path.dirname(__file__))

@app.route('/multiply', methods=['POST'])
def multiply():
    encoded_input = base64.decodestring(request.json['array'])
    decoded_input = numpy.frombuffer(encoded_input, dtype=numpy.float32)
    
    print encoded_input, "is smaller than", json.dumps(decoded_input.tolist())

    twice_the_array = 2. * decoded_input

    encoded_twice_the_array = base64.b64encode(twice_the_array.astype(dtype=numpy.float32))

    return json.dumps({
        "twice_the_array": encoded_twice_the_array
    })

if __name__ == '__main__':
    app.run(debug=True)