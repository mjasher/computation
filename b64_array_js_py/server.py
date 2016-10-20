from flask import Flask, request
import os
import json
import numpy
import np2js

app = Flask(__name__, static_path='/static', static_folder=os.path.dirname(__file__))

@app.route('/multiply', methods=['POST'])
def multiply():
    decoded_input = np2js.js2np(request.json)
    
    print request.json['__array__'], "is smaller than", json.dumps(decoded_input.tolist())

    twice_the_array = 2. * decoded_input

    encoded_twice_the_array = np2js.np2js(twice_the_array.astype(dtype=numpy.float32))

    return json.dumps({
        "twice_the_array": encoded_twice_the_array
    })

if __name__ == '__main__':
    app.run(debug=True)