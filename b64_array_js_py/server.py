from flask import Flask, request, jsonify
import os
import json
import numpy
import np2js

app = Flask(__name__, static_path='/static', static_folder=os.path.dirname(__file__))

app.json_decoder = np2js.DateNumpyDecoder
app.json_encoder = np2js.DateNumpyEncoder


@app.route('/multiply', methods=['POST'])
def multiply():
    decoded_input = request.json
    
    print np2js.np2js(decoded_input)['data'], "is smaller than", json.dumps(decoded_input.tolist())

    twice_the_array = 2. * decoded_input

    return jsonify({
        "twice_the_array": twice_the_array.astype(dtype=numpy.float32)
    })

if __name__ == '__main__':
    app.run(debug=True)