_nipy_spectral_data = {
    'red': [[0.0, 0.1, 0.1], [0.05, 0.4667, 0.4667],
    [0.10, 0.5333, 0.5333], [0.15, 0.0, 0.0],
    [0.20, 0.0, 0.0], [0.25, 0.0, 0.0],
    [0.30, 0.0, 0.0], [0.35, 0.0, 0.0],
    [0.40, 0.0, 0.0], [0.45, 0.0, 0.0],
    [0.50, 0.0, 0.0], [0.55, 0.0, 0.0],
    [0.60, 0.0, 0.0], [0.65, 0.7333, 0.7333],
    [0.70, 0.9333, 0.9333], [0.75, 1.0, 1.0],
    [0.80, 1.0, 1.0], [0.85, 1.0, 1.0],
    [0.90, 0.8667, 0.8667], [0.95, 0.80, 0.80],
    [1.0, 0.80, 0.80]],
    'green': [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0],
    [0.10, 0.0, 0.0], [0.15, 0.0, 0.0],
    [0.20, 0.0, 0.0], [0.25, 0.4667, 0.4667],
    [0.30, 0.6000, 0.6000], [0.35, 0.6667, 0.6667],
    [0.40, 0.6667, 0.6667], [0.45, 0.6000, 0.6000],
    [0.50, 0.7333, 0.7333], [0.55, 0.8667, 0.8667],
    [0.60, 1.0, 1.0], [0.65, 1.0, 1.0],
    [0.70, 0.9333, 0.9333], [0.75, 0.8000, 0.8000],
    [0.80, 0.6000, 0.6000], [0.85, 0.0, 0.0],
    [0.90, 0.0, 0.0], [0.95, 0.0, 0.0],
    [1.0, 0.80, 0.80]],
    'blue': [[0.0, 0.1, 0.1], [0.05, 0.5333, 0.5333],
    [0.10, 0.6000, 0.6000], [0.15, 0.6667, 0.6667],
    [0.20, 0.8667, 0.8667], [0.25, 0.8667, 0.8667],
    [0.30, 0.8667, 0.8667], [0.35, 0.6667, 0.6667],
    [0.40, 0.5333, 0.5333], [0.45, 0.0, 0.0],
    [0.5, 0.0, 0.0], [0.55, 0.0, 0.0],
    [0.60, 0.0, 0.0], [0.65, 0.0, 0.0],
    [0.70, 0.0, 0.0], [0.75, 0.0, 0.0],
    [0.80, 0.0, 0.0], [0.85, 0.0, 0.0],
    [0.90, 0.0, 0.0], [0.95, 0.0, 0.0],
    [1.0, 0.80, 0.80]],
};

function toHex(x) {             //For each array element
    x = parseInt(x * 255).toString(16);      //Convert to a base16 string
    return (x.length == 1) ? "0" + x : x;  //Add zero if we get only one character
}

function getColorValue(x) {
    // var r =1; g=0; b = 0;

    idx = Math.floor(x * 20);
    mod = (x * 20) - idx;


    if (idx >= _nipy_spectral_data.red.length-1) {
        r = _nipy_spectral_data.red[idx][2]
        g = _nipy_spectral_data.green[idx][2]
        b = _nipy_spectral_data.blue[idx][2]
    }
    else {
        r = (_nipy_spectral_data.red[idx][2] * (1 - mod) + _nipy_spectral_data.red[idx + 1][1] * (mod));
        g = (_nipy_spectral_data.green[idx][2] * (1 - mod) + _nipy_spectral_data.green[idx + 1][1] * (mod));
        b = (_nipy_spectral_data.blue[idx][2] * (1 - mod) + _nipy_spectral_data.blue[idx + 1][1] * (mod));
    }

    var color = '#' + [r, g, b].map(toHex).join("");
    return color;
}

