import os
import pathlib
import urllib.request
base_path=pathlib.Path(__file__).parent.absolute()
data_directory=pathlib.PurePath(base_path, "app", "data", "genetics")
js_directory=pathlib.PurePath(base_path, "app", "static", "pkgs", "js")
bootstrap_directory=pathlib.PurePath(base_path, "app", "static", "pkgs", "bootstrap")

# create directories if they don't exist
pathlib.Path(js_directory).mkdir(parents=True, exist_ok=True)
pathlib.Path(bootstrap_directory).mkdir(parents=True, exist_ok=True)
pathlib.Path(data_directory).mkdir(parents=True, exist_ok=True)

urllib.request.urlretrieve("https://code.jquery.com/jquery-3.4.1.min.js", pathlib.PurePath(js_directory,"jquery-3.4.1.min.js" ))
urllib.request.urlretrieve("https://unpkg.com/@popperjs/core@2.10.1/dist/umd/popper.min.js", pathlib.PurePath(js_directory,"popper.min.js" ))
urllib.request.urlretrieve("https://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js", pathlib.PurePath(js_directory,"jquery.min.js" ))


urllib.request.urlretrieve("https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css", pathlib.PurePath(bootstrap_directory,"bootstrap.min.css" ))
urllib.request.urlretrieve("https://cdnjs.cloudflare.com/ajax/libs/bootstrap-select/1.13.1/css/bootstrap-select.css", pathlib.PurePath(bootstrap_directory,"bootstrap-select.css" ))
urllib.request.urlretrieve("https://cdnjs.cloudflare.com/ajax/libs/bootstrap-multiselect/0.9.13/css/bootstrap-multiselect.css", pathlib.PurePath(bootstrap_directory,"bootstrap-multiselect.css" ))



urllib.request.urlretrieve("https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/js/bootstrap.bundle.min.js",pathlib.PurePath(js_directory,"bootstrap.bundle.min.js" ))
urllib.request.urlretrieve("https://cdnjs.cloudflare.com/ajax/libs/bootstrap-select/1.13.1/js/bootstrap-select.min.js",pathlib.PurePath(js_directory,"bootstrap-select.min.js" ))

urllib.request.urlretrieve("https://cdn.jsdelivr.net/npm/@tensorflow/tfjs@2.0.0/dist/tf.min.js",pathlib.PurePath(js_directory,"tf.min.js" ))
urllib.request.urlretrieve("https://cdn.plot.ly/plotly-1.58.4.min.js",pathlib.PurePath(js_directory,"plotly-1.58.4.min.js" ))


print("All files retrieved \n\n\t You can now start the server with command")