from app import app
from flask import render_template, make_response, jsonify, request
import pandas as pd
import numpy as np
import json
import pathlib
from numba import jit, float64, int32, char, prange




#sheet = pd.read_csv(pathlib.PurePath(app.root_path,'data', 'genetics', 'new_sheet.csv'))

@app.route("/")
def index():
    return render_template("indexMain.html")

#@jit(nopython=True)
@jit(nopython=False)
#@jit(float64[:,:,:](float64[:], float64[:], char[:],char[:], int32, int32, float64, int32, int32, int32), nopython=False, forceobj=True, parallel=True)
#@jit(float64[:,:,:](float64[:], float64[:], char[:],char[:], float64, float64, int32, int32, int32, int32), nopython=True, forceobj=False)
def runKDE(X, Y, Zgenes, genes, xmax, ymax, sigma, width, height, nStds):
    vfBuffer = np.zeros(shape=(height, width, len(genes)), dtype=np.float64)
    x = 0
    y = 0
    z = 0
    val = np.float64(0.0)
    n_steps = int(np.ceil(sigma*nStds))
    normalization = 1/(np.power(np.pi*np.power(sigma,2), 0.5))
    #print("Starting LOOP")

    for i in range(0, len(Zgenes)):
        x = int(np.round(X[i] * (height - 1) / xmax))
        y = int(np.round(Y[i] * (width - 1) / ymax))
        if Zgenes[i] in genes:
            z = np.where(genes ==Zgenes[i])
            z=int(z[0][0])
            for m in range(n_steps*-1, n_steps+1):
                for n in range(n_steps*-1, n_steps+1):
                    if ((x + m) > 0) and ((y + n) > 0) and ((x + m) < height - 1) and ((y + n) < width - 1):
                        #val = Math.exp(-(Math.pow(n, 2) + Math.pow(m, 2)) / Math.pow(sigma, 2))*normalization;
                        val=(np.exp(((np.power(n, 2)+np.power(m,2))*-1)/np.power(sigma, 2)))*normalization
                        x1=x+m
                        y1=y+n
                        val1=vfBuffer[x1, y1, z]
                        vfBuffer[x1, y1, z]=val1 + val
                        
    print(vfBuffer.shape)
    return(vfBuffer)

@app.route("/runKDE", methods=["POST"])
def handleKDE():
    data = request.get_json()
    if data is None:
        return jsonify(status="error")
    else:
        data=dict(data)
    X=np.array(data['paramX'], dtype=np.float64)
    Y=np.array(data['paramY'], dtype=np.float64)
    Zgenes=np.array(data['paramZGenes'], dtype=np.string_)
    genes=np.array(data['paramGenes'], dtype=np.string_)
    xmax=np.float64(data['paramXmax'])
    ymax=np.float64(data['paramYmax'])
    sigma=np.float64(data['paramSigma'])
    width=np.int32(data['paramWidth'])
    height=np.int32(data['paramHeight'])

    if 'paramnStds' in data:
        nStds=data['paramnStds']
    else:
        nStds=np.int32(2)
    vfBuffer=runKDE(X, Y, Zgenes, genes, xmax, ymax, sigma, width, height, nStds)
    return jsonify(vfBuffer=vfBuffer.tolist())

@app.route("/assignCelltypes", methods=["POST"])
def assignCelltypes():
    data = request.get_json()
    if data is None:
        return jsonify(status="error")
    vf=np.array(data['paramvf'])
    vfNorm=vf.sum(2)
    signatureMatrix=np.array(data['paramSigMat'])
    
    threshold=data['paramthreshold']
    intermediates=[]
    for i in range(0, vf.shape[0]):
        inter = vf[i,:,:]
        inter = np.log(np.add(inter, 1))
        inter = np.transpose(np.divide(np.transpose(inter), np.sum(inter,1)))
        inter = np.matmul(inter, np.transpose(signatureMatrix))
        intermediates.append(np.expand_dims(np.argmax(inter, 1), 0))
    intermediates=np.concatenate(intermediates)
    inters = np.multiply(intermediates, np.greater(vfNorm, threshold))
    celltypeMap=np.subtract(inters, np.less(vfNorm, threshold))
    return jsonify(celltypeMap=celltypeMap.tolist())
    