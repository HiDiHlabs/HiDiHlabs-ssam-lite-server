from app import app
from flask import render_template, make_response, jsonify, request
import pandas as pd
import numpy as np
import json

sheet = pd.read_csv('../../data/genetics/new_sheet.csv')

@app.route("/")
def index():
    return render_template("indexMain.html")

@app.route("/genes", methods=["POST"])
def get_celltypes():    

    req = request.get_json(force=True)
    genes = req['genes']

    # response = {}
    all_cells = {}#[]

    for g,gene in enumerate(genes):
        sublist = (sheet.loc[sheet['cellMarker'] == gene])
        cell_types = list(sublist['cellName'])
        tissue_types = list(sublist['tissueType'])
        species = list(sublist['speciesType'])

        labels = []
        
        for c in cell_types:
            for t in tissue_types:
                for s in species:
                    labels.append('|'.join([c,t,s]))

        for t,label in enumerate(labels):

            if label in all_cells.keys():
                (all_cells[label])[g]+=1
                # print(all_cells)
            else:
                all_cells[label]=np.zeros((len(genes),))
                (all_cells[label])[g]=1

    tissues = list(all_cells.keys())

    exp_mat = np.array([all_cells[t] for t in all_cells.keys()])

    exp_mat = exp_mat/(exp_mat.sum(0)+0.1)
    
    exp_mat/=exp_mat.max()
    exp_mat = exp_mat[exp_mat.sum(1)>0.001]

    idcs = exp_mat.sum(1).argsort()[::-1]

    print(idcs.mean())
    
    response = []
    for i,idx in enumerate(idcs):
        # if idcs[i]:
        response.append({'tissue' : tissues[idx],
        'expressions' : list(exp_mat[idx]) })

    res = make_response(jsonify(response), 201)
    # for r in list(genes):
    return(res)


@app.route("/runKDE", methods=["POST"])
def runKDE():
    #print('hello')
    data = request.get_json()
    if data is None:
        return jsonify(status="error")
    X=data['paramX'][1:]
    Y=data['paramY'][1:]
    Zgenes=data['paramZGenes'][1:]
    genes=data['paramGenes']
    xmax=data['paramXmax']
    ymax=data['paramYmax']
    sigma=data['paramSigma']
    width=data['paramWidth']
    height=data['paramHeight']
    if 'paramnStds' in data:
        nStds=data['paramnStds']
    else:
        nStds=3
    


    vfBuffer = np.zeros(shape=(height, width, len(genes)))

    x = 0
    y = 0
    z = 0
    val = 0
    counter = 0
    n_steps = int(np.ceil(sigma*nStds))

    normalization = 1/(1*np.pi*sigma**2)**0.5
    #print(type(X))
    for i in range(0, len(Zgenes)):
        #print(i)
        #print(X[i+1])
        x = np.round(X[i] * (height - 1) / xmax)
        y = np.round(Y[i] * (width - 1) / ymax)
        z = genes.index(Zgenes[i])
        if (z >= 0):
            for m in range(n_steps*-1, n_steps+1):
                for n in range(n_steps*-1, n_steps+1):
                    if (((x + m) > 0) and ((y + n) > 0) and ((x + m) < height - 1) and ((y + n) < width - 1)):
                        val = np.exp(((np.power(n, 2) + np.power(m, 2))*-1) / np.power(sigma, 2))*normalization
                        vfBuffer[int(x + m), int(y + n), z]=vfBuffer[int(x + m), int(y + n), int(z)] + val
            counter+=1
        else:
            print(Zgenes[i], i, z)
    vfNorm=vfBuffer.sum(2)
    vfBuffer=vfBuffer.tolist()
    vfNorm=vfNorm.tolist()
    return jsonify(vfBuffer=vfBuffer, vfNorm=vfNorm)

@app.route("/assignCelltypes", methods=["POST"])
def assignCelltypes():
    data = request.get_json()
    if data is None:
        return jsonify(status="error")
    print(data['paramVF'])
    vf=np.asarray(data['paramVF'])
    vfNorm=np.asarray(data['paramVFNorm'])
    signatureMatrix=data['paramZSigMat']
    threshold=data['paramthreshold']
    
    intermediates = []
    inter=[]
    print(vf.shape)
    for i in range(0, vf.shape[0]):
        inter = vf[i]
        inter = np.log(inter+1)
        inter = np.divide(np.transpose(inter.transpose()), np.transpose(inter.sum(1)))
        inter = np.matmul(inter, np.transpose(signatureMatrix))
        intermediates.append(np.expand_dims(np.argmax(inter, 1), 0))
        

    intermediates = np.concat(intermediates)
    inters = np.multiply(intermediates, (np.greater(vfNorm,threshold)))

    celltypeMap=np.substract(inters, (np.less(vfNorm,threshold)))

    return celltypeMap