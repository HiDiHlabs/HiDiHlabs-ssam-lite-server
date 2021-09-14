from app import app
from flask import render_template, make_response, jsonify, request
import pandas as pd
import numpy as np

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
