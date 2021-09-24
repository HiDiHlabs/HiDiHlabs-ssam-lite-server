
function readFileAsync(file) {
    return new Promise((resolve, reject) => {
        let reader = new FileReader();

        reader.onload = () => {
            resolve(reader.result);
        };

        reader.onerror = reject;

        reader.readAsText(file, "UTF-8");
    })
};

function createSignatureMatrix(rows, current_genes) {

    // console.log(rows);
    // genes = current_genes;
    var nGenes = current_genes.length;
    var nClusters = rows.length;
    var nChecked = 0;

    for (var i = 0; i < nClusters; i++) {
        nChecked += ((rows[i]).visible & (rows[i]).checked)
    }

    console.log(nChecked);

    signatureBuffer = tf.buffer([nChecked, nGenes]);

    clusterLabels = [];

    var n = 0;

    for (var i = 0; i < nClusters; i++) {
        if ((rows[i]).visible & (rows[i]).checked) {
            line = ((rows[i]).expressions)

            for (var j = 0; j < nGenes; j++) {
                signatureBuffer.set(parseFloat(line[j]), n, j);
            }
            n++;
            clusterLabels.push(((rows[i]).celltype));
        }
    }
    signatureBuffer = signatureBuffer.toTensor();

    return [signatureBuffer, clusterLabels];

};

function processSignatures(allText) {

    console.log('...loading...');

    var genes;
    var clusterLabels;
    var signatureBuffer;

    var allTextLines = allText.split(/\r\n|\n/);
    var nClusters = allTextLines.length - 2;
    genes = allTextLines[0].split(',').slice(1);
    var nGenes = genes.length;

    signatureBuffer = tf.buffer([nClusters, nGenes]);

    clusterLabels = [];

    for (var i = 0; i < nClusters; i++) {
        line = allTextLines[i + 1].split(',');

        for (var j = 0; j < nGenes; j++) {
            signatureBuffer.set(parseFloat(line[j + 1]), i, j);
        }

        clusterLabels.push(line[0]);
    }

    signatureBuffer = signatureBuffer.toTensor();

    return [signatureBuffer, clusterLabels, genes];

};

function processCoordinates(allText) {

    // var t0 = performance.now()
    var allTextLines = allText.split(/\r\n|\n/);

    var len = allTextLines.length;
    ZGenes = Array(len);
    X = Array(len);
    Y = Array(len);

    allGenes = [];

    // ZGenes = []; //Array(len);
    // X = [];//Array(len);
    // Y = []//Array(len);
    var x = 0;
    var y = 0;
    var line;
    xmax = 0;
    ymax = 0;

    for (var i = 1; i < allTextLines.length; i++) {
        line = allTextLines[i].split(',');

        x = parseFloat(line[1]) || 0;
        y = parseFloat(line[2]) || 0;

        xmax = Math.max(x, xmax);
        ymax = Math.max(y, ymax);

        X[i] = x;//.push(x);
        Y[i] = y;//.push(y);
        ZGenes[i] = line[0];//.push(line[0]);

        if (!allGenes.includes(line[0])) {
            allGenes.push(line[0]);
        };
    }
    // console.log(performance.now() -t0);

    var edgeRatio = Math.ceil(xmax / ymax);
    var width = Math.ceil(500);
    var height = Math.ceil(width / edgeRatio);

    return [X, Y, ZGenes, allGenes, xmax, ymax, edgeRatio, width, height];
};



function main() {

    {
        var genes = [];
        var allTissues = [];
        var allSpecies = [];
        var signatureContent = [];
        var clusterLabels;
        var signatureMatrix;

        var X;          // mRNA x coordinates
        var Y;          // mRNA y coordinates
        var ZGenes;     // gene information
        var xmax;       // highest coordinate
        var ymax;       // lowest coordinate
        var sigma = 1;    // KDE kernel width

        var height = 300; // vf height (pixels)
        var width;      // vf width (pixels)
        var edgeRatio = 1;// radion height/width
        var threshold = 10;
        var vf;         // tensor vectorfield
        var vfNorm;     // tensor vfNorm

        var parameterWindow = [250, 250];
        var parameterWidth = 50;
        var parameterX = [];
        var parameterY = [];
        var parameterZ = [];
        var vfParameter;
        var vfNormParameter;

        var pointerCoordinates = parameterWindow;

        function getClusterLabel(i) {
            return clusterLabels[i];
        }
    }


    async function importSignatures(path) {

        // console.log('loading signatures');
        document.getElementById("signature-loader").style.display = "block";   //display waiting symbol

        var fileToLoad = document.getElementById("btn-signatures-hidden").files[0];

        [signatureMatrix, clusterLabels, genes] = (processSignatures(await readFileAsync(fileToLoad)));


        plotSignatures('signatures-preview', genes, clusterLabels, signatureMatrix.arraySync()).then(function () {
            document.getElementById("signature-loader").style.display = "none";
        });


    };


    async function downloadSignatures(current_genes) {


        $.ajax({
            type: "POST",
            url: "http://127.0.0.1:5000/genes",
            contentType: "application/json",
            data: JSON.stringify({ "genes": current_genes }),
            dataType: "json",
            success: function (data) {
                console.log('downloaded...');
                current_genes;
                signatureContent = (processDownloadedSignatures(data));
                // console.log(signatureContent);
                createTableHeader(current_genes, allTissues, allSpecies);
                updateTable(signatureContent, allTissues, allSpecies);
                [signatureMatrix, clusterLabels] = createSignatureMatrix(data, current_genes);
                //console.log(signatureMatrix);
            },
            error: function (err) {
                // console.log(err);
                return (err);
            }
        });


    }
    function processDownloadedSignatures(data) {
        allTissues = [];
        allSpecies = [];
        signatureContent = []

        for (var i = 0; i < data.length; i++) {
            entry = data[i];

            [cellType, tissueType, species] = (data[i]).tissue.split('|');
            // console.log(tissueType, allTissues)
            if (!allTissues.includes(tissueType)) { allTissues.push(tissueType) };
            if (!allSpecies.includes(species)) { allSpecies.push(species) };

            entry.label = (data[i]).tissue;
            entry.tissue = tissueType;
            entry.celltype = cellType;
            entry.species = species;
            entry.checked = true;
            entry.visible = true;
            signatureContent.push(entry);
        }

        return signatureContent;
    }

    function createTableHeader(genes, tissues, species) {
        var header = document.getElementById('sheet-header')
        // row = header.insertRow(-1);
        typeCell = header.insertCell(0);

        typeCell = header.insertCell(1);
        typeCell.innerHTML = "Cell Type";
        typeCell.classList.add('header-cell');

        tissueCell = header.insertCell(2);
        tissueCell.innerHTML = "Tissue"
        tissueCell.classList.add('header-cell');

        var selTissue = document.createElement("select");
        selTissue.id = "tissueSelect"
        selTissue.classList.add('selectpicker')
        selTissue.setAttribute('multiple', '');
        selTissue.onchange = updateSignatureContent

        for (var i = 0; i < tissues.length; i++) {
            //Add the options
            // console.log(selTissue);
            selTissue.options[selTissue.options.length] = new Option(tissues[i]);
        }
        tissueCell.appendChild(selTissue);

        speciesCell = header.insertCell(3);
        speciesCell.innerHTML = "Species"
        speciesCell.classList.add('header-cell');


        var selSpecies = document.createElement("select");
        selSpecies.id = "speciesSelect"
        selSpecies.classList.add('selectpicker')
        selSpecies.setAttribute('multiple', '');
        selSpecies.onchange = updateSignatureContent

        for (var i = 0; i < species.length; i++) {
            //Add the options
            // console.log(selSpecies);
            selSpecies.options[selSpecies.options.length] = new Option(species[i]);
        }
        speciesCell.appendChild(selSpecies)

        // console.log(genes);
        for (var g = 0; g < genes.length; g++) {
            cell = header.insertCell(g + 4);
            cell.innerHTML = genes[g]
            cell.classList.add('header-cell');
            cell.classList.add('rotate');
        }

        $('select').selectpicker();
        $("#tissueSelect").selectpicker('val', tissues[0])
        $("#speciesSelect").selectpicker('val', species[0])

    }


    function updateSignatureContent(event) {
        markedTissues = ($("#tissueSelect").val());
        markedSpecies = ($("#speciesSelect").val());
        if (markedTissues != null & markedSpecies != null) {

            for (var i = 0; i < signatureContent.length; i++) {
                // console.log(signatureContent[i]);
                if ((markedSpecies.includes((signatureContent[i]).species)) &
                    (markedTissues.includes((signatureContent[i]).tissue))) {
                    (signatureContent[i]).visible = true;
                }
                else {
                    (signatureContent[i]).visible = false;
                    (signatureContent[i]).checked = true;
                }
            }
            updateTable(signatureContent);

            [signatureMatrix, clusterLabels] = createSignatureMatrix(signatureContent, genes)
        }
    }

    function onlyUnique(value, index, self) {
        return self.indexOf(value) === index;
    }

    async function importCoordinates(path) {

        // console.log(ZGenes);
        // console.log('loading coordinates');
        document.getElementById("coordinate-loader").style.display = "block";   //display waiting symbol

        var fileToLoad = document.getElementById("btn-coordinates-hidden").files[0];

        [X, Y, ZGenes, genes, xmax, ymax, edgeRatio, width, height] = processCoordinates(await readFileAsync(fileToLoad));
        // console.log(genes);

        // var current_genes = ZGenes.filter(onlyUnique);

        downloadSignatures(genes);

        edgeRatio = xmax / ymax;
        width = Math.ceil(height * edgeRatio);
        console.log('plotting...');
        plotCoordinates('coordinates-preview', X, Y, ZGenes).then(function () {
            document.getElementById("coordinate-loader").style.display = "none";
        });
    };

    function runFullKDE() {
        //[vf, vfNorm] = runKDE(X, Y, ZGenes, genes, xmax, ymax, sigma, width, height);
        runFullKDEviaPOST()
        //console.log(X)
        //plotVfNorm('vf-norm-preview', vfNorm.arraySync());
    };

    function runFullKDEviaPOST() {
        $(document).on({
            ajaxStart: function(){
                $("#vf-norm-load").show();
                $("#vf-norm-preview").hide(); 
            },
            ajaxStop: function(){ 
                $("#vf-norm-load").hide();
                $("#vf-norm-preview").show();
            }    
        });
        payload = JSON.stringify({
            paramX: X,
            paramY: Y,
            paramZGenes: ZGenes,
            paramGenes: genes,
            paramXmax: xmax,
            paramYmax: ymax,
            paramSigma: sigma,
            paramWidth: width,
            paramHeight: height
        });
        //console.log(payload);

        $.ajax({
            type: "POST",
            url: "/runKDE",
            contentType: "application/json",
            data: payload,
            dataType: "json",
            success: function (response) {
                vfNorm=response.vfNorm;
                vf=response.vfBuffer;
                plotVfNorm('vf-norm-preview', tf.tensor(vfNorm).arraySync());
                console.log(response.vfNorm);
            },
            error: function (err) {
                console.log(err);
            }
        });
        
    };

    function runCelltypeAssignments() {
        [signatureMatrix, clusterLabels] = createSignatureMatrix(signatureContent, genes);
        celltypeMap = assignCelltypes(tf.tensor(vf), tf.tensor(vfNorm), signatureMatrix, threshold);
        assignCelltypesviaPOST(vf, vfNorm, signatureMatrix, threshold)
        plotCelltypeMap('celltypes-preview', celltypeMap.arraySync(), clusterLabels, getClusterLabel);

    };

    function assignCelltypesviaPOST(vf, vfNorm, signatureMatrix, threshold){
        // $(document).on({
        //     ajaxStart: function(){
        //         $("#vf-norm-load").show();
        //         $("#vf-norm-preview").hide(); 
        //     },
        //     ajaxStop: function(){ 
        //         $("#vf-norm-load").hide();
        //         $("#vf-norm-preview").show();
        //     }    
        // });
        payload = JSON.stringify({
            paramVF: vf,
            paramVFNorm: vfNorm,
            paramZSigMat: signatureMatrix,
            paramthreshold: threshold,
        });
        //console.log(payload);

        $.ajax({
            type: "POST",
            url: "/assignCelltypes",
            contentType: "application/json",
            data: payload,
            dataType: "json",
            success: function (response) {
                celltypeMap=tf.tensor(response.vfNorm);
                //plotCelltypeMap('celltypes-preview', celltypeMap.arraySync(), clusterLabels, getClusterLabel);
                console.log(response.celltypeMap);
            },
            error: function (err) {
                console.log(err);
            }
        });
    }

    function updateVfShape() {
        height = parseInt(document.getElementById('vf-width').value);
        width = Math.ceil(height * edgeRatio);
        document.getElementById("vf-size-information").innerHTML =
            "total size: (" + width + "," + height + "," + genes.length + "); " + Math.ceil(width * height * 32 * genes.length / 1024 ** 2) + " gB"
        updateParameterVf();
    };

    function updateSigma() {
        sigma = parseFloat(document.getElementById('KDE-bandwidth').value);
        updateParameterVf();
    };

    function updateThreshold() {
        threshold = parseFloat(document.getElementById('threshold').value);
        updateParameterCelltypes();
    };

    function toggleParameterGenerator() {
        var previewGenerator = document.getElementById('preview-generator');

        if (previewGenerator.style.display === "none") {
            displayParameterGenerator();
            createParameterCoodinatesPlot();
        } else {
            hideParameterGenerator();
        }

    };



    function updateParameterCelltypes() {
        parameterCelltypeMap = assignCelltypes(vfParameter, vfNormParameter, signatureMatrix, threshold);
        plotCelltypeMap('parameter-celltypes', parameterCelltypeMap.arraySync(), clusterLabels);
    }

    function updateParameterVf() {
        [vfParameter, vfNormParameter] = runKDE(parameterX, parameterY, parameterZ,
            genes, parameterWidth * xmax / width * 2, parameterWidth * ymax / height * 2,
            sigma, parameterWidth * 2, parameterWidth * 2);
        plotVfNorm('parameter-vf', vfNormParameter.arraySync());
        updateParameterCelltypes();
    };

    function updateParameterCoordinates() {
        parameterX = []
        parameterY = []

        var rectCenter = [parameterWindow[0] / width * xmax, parameterWindow[1] / height * ymax];
        var rectEdge = parameterWidth / width * xmax;

        for (var i = 0; i < X.length; i++) {
            if (X[i] > rectCenter[0] - rectEdge &&
                Y[i] > rectCenter[1] - rectEdge &&
                X[i] < rectCenter[0] + rectEdge &&
                Y[i] < rectCenter[1] + rectEdge) {
                parameterX.push(X[i] - rectCenter[0] + rectEdge);
                parameterY.push(Y[i] - rectCenter[1] + rectEdge);
                parameterZ.push(ZGenes[i]);
            }
        }
    };

    function createParameterCoodinatesPlot() {
        updateParameterCoordinates();
        var rectCenter = [parameterWindow[0] / width * xmax, parameterWindow[1] / height * ymax];
        var rectEdge = parameterWidth / width * xmax;
        var layoutRect = {
            shapes: [
                //Unfilled Rectangle
                {
                    type: 'rect',
                    x0: rectCenter[0] - rectEdge,
                    y0: rectCenter[1] - rectEdge,
                    x1: rectCenter[0] + rectEdge,
                    y1: rectCenter[1] + rectEdge,
                    line: {
                        color: 'rgba(255, 255, 255, 1)'
                    },
                },],

        }

        plotCoordinates('parameter-coordinates', X, Y, ZGenes, layoutRect);
        document.getElementById('parameter-coordinates')
            .on('plotly_hover', updatePointerCoordinates);
        document.getElementById('parameter-coordinates')
            .on('plotly_click', updateRectangle);
        updateParameterVf();
    }

    function updatePointerCoordinates(eventData) {
        pointerCoordinates = [eventData.xvals[0], eventData.yvals[0]];
    };

    function updateRectangle(eventData) {
        updateParameterRectangle(pointerCoordinates, parameterWidth * xmax / width);
        parameterWindow = [Math.ceil(pointerCoordinates[1] / xmax * width),
        Math.ceil(pointerCoordinates[0] / xmax * width)];
        updateParameterCoordinates();
        updateParameterVf();
        console.log(pointerCoordinates);
    };

    function initiateButtons() {


        document.getElementById('btn-signatures-hidden')
            .addEventListener('change', importSignatures);
        document.getElementById('btn-coordinates-hidden')
            .addEventListener('change', importCoordinates);
        document.getElementById('btn-KDE')
            .addEventListener('click', runFullKDE);
        document.getElementById('btn-types')
            .addEventListener('click', runCelltypeAssignments);

        document.getElementById('btn-parameters')
            .addEventListener('click', toggleParameterGenerator);
        document.getElementById('vf-width')
            .addEventListener('change', updateVfShape);
        document.getElementById('KDE-bandwidth')
            .addEventListener('change', updateSigma);
        document.getElementById('threshold')
            .addEventListener('change', updateThreshold);
        document.getElementById('button-tutorial')
            .addEventListener('click', runTutorial);

    };

    initiateButtons();

}
main();