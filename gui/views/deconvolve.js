var assignments;    // array of Assignment objects
var spotCounts;     // list
var minPositiveSpotCount;
var minPositiveSpotCountSaved;
var minCoverage;
var minCoverageSaved;
var statisticalDeconvolutionMethod = 'cem';
var statisticalDeconvolutionMethodSaved = 'cem';
var spotCountChart;
var deconvolutionResults;
var positiveWells;
var inputSpotCountsBarPlotSorted = false;

function onLoadBody() {}

function onChangeMinCoverage() {
    minCoverage = document.getElementById('input-min-coverage').value;
}

function onChangeMinPositivePoolSpotCount() {
    minPositiveSpotCount = document.getElementById('input-min-spot-count').value;
}

function onChangeStatisticalDeconvolutionMethod() {
    statisticalDeconvolutionMethod = document.getElementById('input-statistical-deconvolution-mode').value;
}

function renderSpotCountsBarPlot(sort) {
    document.getElementById('div-chart').style.display = 'block';
    var plateWellIds = [];
    var spotCountValues = [];
    var spotCounts_ = [];
    for (var i = 0; i < spotCounts.length; i++) {
        spotCounts_.push(spotCounts[i]);
    }
    if (sort == true) {
        spotCounts_ = spotCounts_.sort(({spot_count:a}, {spot_count:b}) => b-a);
    }
    for (var i = 0; i < spotCounts_.length; i++) {
        plateWellIds.push(String(spotCounts_[i].plate_id) + '-' + String(spotCounts_[i].well_id));
        spotCountValues.push(spotCounts_[i].spot_count);
    }
    const ctx = document.getElementById("chart").getContext('2d');
    if (spotCountChart != null) {
        spotCountChart.destroy();
    }
    spotCountChart = new Chart(ctx, {
        type: 'bar',
        data: {
          labels: plateWellIds,
          datasets: [{
            label: 'Spot Count',
            backgroundColor: 'rgb(238, 42, 123, 1)',
            borderColor: 'rgb(238, 42, 123)',
            data: spotCountValues,
          }]
        },
        options: {
          indexAxis: 'y',
          plugins: {
            legend: {
              display: false
            }
          },
          scales: {
            yAxes: [{
              ticks: {
                beginAtZero: true,
              }
            }]
          },
        },
      });
}

function addConfigRow(peptide_id, peptide_sequence, plate_id, well_id) {
    assignments.push(new Assignment(peptide_id, peptide_sequence, -1, -1, plate_id, well_id));
    var table = document.getElementById("configuration-table");
    var rowCount = table.rows.length;
    var row = table.insertRow(rowCount);

    // Peptide index
    var newCell	= row.insertCell(0);
    newCell.innerHTML = "<p name='assignment-index' class=\"text-center my-auto\">" + rowCount + "</p>";

    // Peptide ID
    var newCell	= row.insertCell(1);
    newCell.innerHTML = "<p name='peptide-id' class=\"text-center my-auto\">" + peptide_id + "</p>";

    // Peptide sequence
    var newCell	= row.insertCell(2);
    newCell.innerHTML = "<p name='peptide-sequence' class=\"text-center my-auto font-monospace\">" + peptide_sequence + "</p>";

    // Plate ID
    var newCell	= row.insertCell(3);
    newCell.innerHTML = "<p name='plate-id' class=\"text-center my-auto\">" + plate_id + "</p>";

    // Well ID
    var newCell	= row.insertCell(4);
    newCell.innerHTML = "<p name='well-id' class=\"text-center my-auto\">" + well_id + "</p>";
}

function addSpotCountRow(plate_id, well_id, spot_count) {
    spotCounts.push(new SpotCount(plate_id, well_id, spot_count));
    var table = document.getElementById("spot-counts-table");
    var rowCount = table.rows.length;
    var row = table.insertRow(rowCount);

    // Index
    var newCell	= row.insertCell(0);
    newCell.innerHTML = "<p name='assignment-index' class=\"text-center my-auto\">" + rowCount + "</p>";
    newCell.style = "vertical-align: middle;";

    // Plate ID
    var newCell	= row.insertCell(1);
    newCell.innerHTML = "<p name='plate-id' class=\"text-center my-auto\">" + plate_id + "</p>";

    // Well ID
    var newCell	= row.insertCell(2);
    newCell.innerHTML = "<p name='plate-id' class=\"text-center my-auto\">" + well_id + "</p>";

    // Spot Count
    var newCell	= row.insertCell(3);
    newCell.innerHTML = "<p name='plate-id' class=\"text-center my-auto\">" + spot_count + "</p>";
}

function loadExampleAssignments() {
    assignments = [];
    let filePath = './res/25peptides_5perpool_3x_example_configuration.xlsx';
    fetch(filePath)
        .then(response => response.arrayBuffer())
        .then(data => {
            const workbook = XLSX.read(data, { type: "array" });
            const worksheet = workbook.Sheets['assignment'];

            // Convert worksheet to JSON format
            const excelData = XLSX.utils.sheet_to_json(worksheet, { raw: true });
            for (var i = 0; i < excelData.length; i++) {
                let currPeptideId = excelData[i]['peptide_id'];
                let currPeptideSequence = excelData[i]['peptide_sequence'];
                let currPlateId = excelData[i]['plate_id'];
                let currWellId = excelData[i]['well_id'];
                addConfigRow(currPeptideId, currPeptideSequence, currPlateId, currWellId);
            }
        })
        .catch(error => {
            console.error('Error:', error);
        });
}

function loadExampleSpotCounts() {
    spotCounts = [];
    let filePath = './res/25peptides_5perpool_3x_example_spot_counts.xlsx';
    fetch(filePath)
        .then(response => response.arrayBuffer())
        .then(data => {
            const workbook = XLSX.read(data, { type: "array" });
            const worksheet = workbook.Sheets['Sheet1'];

            // Convert worksheet to JSON format
            const excelData = XLSX.utils.sheet_to_json(worksheet, { raw: true });
            for (var i = 0; i < excelData.length; i++) {
                let currPlateId = excelData[i]['plate_id'];
                let currWellId = excelData[i]['well_id'];
                let currSpotCountValue = excelData[i]['spot_count'];
                addSpotCountRow(currPlateId, currWellId, currSpotCountValue);
                renderSpotCountsBarPlot(false);
            }
        })
        .catch(error => {
            console.error('Error:', error);
        });
}

function renderInputParamsDiv() {
    if (minPositiveSpotCount == 0) {
        document.getElementById('input-min-spot-count').value = '';
    } else {
        document.getElementById('input-min-spot-count').value = minPositiveSpotCount;
    }
    if (minCoverage == 0) {
        document.getElementById('input-min-coverage').value = '';
    } else {
        document.getElementById('input-min-coverage').value = minCoverage;
    }
    document.getElementById('input-statistical-deconvolution-mode').value = statisticalDeconvolutionMethod;
}

function clearConfigTable() {
    var table = document.getElementById('configuration-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }
}

function clearSpotCountsTable() {
    var table = document.getElementById('spot-counts-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }
}

function hideInputSpotCountsBarPlot() {
    sort = false;
    document.getElementById('div-chart').style.display = 'none';
}

function hideResultsDiv() {
    document.getElementById('div-results').style.display = 'none';
}

function reset() {
    assignments = [];
    spotCounts = [];
    minPositiveSpotCount = 0;
    minCoverage = 0;
    statisticalDeconvolutionMethod = 'cem';
    resetConfigFile();
    resetSpotCountsFile();
    clearConfigTable();
    clearSpotCountsTable();
    hideInputSpotCountsBarPlot();
    hideResultsDiv();
    renderInputParamsDiv();
}

function loadExample() {
    reset();
    minPositiveSpotCount = 200;
    minCoverage = 3;
    statisticalDeconvolutionMethod = 'cem';
    loadExampleAssignments();
    loadExampleSpotCounts();
    renderInputParamsDiv();
}

async function downloadExampleConfigFile() {
    if ( window.showSaveFilePicker ) {
        const opts = {
            types: [{
                description: 'ACE_ELISpot_Configuration',
            }],
            suggestedName: 'ace_example_configuration.xlsx',
        };
        var handle = await showSaveFilePicker(opts);
        var writable = await handle.createWritable();
        const downloadLink = document.getElementById('download-example-config-file');
        const fileUrl = downloadLink.getAttribute('href');
        const response = await fetch(fileUrl);
        const fileBlob = await response.blob();
        await writable.write( fileBlob );
        writable.close();
    }
}

async function downloadExampleSpotCountsFile() {
    if ( window.showSaveFilePicker ) {
        const opts = {
            types: [{
                description: 'ACE_ELISpot_Configuration',
            }],
            suggestedName: 'ace_example_spot_counts.xlsx',
        };
        var handle = await showSaveFilePicker(opts);
        var writable = await handle.createWritable();
        const downloadLink = document.getElementById('download-example-spotcounts-file');
        const fileUrl = downloadLink.getAttribute('href');
        const response = await fetch(fileUrl);
        const fileBlob = await response.blob();
        await writable.write( fileBlob );
        writable.close();
    }
}

function onClickInputSpotCountsBarPlotSort() {
    if (inputSpotCountsBarPlotSorted == false) {
        document.getElementById('btn-input-sort-spot-counts-barplot').style.background = 'black';
        document.getElementById('btn-input-sort-spot-counts-barplot').style.color = 'white';
        inputSpotCountsBarPlotSorted = true;
        renderSpotCountsBarPlot(inputSpotCountsBarPlotSorted);
    } else {
        document.getElementById('btn-input-sort-spot-counts-barplot').style.background = null;
        document.getElementById('btn-input-sort-spot-counts-barplot').style.color = null;
        inputSpotCountsBarPlotSorted = false;
        renderSpotCountsBarPlot(inputSpotCountsBarPlotSorted);
    }
}

function resetConfigFile() {
    document.getElementById('configuration-file').value = '';
}

function clearConfigTable() {
    var table = document.getElementById('configuration-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }
}

function loadConfigFile() {
    assignments = [];
    var inputFile = document.getElementById('configuration-file');
    if (inputFile.files.length > 0) {
        clearConfigTable();
        var file = inputFile.files[0];
        let extension = file.name.substring(file.name.lastIndexOf(".") + 1);
        if (extension.toLowerCase() == 'xlsx') {
            var proceed = false;
            try {
                const targetColumnNames1 = ['peptide_id', 'peptide_sequence', 'plate_id', 'well_id'];
                readXlsxFile(file, { sheet: 'assignment' }).then(function(rows) {
                    // Assume that the first row contains column headers.
                    const columnHeaders = rows[0];
                    const columnIndexMap = {};
                    columnHeaders.forEach((columnName, columnIndex) => {
                        columnIndexMap[columnName] = columnIndex;
                    });
                    const targetColumnIndexes = targetColumnNames1.map((targetColumnName) => {
                        return columnIndexMap[targetColumnName];
                    });
                    if (targetColumnIndexes[0] == undefined ||
                        targetColumnIndexes[1] == undefined ||
                        targetColumnIndexes[2] == undefined ||
                        targetColumnIndexes[3] == undefined) {
                        throw new Error('Expected columns do not exist.')
                    }

                    var i = 0;
                    rows.map((row, index) => {
                        if (i > 0) { // data
                            let currPeptideId = row[targetColumnIndexes[0]];
                            let currPeptideSequence = row[targetColumnIndexes[1]];
                            let currPlateId = row[targetColumnIndexes[2]];
                            let currWellId = row[targetColumnIndexes[3]];
                            addConfigRow(
                                currPeptideId,
                                currPeptideSequence,
                                currPlateId,
                                currWellId
                            )
                        }
                        i++;
                    });
                }).catch((error) => {
                    resetConfigFile();
                    document.getElementById('alert-message').innerHTML = 'Please choose a valid ACE configuration file (.xlsx) file. ' +
                        'Expected sheets: "assignment" and "parameters". ' +
                        'Expected columns in "assignment": "peptide_id", "peptide_sequence", "plate_id", "well_id". ' +
                        'Expected columns in "parameters": "num_coverage"';
                    new bootstrap.Modal(document.querySelector("#alert-modal")).show();
                });

                const targetColumnNames2 = ['num_coverage'];
                readXlsxFile(file, { sheet: 'parameters' }).then(function(rows) {
                    // Assume that the first row contains column headers.
                    const columnHeaders = rows[0];
                    const columnIndexMap = {};
                    columnHeaders.forEach((columnName, columnIndex) => {
                        columnIndexMap[columnName] = columnIndex;
                    });
                    const targetColumnIndexes = targetColumnNames2.map((targetColumnName) => {
                        return columnIndexMap[targetColumnName];
                    });
                    var i = 0;
                    rows.map((row, index) => {
                        if (i > 0) { // data
                            minCoverage = row[targetColumnIndexes[0]];
                            renderInputParamsDiv();
                        }
                        i++;
                    });
                }).catch((error) => {
                    resetConfigFile();
                    document.getElementById('alert-message').innerHTML = 'Please choose a valid ACE configuration file (.xlsx) file. ' +
                        'Expected sheets: "assignment" and "parameters". ' +
                        'Expected columns in "assignment": "peptide_id", "peptide_sequence", "plate_id", "well_id". ' +
                        'Expected columns in "parameters": "num_coverage"';
                    new bootstrap.Modal(document.querySelector("#alert-modal")).show();
                });
            } catch (exception) {
                document.getElementById('alert-message').innerHTML = 'Please choose a valid ACE configuration file (.xlsx) file. ' +
                    'Expected sheets: "assignment" and "parameters". ' +
                    'Expected columns in "assignment": "peptide_id", "peptide_sequence", "plate_id", "well_id". ' +
                    'Expected columns in "parameters": "num_coverage"';
                new bootstrap.Modal(document.querySelector("#alert-modal")).show();
            }
        } else {
            document.getElementById('alert-message').innerHTML = 'Please choose a valid ACE configuration file (.xlsx) file. ' +
                'Expected sheets: "assignment" and "parameters". ' +
                'Expected columns in "assignment": "peptide_id", "peptide_sequence", "plate_id", "well_id". ' +
                'Expected columns in "parameters": "num_coverage"';
            new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        }
    }
}

function resetSpotCountsFile() {
    document.getElementById('spot-count-file').value = '';
}

function loadSpotCountsFile() {
    spotCounts = [];
    var inputFile = document.getElementById('spot-count-file');
    if (inputFile.files.length > 0) {
        clearSpotCountsTable();
        var file = inputFile.files[0];
        let extension = file.name.substring(file.name.lastIndexOf(".") + 1);
        if (extension.toLowerCase() == 'xlsx') {
            try {
                const targetColumnNames = ['plate_id', 'well_id', 'spot_count'];
                readXlsxFile(file).then(function(rows) {
                    // Assume that the first row contains column headers.
                    const columnHeaders = rows[0];
                    const columnIndexMap = {};
                    columnHeaders.forEach((columnName, columnIndex) => {
                        columnIndexMap[columnName] = columnIndex;
                    });
                    const targetColumnIndexes = targetColumnNames.map((targetColumnName) => {
                        return columnIndexMap[targetColumnName];
                    });
                    if (targetColumnIndexes[0] == undefined ||
                        targetColumnIndexes[1] == undefined ||
                        targetColumnIndexes[2] == undefined) {
                        throw new Error('Expected columns do not exist.')
                    }
                    var i = 0;
                    rows.map((row, index) => {
                        if (i > 0) { // data
                            let currPlateId = row[targetColumnIndexes[0]];
                            let currWellId = row[targetColumnIndexes[1]];
                            let currSpotCount = row[targetColumnIndexes[2]];
                            addSpotCountRow(
                                currPlateId,
                                currWellId,
                                currSpotCount
                            )
                        }
                        i++;
                        renderSpotCountsBarPlot(false);
                    });
                }).catch((error) => {
                    resetSpotCountsFile();
                    document.getElementById('alert-message').innerHTML = 'Please choose a valid pool spot counts (.xlsx) file (expected columns: "plate_id", "well_id", "spot_count").';
                    new bootstrap.Modal(document.querySelector("#alert-modal")).show();
                });
            } catch (exception) {
                document.getElementById('alert-message').innerHTML = 'Please choose a valid spot counts file (.xlsx) file (expected columns: "plate_id", "well_id", "spot_count").';
                new bootstrap.Modal(document.querySelector("#alert-modal")).show();
            }
        } else if (extension.toLowerCase() == 'csv') {
            try {
                const fileReader = new FileReader();
                fileReader.onload = function (event) {
                    const csvData = event.target.result;
                    const rows = csvData.split('\n');
                    var firstRow = true;
                    for (const row of rows) {
                        if (firstRow == true) {
                            firstRow = false;
                        } else {
                            const columns = row.split(',');
                            addSpotCountRow(columns[0], columns[1], columns[2]);
                        }
                    }
                    renderSpotCountsBarPlot(false);
                };
                fileReader.readAsText(file);
            } catch (exception) {
                document.getElementById('alert-message').innerHTML = 'Please choose a valid spot counts file (.csv) file (expected columns: "plate_id", "well_id", "spot_count").';
                new bootstrap.Modal(document.querySelector("#alert-modal")).show();
            }
        }
    }
}

function deconvolve() {
    if (isReady() == true) {
        minPositiveSpotCountSaved = minPositiveSpotCount;
        minCoverageSaved = minCoverage;
        statisticalDeconvolutionMethodSaved = statisticalDeconvolutionMethod;
        positiveWells = [];
        for (var i = 0; i < spotCounts.length; i++) {
            if (spotCounts[i].spot_count >= minPositiveSpotCount) {
                positiveWells.push(spotCounts[i]);
            }
        }
        eel.deconvolve(
            assignments,
            spotCounts,
            statisticalDeconvolutionMethod,
            minCoverage,
            minPositiveSpotCount
        )(renderOutputDiv)
    }
}

function renderOutputDiv(output) {
    document.getElementById('div-results').style.display = "block";
    document.getElementById('div-results').focus();

    // Save deconvolution results
    const peptideIdsMap = new Map(Object.entries(output['peptide_id']));
    const peptideSequencesMap = new Map(Object.entries(output['peptide_sequence']));
    const peptideSpotCountMap = new Map(Object.entries(output['estimated_peptide_spot_count']));
    const hitWellIdsMap = new Map(Object.entries(output['hit_well_ids']));
    const hitWellIdsCountMap = new Map(Object.entries(output['hit_well_ids_count']));
    const deconvolutionResultsMap = new Map(Object.entries(output['deconvolution_result']));
    deconvolutionResults = [];
    for (const [key, value] of peptideIdsMap.entries()) {
        let currPeptideId = value;
        const currPeptideSequence = peptideSequencesMap.get(key);
        const currPeptideSpotCount = peptideSpotCountMap.get(key);
        const currWellIdsCount = hitWellIdsCountMap.get(key);
        const currWellIds = hitWellIdsMap.get(key);
        const currResult = deconvolutionResultsMap.get(key);
        var deconvolutionResult = new DeconvolutionResult(
            currPeptideId,
            currPeptideSequence,
            currWellIds,
            currWellIdsCount,
            currPeptideSpotCount,
            currResult
        );
        deconvolutionResults.push(deconvolutionResult);
    }

    // Clear hit wells table
    var table = document.getElementById('positive-wells-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }

    // Add hit wells to table
    for (var i = 0; i < positiveWells.length; i++) {
        rowCount = table.rows.length;
        var row = table.insertRow(rowCount);

        // Index
        var newCell	= row.insertCell(0);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + rowCount + "</p>";
        newCell.style = "vertical-align: middle;";

        // Plate ID
        var newCell	= row.insertCell(1);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + positiveWells[i].plate_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Well ID
        var newCell	= row.insertCell(2);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + positiveWells[i].well_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Spot Count
        var newCell	= row.insertCell(3);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + positiveWells[i].spot_count + "</p>";
        newCell.style = "vertical-align: middle;";
    }

    // Clear hit peptides table
    var table = document.getElementById('positive-peptides-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }

    // Add hit wells to table
    for (var i = 0; i < deconvolutionResults.length; i++) {
        if (deconvolutionResults[i].result == 'confident_hit' || deconvolutionResults[i].result == 'candidate_hit') {
            rowCount = table.rows.length;
            var row = table.insertRow(rowCount);

            // Index
            var newCell	= row.insertCell(0);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + rowCount + "</p>";
            newCell.style = "vertical-align: middle;";

            // Peptide ID
            var newCell	= row.insertCell(1);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + deconvolutionResults[i].peptide_id + "</p>";
            newCell.style = "vertical-align: middle;";

            // Peptide sequence
            var newCell	= row.insertCell(2);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + deconvolutionResults[i].peptide_sequence + "</p>";
            newCell.style = "vertical-align: middle;";

            // Hit well IDs
            var newCell	= row.insertCell(3);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + deconvolutionResults[i].hit_well_ids + "</p>";
            newCell.style = "vertical-align: middle;";

            // Hit well IDs count
            var newCell	= row.insertCell(4);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + deconvolutionResults[i].hit_well_ids_count + "</p>";
            newCell.style = "vertical-align: middle;";

            // Estimated peptide spot count
            var newCell	= row.insertCell(5);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + deconvolutionResults[i].peptide_spot_count.toFixed(0) + "</p>";
            newCell.style = "vertical-align: middle;";

            // Deconvolution result
            var resultLabel = ''
            if (deconvolutionResults[i].result == 'confident_hit') {
                resultLabel = 'confident hit';
            } else {
                resultLabel = 'candidate hit';
            }

            var newCell	= row.insertCell(6);
            newCell.innerHTML = "<p class=\"text-center my-auto\">" + resultLabel + "</p>";
            newCell.style = "vertical-align: middle;";
        }
    }
}

async function saveResultsFile() {
    if ( window.showSaveFilePicker ) {
        var zip = new JSZip();

        // Step 1. Create deconvolution results CSV file
        var csvRows1 = ["peptide_id,peptide_sequence,hit_well_ids,hit_well_ids_count,estimated_peptide_spot_count,deconvolution_result"];
        for (var i = 0; i < deconvolutionResults.length; i++) {
            var values = [];
            values.push(deconvolutionResults[i].peptide_id);
            values.push(deconvolutionResults[i].peptide_sequence);
            values.push(deconvolutionResults[i].hit_well_ids);
            values.push(deconvolutionResults[i].hit_well_ids_count);
            values.push(deconvolutionResults[i].peptide_spot_count);
            values.push(deconvolutionResults[i].result);
            csvRows1.push(values.join(','));
        }
        let csvFile1Content = csvRows1.join('\n');
        zip.file("ace_elispot_deconvolution_results.csv", csvFile1Content);

        // Step 2. Create positive well IDs CSV file
        var csvRows2 = ["plate_id,well_id,spot_count"];
        for (var i = 0; i < positiveWells.length; i++) {
            var values = [];
            values.push(positiveWells[i].plate_id);
            values.push(positiveWells[i].well_id);
            values.push(positiveWells[i].spot_count);
            csvRows2.push(values.join(','));
        }
        let csvFile2Content = csvRows2.join('\n');
        zip.file("ace_elispot_positive_wells.csv", csvFile2Content);

        // Step 3. Create ELISpot configuration parameters CSV file
        var csvHeader = [
            'min_positive_well_spot_count',
            'min_coverage',
            'statistical_method'
        ]
        var csvRowContent = [
            minPositiveSpotCountSaved,
            minCoverageSaved,
            statisticalDeconvolutionMethodSaved
        ]
        var csvRows3 = [
            csvHeader.join(','),
            csvRowContent.join(',')
        ]
        let csvFile3Content = csvRows3.join('\n');
        zip.file("ace_elispot_deconvolution_parameters.csv", csvFile3Content);

        // Step 4. Create ELISpot configuration XLSX file
        var wb = XLSX.utils.book_new();

        // Sheet 1 - Deconvolution
        var sheetData1 = [
            ['peptide_id',
             'peptide_sequence',
             'hit_well_ids',
             'hit_well_ids_count',
             'estimated_peptide_spot_count',
             'deconvolution_result']
        ];
        for (var i = 0; i < deconvolutionResults.length; i++) {
            sheetData1.push([
                deconvolutionResults[i].peptide_id,
                deconvolutionResults[i].peptide_sequence,
                deconvolutionResults[i].hit_well_ids,
                deconvolutionResults[i].hit_well_ids_count,
                deconvolutionResults[i].peptide_spot_count,
                deconvolutionResults[i].result
            ]);
        }

        // Sheet 2 - Positive wells
        var sheetData2 = [
            ['plate_id',
             'well_id',
             'spot_count']
        ];
        for (var i = 0; i < positiveWells.length; i++) {
            sheetData2.push([
                positiveWells[i].plate_id,
                positiveWells[i].well_id,
                positiveWells[i].spot_count,
            ]);
        }

        // Sheet 3 - Parameters
        var sheetData3 = [
            ['min_positive_well_spot_count',
             'min_coverage',
             'statistical_method'
            ],
            [
                minPositiveSpotCountSaved,
                minCoverageSaved,
                statisticalDeconvolutionMethodSaved
            ]
        ]

        // Create worksheet instances
        const ws1 = XLSX.utils.aoa_to_sheet(sheetData1);
        const ws2 = XLSX.utils.aoa_to_sheet(sheetData2);
        const ws3 = XLSX.utils.aoa_to_sheet(sheetData3);

        // Add worksheets to the workbook
        XLSX.utils.book_append_sheet(wb, ws1, 'deconvolution_results');
        XLSX.utils.book_append_sheet(wb, ws2, 'positive_wells');
        XLSX.utils.book_append_sheet(wb, ws3, 'parameters');

        // Generate a blob from the workbook
        const wbArrayBuffer = XLSX.write(wb, { bookType: 'xlsx', type: 'array' });

        // Convert the ArrayBuffer to a Blob
        const wbBlob = new Blob([wbArrayBuffer], {
            type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        });

        zip.file("ace_elispot_deconvolution_results.xlsx", wbBlob);

        var zipContent = await zip.generateAsync({ type: "blob" });
        const opts = {
            types: [{
                description: 'ACE_ELISpot_Deconvolution',
                accept: {'text/plain': ['.zip']},
            }],
            suggestedName: 'ace_elispot_deconvolution_results',
        };
        var handle = await showSaveFilePicker(opts);
        var writable = await handle.createWritable();
        await writable.write(zipContent);
        writable.close();
    }
}

function isReady() {
    if (assignments == null) {
        document.getElementById('alert-message').innerHTML = 'Please load a ELISpot configuration.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (assignments.length == 0) {
        document.getElementById('alert-message').innerHTML = 'Please load a ELISpot configuration.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (minCoverage == null || minCoverage == '') {
        document.getElementById('alert-message').innerHTML = 'Please specify the minimum coverage.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (minCoverage == 0) {
        document.getElementById('alert-message').innerHTML = 'Minimum coverage cannot be zero. Please specify a non-zero minimum coverage.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (spotCounts == null) {
        document.getElementById('alert-message').innerHTML = 'Please load pool spot counts.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (spotCounts.length == 0) {
        document.getElementById('alert-message').innerHTML = 'Please load pool spot counts.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (minPositiveSpotCount == null) {
        document.getElementById('alert-message').innerHTML = 'Please specify the minimum positive pool (well) spot count.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (minPositiveSpotCount == 0) {
        document.getElementById('alert-message').innerHTML = 'Minimum positive pool (well) cannot be zero. Please specify a non-zero minimum positive pool (well) spot count.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    return true;
}
