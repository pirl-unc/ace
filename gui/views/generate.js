var peptides;
var numPeptides = 0;
var numPeptidesPerPool = 0;
var numCoverage = 0;
var plateSize = 96;
var clusterPeptides = 'yes';
var sequencesAvailable = false;
var sequenceSimilarityFxn = 'euclidean';
var sequenceSimilarityThreshold = 0.7;
var maxIterations = 2000;
var initializationStrategy = 'greedy';
var allowExtraPools = 'no';
var numTotalPools = 0;
var numTotalPlates = 0;
var expandedAccordion = 'sequences';
var assignments;
var assignmentsBenchReady;
var preferredPeptidePairs;


function onClickSpecifySequences() {
    expandedAccordion = 'sequences';
}

function onClickSpecifyCount() {
    expandedAccordion = 'count';
}

function onChangePeptideCount() {
    numPeptides = document.getElementById('num-peptides').value;
    updateTotalCounts();
}

function onChangePeptidesPerPool() {
    numPeptidesPerPool = document.getElementById('num-peptides-per-pool').value;
    updateTotalCounts();
}

function onChangeCoverage() {
    numCoverage = document.getElementById('num-coverage').value;
    updateTotalCounts();
}

function onChangePlateSize() {
    plateSize = document.getElementById('num-wells').value;
    updateTotalCounts();
}

function addPeptideRow(peptideId, peptideSequence) {
    peptides.push(new Peptide(peptideId, peptideSequence));

    var table = document.getElementById("peptides-table");
    var rowCount = table.rows.length;
    var row = table.insertRow(rowCount);

    // Peptide index
    var newCell	= row.insertCell(0);
    newCell.innerHTML = "<p name='peptide-index' class=\"text-center my-auto\">" + rowCount + "</p>";

    // Peptide ID
    var newCell	= row.insertCell(1);
    newCell.innerHTML = "<p name='peptide-id' class=\"text-center my-auto\">" + peptideId + "</p>"

    // Peptide sequence
    var newCell	= row.insertCell(2);
    newCell.innerHTML = "<p name='peptide-sequence' class=\"text-center my-auto font-monospace\">" + peptideSequence + "</p>"

    updateTotalCounts();
}

function clearTable() {
    var table = document.getElementById('peptides-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }
}

function expandSequencesAccordion() {
    expandedAccordion = 'sequences';
    document.getElementById('collapseOne').classList.add('show');
    document.getElementById('collapseOne').classList.remove('collapsed');
    document.getElementById('accordion-1-btn').classList.add('show');
    document.getElementById('accordion-1-btn').classList.remove('collapsed');

    document.getElementById('collapseTwo').classList.add('collapsed');
    document.getElementById('collapseTwo').classList.remove('show');
    document.getElementById('accordion-2-btn').classList.remove('show');
    document.getElementById('accordion-2-btn').classList.add('collapsed');
}

function loadExamplePeptides() {
    peptides = [];
    clearTable();
    let filePath = './res/25peptides_example.xlsx';
    fetch(filePath)
        .then(response => response.arrayBuffer())
        .then(data => {
            const workbook = XLSX.read(data, { type: "array" });
            const worksheet = workbook.Sheets['Sheet1'];

            // Convert worksheet to JSON format
            const excelData = XLSX.utils.sheet_to_json(worksheet, { raw: true });
            for (var i = 0; i < excelData.length; i++) {
                let currId = excelData[i]['peptide_id'];
                let currSequence = excelData[i]['peptide_sequence'];
                addPeptideRow(currId, currSequence);
            }
        })
        .catch(error => {
            console.error('Error:', error);
        });
}

function renderInputParamsDiv() {
    if (numPeptides == 0) {
        document.getElementById("num-peptides").value = '';
    }
    if (numPeptidesPerPool == 0) {
        document.getElementById("num-peptides-per-pool").value = '';
    } else {
        document.getElementById("num-peptides-per-pool").value = numPeptidesPerPool;
    }
    if (numCoverage == 0) {
        document.getElementById("num-coverage").value = '';
    } else {
        document.getElementById("num-coverage").value = numCoverage;
    }
    document.querySelector('#num-wells').value = plateSize;
    document.querySelector('#cluster-peptides').value = clusterPeptides;
    document.querySelector('#sequence-similarity-fxn').value = sequenceSimilarityFxn;
    document.getElementById("sequence-similarity-threshold").value = sequenceSimilarityThreshold;
    document.getElementById("threshold-label").innerHTML = sequenceSimilarityThreshold;
    document.getElementById("num-iterations").value = maxIterations;
    document.querySelector('#strategy').value = initializationStrategy;
    document.querySelector('#allow-extra-pools').value = allowExtraPools;
}

function updateTotalCounts() {
    if (expandedAccordion == 'sequences') {
        numTotalPools = Math.ceil(peptides.length / numPeptidesPerPool) * numCoverage;
    } else {
        numTotalPools = Math.ceil(numPeptides / numPeptidesPerPool) * numCoverage;
    }
    document.getElementById('num-total-pools').innerHTML = numTotalPools;

    if (plateSize != 'do_not_assign') {
        numTotalPlates = Math.ceil(numTotalPools / plateSize);
        document.getElementById('num-total-plates').innerHTML = numTotalPlates;
    } else {
        document.getElementById('num-total-plates').innerHTML = 'NaN';
    }
}

function loadExample() {
    numPeptidesPerPool = 5;
    numCoverage = 3;
    plateSize = 96;
    clusterPeptides = 'yes';
    sequenceSimilarityFxn = 'euclidean';
    sequenceSimilarityThreshold = 0.7;
    initializationStrategy = 'greedy';
    allowExtraPools = 'no';
    expandSequencesAccordion();
    loadExamplePeptides();
    renderInputParamsDiv();
    updateTotalCounts();
}

function reset() {
    numPeptides = 0;
    numPeptidesPerPool = 0;
    numCoverage = 0;
    plateSize = 96;
    clusterPeptides = 'yes';
    sequenceSimilarityFxn = 'euclidean';
    sequenceSimilarityThreshold = 0.7;
    maxIterations = 2000;
    initializationStrategy = 'greedy';
    allowExtraPools = 'no';
    clearTable();
    renderInputParamsDiv();
}

function loadPeptidesFile() {
    peptides = [];
    var inputFile = document.getElementById('peptide-sequences-file');
    if (inputFile.files.length > 0) {
        clearTable();
        var file = inputFile.files[0];
        let extension = file.name.substring(file.name.lastIndexOf(".") + 1);
        if (extension.toLowerCase() == 'xlsx') {
            readXlsxFile(file).then(function(data) {
                var i = 0;
                data.map((row, index)=>{
                    if (i > 0) { // data
                        addPeptideRow(row[0], row[1]);
                    }
                    i++;
                });
            });
        } else if(extension.toLowerCase() == 'csv') {
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
                        addPeptideRow(columns[0], columns[1]);
                    }
                }
            };
            fileReader.readAsText(file);
        }
    }
}

function resetPeptidesFile() {
    document.getElementById('peptide-sequences-file').value = '';
}

function renderProcessingDiv() {
    var div1 = document.getElementById('generate-input');
    var div2 = document.getElementById('generate-processing');
    var div3 = document.getElementById('generate-output');
    div1.style.display = "none";
    div2.style.display = "block";
    div3.style.display = "none";
    document.getElementById('processing-num-peptides').innerHTML = numPeptides;
    document.getElementById('processing-num-peptides-per-pool').innerHTML = numPeptidesPerPool;
    document.getElementById('processing-num-coverage').innerHTML = numCoverage;
    document.getElementById('processing-plate-size').innerHTML = plateSize;
    document.getElementById('processing-cluster-peptides').innerHTML = clusterPeptides;
    document.getElementById('processing-sequence-similarity-fxn').innerHTML = sequenceSimilarityFxn;
    document.getElementById('processing-sequence-similarity-threshold').innerHTML = sequenceSimilarityThreshold;
    document.getElementById('processing-max-iters').innerHTML = maxIterations;
    document.getElementById('processing-init-strategy').innerHTML = initializationStrategy;
    document.getElementById('processing-allow-extra').innerHTML = allowExtraPools;
    document.getElementById('processing-num-total-pools').innerHTML = numTotalPools;
    document.getElementById('processing-num-total-plates').innerHTML = numTotalPlates;
}

function renderOutputDiv(output) {
    let configuration = output[0];
    let configurationBenchReady = output[1];
    let metadata = output[3];
    let peptidePairs = output[4];

    // Save configuration
    const peptideIdsMap = new Map(Object.entries(configuration['peptide_id']));
    const peptideSequencesMap = new Map(Object.entries(configuration['peptide_sequence']));
    const coverageIdsMap = new Map(Object.entries(configuration['coverage_id']));
    const poolIdsMap = new Map(Object.entries(configuration['pool_id']));
    const plateIdsMap = new Map(Object.entries(configuration['plate_id']));
    const wellIdsMap = new Map(Object.entries(configuration['well_id']));
    assignments = [];
    for (const [key, value] of peptideIdsMap.entries()) {
        let currPeptideId = value;
        const currPeptideSequence = peptideSequencesMap.get(key);
        const currCoverageId = coverageIdsMap.get(key);
        const currPoolId = poolIdsMap.get(key);
        const currPlateId = plateIdsMap.get(key);
        const currWellId = wellIdsMap.get(key);
        var assignment = new Assignment(
            currPeptideId,
            currPeptideSequence,
            currPoolId,
            currCoverageId,
            currPlateId,
            currWellId
        )
        assignments.push(assignment);
    }

    // Save configuration (bench ready)
    const benchReadyPlateIdsMap = new Map(Object.entries(configurationBenchReady['plate_id']));
    const benchReadyWellIdsMap = new Map(Object.entries(configurationBenchReady['well_id']));
    const benchReadyPeptideIdsMap = new Map(Object.entries(configurationBenchReady['peptide_ids']));
    const benchReadyPeptideSequencesMap = new Map(Object.entries(configurationBenchReady['peptide_sequences']));
    assignmentsBenchReady = [];
    for (const [key, value] of benchReadyPlateIdsMap.entries()) {
        let currPlateId = value;
        const currWellId = benchReadyWellIdsMap.get(key);
        const currPeptideIds = benchReadyPeptideIdsMap.get(key);
        const currPeptideSequences = benchReadyPeptideSequencesMap.get(key);
        var assignmentBenchReady = new AssignmentBenchReady(
            currPlateId,
            currWellId,
            currPeptideIds,
            currPeptideSequences
        )
        assignmentsBenchReady.push(assignmentBenchReady);
    }

    // Save preferred peptide pairs
    const pairsPeptide1IdsMap = new Map(Object.entries(peptidePairs['peptide_1_id']));
    const pairsPeptide2IdsMap = new Map(Object.entries(peptidePairs['peptide_2_id']));
    const pairsPeptide1SequencesMap = new Map(Object.entries(peptidePairs['peptide_1_sequence']));
    const pairsPeptide2SequencesMap = new Map(Object.entries(peptidePairs['peptide_2_sequence']));
    const pairsSimilarityScoresMap = new Map(Object.entries(peptidePairs['similarity_score']));
    preferredPeptidePairs = [];
    for (const [key, value] of pairsPeptide1IdsMap.entries()) {
        let currPeptideId1 = value;
        const currPeptideId2 = pairsPeptide2IdsMap.get(key);
        const currPeptideSequence1 = pairsPeptide1SequencesMap.get(key);
        const currPeptideSequence2 = pairsPeptide2SequencesMap.get(key);
        const currSimilarityScore = pairsSimilarityScoresMap.get(key);
        var peptidePair = new PreferredPeptidePair(
            currPeptideId1,
            currPeptideSequence1,
            currPeptideId2,
            currPeptideSequence2,
            currSimilarityScore
        )
        preferredPeptidePairs.push(peptidePair);
    }

    // Update UI
    var div1 = document.getElementById('generate-input');
    var div2 = document.getElementById('generate-processing');
    var div3 = document.getElementById('generate-output');
    div1.style.display = "none";
    div2.style.display = "none";
    div3.style.display = "block";
    document.getElementById('output-num-peptides').innerHTML = numPeptides;
    document.getElementById('output-num-peptides-per-pool').innerHTML = numPeptidesPerPool;
    document.getElementById('output-num-coverage').innerHTML = numCoverage;
    document.getElementById('output-plate-size').innerHTML = plateSize;
    document.getElementById('output-cluster-peptides').innerHTML = clusterPeptides;
    document.getElementById('output-sequence-similarity-fxn').innerHTML = sequenceSimilarityFxn;
    document.getElementById('output-sequence-similarity-threshold').innerHTML = sequenceSimilarityThreshold;
    document.getElementById('output-max-iters').innerHTML = maxIterations;
    document.getElementById('output-init-strategy').innerHTML = initializationStrategy;
    document.getElementById('output-allow-extra').innerHTML = allowExtraPools;
    document.getElementById('output-num-total-pools').innerHTML = numTotalPools;
    document.getElementById('output-num-total-plates').innerHTML = numTotalPlates;

    // Clear assignments table
    var table = document.getElementById('assignment-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }

    // Add assignments to table
    for (var i = 0; i < assignments.length; i++) {
        rowCount = table.rows.length;
        var row = table.insertRow(rowCount);

        // Peptide ID
        var newCell	= row.insertCell(0);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + assignments[i].peptide_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Peptide sequence
        var newCell	= row.insertCell(1);
        newCell.innerHTML = "<p class=\"text-center my-auto\" style=\"font-family:monospace;\">" + assignments[i].peptide_sequence + "</p>";
        newCell.style = "vertical-align: middle;";

        // Coverage ID
        var newCell	= row.insertCell(2);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + assignments[i].coverage_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Plate ID
        var newCell	= row.insertCell(3);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + assignments[i].plate_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Well ID
        var newCell	= row.insertCell(4);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + assignments[i].well_id + "</p>";
        newCell.style = "vertical-align: middle;";
    }

    // Clear peptide pairs table
    var table = document.getElementById('peptide-pairs-table');
    var rowCount = table.rows.length;
    while (rowCount > 1) {
        table.deleteRow(rowCount - 1);
        rowCount = table.rows.length;
    }

    // Add preferred peptides to table
    for (var i = 0; i < preferredPeptidePairs.length; i++) {
        rowCount = table.rows.length;
        var row = table.insertRow(rowCount);

        // Peptide ID 1
        var newCell	= row.insertCell(0);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + preferredPeptidePairs[i].peptide_1_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Peptide sequence 1
        var newCell	= row.insertCell(1);
        newCell.innerHTML = "<p class=\"text-center my-auto\" style=\"font-family:monospace;\">" + preferredPeptidePairs[i].peptide_1_sequence + "</p>";
        newCell.style = "vertical-align: middle;";

        // Peptide ID 2
        var newCell	= row.insertCell(2);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + preferredPeptidePairs[i].peptide_2_id + "</p>";
        newCell.style = "vertical-align: middle;";

        // Peptide sequence 2
        var newCell	= row.insertCell(3);
        newCell.innerHTML = "<p class=\"text-center my-auto\" style=\"font-family:monospace;\">" + preferredPeptidePairs[i].peptide_2_sequence + "</p>";
        newCell.style = "vertical-align: middle;";

        // Similarity score
        var newCell	= row.insertCell(4);
        newCell.innerHTML = "<p class=\"text-center my-auto\">" + preferredPeptidePairs[i].similarity_score.toFixed(3) + "</p>";
        newCell.style = "vertical-align: middle;";
    }
}

function generate() {
    if (isReady()) {
        if (expandedAccordion == 'sequences') {
            sequencesAvailable = true;
            numPeptides = peptides.length;
        } else {
            peptides = [];
        }

        clusterPeptides = document.getElementById('cluster-peptides').value;
        sequenceSimilarityFxn = document.getElementById('sequence-similarity-fxn').value;
        sequenceSimilarityThreshold = document.getElementById('sequence-similarity-threshold').value;
        maxIterations = document.getElementById('num-iterations').value;
        initializationStrategy = document.getElementById('strategy').value;
        allowExtraPools = document.getElementById('allow-extra-pools').value;
        renderProcessingDiv();
        eel.generate_configuration(
            sequencesAvailable,
            peptides,
            numPeptides,
            numPeptidesPerPool,
            numCoverage,
            plateSize,
            clusterPeptides,
            sequenceSimilarityFxn,
            sequenceSimilarityThreshold,
            maxIterations,
            initializationStrategy,
            allowExtraPools
        )(renderOutputDiv)
    }
}

function s2ab(s) {
  var buf = new ArrayBuffer(s.length);
  var view = new Uint8Array(buf);
  for (var i=0; i!=s.length; ++i) {
      view[i] = s.charCodeAt(i) & 0xFF;
  }
  return buf;
}

async function saveConfigurationFile() {
    if ( window.showSaveFilePicker ) {
        var zip = new JSZip();

        // Step 1. Create ELISpot configuration CSV file
        var csvRows1 = ["peptide_id,peptide_sequence,plate_id,well_id"];
        for (var i = 0; i < assignments.length; i++) {
            var values = [];
            values.push(assignments[i].peptide_id);
            values.push(assignments[i].peptide_sequence);
            values.push(assignments[i].plate_id);
            values.push(assignments[i].well_id);
            csvRows1.push(values.join(','));
        }
        let csvFile1Content = csvRows1.join('\n');
        zip.file("ace_elispot_assignment.csv", csvFile1Content);

        // Step 2. Create ELISpot configuration (bench-ready) CSV file
        var csvRows2 = ["plate_id,well_id,peptide_ids,peptide_sequences"];
        for (var i = 0; i < assignmentsBenchReady.length; i++) {
            var values = [];
            values.push(assignmentsBenchReady[i].plate_id);
            values.push(assignmentsBenchReady[i].well_id);
            values.push(assignmentsBenchReady[i].peptide_ids);
            values.push(assignmentsBenchReady[i].peptide_sequences);
            csvRows2.push(values.join(','));
        }
        let csvFile2Content = csvRows2.join('\n');
        zip.file("ace_elispot_assignment_bench_ready.csv", csvFile2Content);

        // Step 3. Create ELISpot configuration parameters CSV file
        var csvHeader = [
            'num_peptides',
            'num_peptides_per_pool',
            'num_coverage',
            'plate_size',
            'cluster_peptides',
            'sequence_similarity_function',
            'sequence_similarity_threshold',
            'maximum_iterations',
            'initialization_strategy',
            'allow_extra_pools',
            'num_total_pools',
            'num_total_plates'
        ]
        var csvRowContent = [
            numPeptides,
            numPeptidesPerPool,
            numCoverage,
            plateSize,
            clusterPeptides,
            sequenceSimilarityFxn,
            sequenceSimilarityThreshold,
            maxIterations,
            initializationStrategy,
            allowExtraPools,
            numTotalPools,
            numTotalPlates
        ]
        var csvRows3 = [
            csvHeader.join(','),
            csvRowContent.join(',')
        ]
        let csvFile3Content = csvRows3.join('\n');
        zip.file("ace_elispot_parameters.csv", csvFile3Content);

        // Step 4. Create ELISpot preferred peptide pairs CSV file
        var csvRows4 = ["peptide_1_id,peptide_1_sequence,peptide_2_id,peptide_2_sequence,similarity_score"];
        for (var i = 0; i < preferredPeptidePairs.length; i++) {
            var values = [];
            values.push(preferredPeptidePairs[i].peptide_1_id);
            values.push(preferredPeptidePairs[i].peptide_1_sequence);
            values.push(preferredPeptidePairs[i].peptide_2_id);
            values.push(preferredPeptidePairs[i].peptide_2_sequence);
            values.push(preferredPeptidePairs[i].similarity_score);
            csvRows4.push(values.join(','));
        }
        let csvFile4Content = csvRows4.join('\n');
        zip.file("ace_elispot_preferred_peptide_pairs.csv", csvFile4Content);

        // Step 5. Create ELISpot peptides CSV file
        var csvRows5 = ["peptide_id,peptide_sequence"];
        for (var i = 0; i < peptides.length; i++) {
            var values = [];
            values.push(peptides[i].id);
            values.push(peptides[i].sequence);
            csvRows5.push(values.join(','));
        }
        let csvFile5Content = csvRows5.join('\n');
        zip.file("ace_elispot_peptides.csv", csvFile5Content);

        // Step 6. Create ELISpot configuration XLSX file
        var wb = XLSX.utils.book_new();

        // Sheet 1 - Assignment
        var sheetData1 = [
            ['peptide_id',
             'peptide_sequence',
             'plate_id',
             'well_id']
        ];
        for (var i = 0; i < assignments.length; i++) {
            sheetData1.push([
                assignments[i].peptide_id,
                assignments[i].peptide_sequence,
                assignments[i].plate_id,
                assignments[i].well_id
            ]);
        }

        // Sheet 2 - Assignment (bench-ready)
        var sheetData2 = [
            ['plate_id',
             'well_id',
             'peptide_ids',
             'peptide_sequences']
        ];
        for (var i = 0; i < assignmentsBenchReady.length; i++) {
            sheetData2.push([
                assignmentsBenchReady[i].plate_id,
                assignmentsBenchReady[i].well_id,
                assignmentsBenchReady[i].peptide_ids,
                assignmentsBenchReady[i].peptide_sequences
            ]);
        }

        // Sheet 3 - Peptides
        var sheetData3 = [
            ['peptide_id', 'peptide_sequence']
        ]
        for (var i = 0; i < peptides.length; i++) {
            sheetData3.push([
                peptides[i].id,
                peptides[i].sequence
            ])
        }

        // Sheet 4 - Preferred Peptide Pairs
        var sheetData4 = [
            ['peptide_1_id', 'peptide_1_sequence', 'peptide_2_id', 'peptide_2_sequence', 'similarity_score']
        ]
        for (var i = 0; i < preferredPeptidePairs.length; i++) {
            sheetData4.push([
                preferredPeptidePairs[i].peptide_1_id,
                preferredPeptidePairs[i].peptide_1_sequence,
                preferredPeptidePairs[i].peptide_2_id,
                preferredPeptidePairs[i].peptide_2_sequence,
                preferredPeptidePairs[i].similarity_score
            ])
        }

        // Sheet 5 - Parameters
        var sheetData5 = [
            ['num_peptides',
             'num_peptides_per_pool',
             'num_coverage',
             'plate_size',
             'cluster_peptides',
             'sequence_similarity_function',
             'sequence_similarity_threshold',
             'maximum_iterations',
             'initialization_strategy',
             'allow_extra_pools',
             'num_total_pools',
             'num_total_plates'
            ],
            [
                numPeptides,
                numPeptidesPerPool,
                numCoverage,
                plateSize,
                clusterPeptides,
                sequenceSimilarityFxn,
                sequenceSimilarityThreshold,
                maxIterations,
                initializationStrategy,
                allowExtraPools,
                numTotalPools,
                numTotalPlates
            ]
        ]

        // Create worksheet instances
        const ws1 = XLSX.utils.aoa_to_sheet(sheetData1);
        const ws2 = XLSX.utils.aoa_to_sheet(sheetData2);
        const ws3 = XLSX.utils.aoa_to_sheet(sheetData3);
        const ws4 = XLSX.utils.aoa_to_sheet(sheetData4);
        const ws5 = XLSX.utils.aoa_to_sheet(sheetData5);

        // Add worksheets to the workbook
        XLSX.utils.book_append_sheet(wb, ws1, 'assignment');
        XLSX.utils.book_append_sheet(wb, ws2, 'assignment_bench_ready');
        XLSX.utils.book_append_sheet(wb, ws3, 'peptides');
        XLSX.utils.book_append_sheet(wb, ws4, 'preferred_peptide_pairs');
        XLSX.utils.book_append_sheet(wb, ws5, 'parameters');

        // Generate a blob from the workbook
        const wbArrayBuffer = XLSX.write(wb, { bookType: 'xlsx', type: 'array' });

        // Convert the ArrayBuffer to a Blob
        const wbBlob = new Blob([wbArrayBuffer], {
            type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        });

        zip.file("ace_elispot_configuration.xlsx", wbBlob);

        var zipContent = await zip.generateAsync({ type: "blob" });
        const opts = {
            types: [{
                description: 'ACE_ELISpot_Configuration',
                accept: {'text/plain': ['.zip']},
            }],
            suggestedName: 'ace_elispot_configuration',
        };
        var handle = await showSaveFilePicker(opts);
        var writable = await handle.createWritable();
        await writable.write(zipContent);
        writable.close();
    }
}

async function downloadExampleFile() {
    if ( window.showSaveFilePicker ) {
        const opts = {
            types: [{
                description: 'ACE_ELISpot_Configuration'
            }],
            suggestedName: 'ace_example_peptides.xlsx',
        };
        var handle = await showSaveFilePicker(opts);
        var writable = await handle.createWritable();
        const downloadLink = document.getElementById('download-example-file');
        const fileUrl = downloadLink.getAttribute('href');
        const response = await fetch(fileUrl);
        const fileBlob = await response.blob();
        await writable.write( fileBlob );
        writable.close();
    }
}

function isValidAminoAcidSequence(sequence) {
  const validAminoAcidPattern = /^[ACDEFGHIKLMNPQRSTVWY]+$/;
  return validAminoAcidPattern.test(sequence);
}

function isReady() {
    if (expandedAccordion == 'sequences') {
        for (var i = 0; i < peptides.length; i++) {
            if (peptides[i]['id'] == '') {
                document.getElementById('alert-message').innerHTML = 'You have selected to specify peptide sequences. Please make sure all peptides have IDs.'
                new bootstrap.Modal(document.querySelector("#alert-modal")).show();
                return false;
            }
            if (peptides[i]['sequence'] == '') {
                document.getElementById('alert-message').innerHTML = 'You have selected to specify peptide sequences. Please make sure all peptides have sequences.'
                new bootstrap.Modal(document.querySelector("#alert-modal")).show();
                return false;
            }
            if (isValidAminoAcidSequence(peptides[i]['sequence']) == false) {
                document.getElementById('alert-message').innerHTML = peptides[i]['sequence'] + ' is not a valid amino acid sequence.';
                new bootstrap.Modal(document.querySelector("#alert-modal")).show();
                return false;
            }
        }
        var peptideIds = [];
        for (var i = 0; i < peptides.length; i++) {
            peptideIds.push(peptides[i]['id'])
        }
        const peptideIdsSet = new Set(peptideIds);
        if (peptideIdsSet.size < peptides.length) {
            document.getElementById('alert-message').innerHTML = 'You have selected to specify peptide sequences. Please make sure all peptides have unique IDs.'
            new bootstrap.Modal(document.querySelector("#alert-modal")).show();
            return false;
        }
    } else {
        if (numPeptides == 0) {
            document.getElementById('alert-message').innerHTML = 'You have selected to specify a peptide count. Please make sure you specify a non-zero peptide count.'
            new bootstrap.Modal(document.querySelector("#alert-modal")).show();
            return false;
        }
    }
    if (numPeptidesPerPool == 0) {
        document.getElementById('alert-message').innerHTML = 'Please specify a non-zero number of peptides per pool.'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    if (numCoverage == 0) {
        document.getElementById('alert-message').innerHTML = 'Please specify a non-zero number of peptide replicates (coverage).'
        new bootstrap.Modal(document.querySelector("#alert-modal")).show();
        return false;
    }
    return true;
}

