var elispotConfiguration = '';
var peptideIdToWellIds = '';    // key = peptide ID, value = set of [plate number, plate well ID]
var wellIdToPeptideIds = '';    // key = [plate number, plate well ID], value = set of peptide IDs
var peptideIdsSequences = '';   // array of [peptide ID, sequence, peptide ID formatted pretty]
var selectedPlateNumber = 1;

window.onload = function() {
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    console.log(elispotConfiguration);
    loadConfigurationData();
    loadPeptideSequencesList();
    loadWells(1);
};

function onclickPeptideSequenceListItem(peptideId) {
    // Step 1. Update the class of current peptide id
    document.getElementById(peptideId).className = "text-peptide-sequence-list-item-selected";

    // Step 2. Select all corresponding wells
    let plateUniqueIdsArr = Array.from(peptideIdToWellIds[peptideId]);
    for (let i = 0; i < plateUniqueIdsArr.length; i++) {
        let currPlateNumber = plateUniqueIdsArr[i][0];
        let currPlateWellId = plateUniqueIdsArr[i][1];
        if (currPlateNumber == selectedPlateNumber) {
            document.getElementById(currPlateWellId).className = "btn-plate-peptide-pool";
        }
    }
}

function onclickWell(wellId) {
    console.log(wellId);
}

function loadConfigurationData() {
    // Step 1. Load peptideIdToWellIds
    // Initialize peptideIdToWellIds keys
    peptideIdToWellIds = Object.create(null);
    for (const [peptideIdIdx, peptideId] of Object.entries(elispotConfiguration['peptide_id'])) {
        peptideIdToWellIds[peptideId] = new Set();
    }
    // Add well IDs to peptideIdToWellIds values
    for (const [peptideIdIdx, peptideId] of Object.entries(elispotConfiguration['peptide_id'])) {
        for (const [plateUniqueIdIdx, plateUniqueId] of Object.entries(elispotConfiguration['plate_unique_id'])) {
            let plateNumber = Number(plateUniqueId.split('_')[0].replace('plate', ''));
            let plateWellId = plateUniqueId.split('_')[1];
            if (plateUniqueIdIdx == peptideIdIdx) {
                peptideIdToWellIds[peptideId].add([plateNumber, plateWellId]);
            }
        }
    }

    // Step 2. Load wellIdToPeptideIds
    // Initialize wellIdToPeptideIds keys
    wellIdToPeptideIds = Object.create(null);
    for (const [plateUniqueIdIdx, plateUniqueId] of Object.entries(elispotConfiguration['plate_unique_id'])) {
        let plateNumber = Number(plateUniqueId.split('_')[0].replace('plate', ''));
        let plateWellId = plateUniqueId.split('_')[1];
        wellIdToPeptideIds[[plateNumber, plateWellId]] = new Set();
    }
    // Add well IDs to wellIdToPeptideIds values
    for (const [plateUniqueIdIdx, plateUniqueId] of Object.entries(elispotConfiguration['plate_unique_id'])) {
        let plateNumber = Number(plateUniqueId.split('_')[0].replace('plate', ''));
        let plateWellId = plateUniqueId.split('_')[1];
        for (const [peptideIdIdx, peptideId] of Object.entries(elispotConfiguration['peptide_id'])) {
            if (peptideIdIdx == plateUniqueIdIdx) {
                wellIdToPeptideIds[[plateNumber, plateWellId]].add(peptideId);
            }
        }
    }

    // Step 3. Load peptideIdsSequences
    // Fetch the peptide IDs in order
    var peptideIdsSet = new Set();
    for (const [key, value] of Object.entries(elispotConfiguration['peptide_id'])) {
        peptideIdsSet.add(value);
    }
    var peptideIdsArr = Array.from(peptideIdsSet).sort(function(x, y) {
        const peptideNum1 = Number(x.split('_')[1]);
        const peptideNum2 = Number(y.split('_')[1]);
        if (peptideNum1 < peptideNum2) {
            return -1;
        }
        if (peptideNum1 > peptideNum2) {
            return 1;
        }
        return 0;
    });

    // Step 3. Fetch the peptide sequences
    peptideIdsSequences = [];
    for (let i = 0; i < peptideIdsArr.length; i++) {
        let currPeptideId = peptideIdsArr[i];
        var currPeptideSequence = '';
        for (const [key, value] of Object.entries(elispotConfiguration['peptide_id'])) {
            if (value == currPeptideId) {
                currPeptideSequence = elispotConfiguration['sequence'][key];
                break;
            }
        }
        let currPeptideIdPretty = 'Peptide ' + String(currPeptideId.split('_')[1]);
        peptideIdsSequences.push([currPeptideId, currPeptideSequence, currPeptideIdPretty + ' (' + currPeptideSequence + ')']);
    }
}

function loadPeptideSequencesList() {
    // Step 1. Clear the existing list
    document.getElementById("peptide-list").innerHTML = "";

    // Step 2. Update the list
    var list = document.getElementById('peptide-list');
    for (let i = 0; i < peptideIdsSequences.length; i++) {
        var entry = document.createElement('li');
        entry.setAttribute("id", peptideIdsSequences[i][0]);
        entry.classList.add('text-peptide-sequence-list-item');
        entry.innerHTML = peptideIdsSequences[i][2];
        entry.setAttribute("onclick", "onclickPeptideSequenceListItem(this.id)");
        list.appendChild(entry);
    }
}

function loadWells(plateNumber) {
    // Change the class of all wells in use
    for (const [key, value] of Object.entries(wellIdToPeptideIds)) {
        let currPlateNumber = key.split(',')[0];
        let currPlateWellId = key.split(',')[1];
        if (currPlateNumber == plateNumber) {
            document.getElementById(currPlateWellId).className = "btn-plate-assigned";
            document.getElementById(currPlateWellId).setAttribute("onclick", "onclickWell(this.id)");
        }
    }
}
