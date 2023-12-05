var elispotConfiguration = '';
var peptideIdToWellIds = '';    // key = peptide ID, value = set of [plate number, plate well ID]
var wellIdToPeptideIds = '';    // key = [plate number, plate well ID], value = set of peptide IDs
var peptideIdsSequences = '';   // array of [peptide ID, sequence, peptide ID formatted pretty]
var selectedPlateNumber = 1;
var previouslySelectedPeptideId = '';
var previouslySelectedWellId = '';
var selectedWellIds = '';

window.onload = function() {
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    localStorage.removeItem("elispot-configuration");
    selectedWellIds = new Set();
    loadConfigurationData();
    loadPeptideSequencesList();
    loadWells(1);
};

async function exportConfiguration() {
	var excel = await eel.download_xls_bttn()();
        let excelFilePath = String(excel);
        if (excelFilePath.length > 0) {
            let returnValue = eel.write_elispot_configuration_excel(
                elispotConfiguration,
                excelFilePath
            )(onfinishExport)
        }
}

function onfinishExport(returnValue) {
    if (returnValue) {
        window.alert("Successfully exported ELIspot configuration file.");
    }
}

function onclickPeptideSequenceListItem(peptideId) {
    // Step 1. Unselect previously selected peptide and wells
    if (previouslySelectedPeptideId) {
        // Unselect peptide
        document.getElementById(previouslySelectedPeptideId).className = "text-peptide-sequence-list-item";

        // Unselect wells
        let wellIdsArr = Array.from(peptideIdToWellIds[previouslySelectedPeptideId]);
        for (let i = 0; i < wellIdsArr.length; i++) {
            let currPlateNumber = wellIdsArr[i][0];
            let currPlateWellId = wellIdsArr[i][1];
            if (currPlateNumber == selectedPlateNumber) {
                document.getElementById(currPlateWellId).className = "btn-plate-assigned";
            }
        }
    }

    // Step 2. Update previously selected peptide id
    previouslySelectedPeptideId = peptideId;

    // Step 3. Update the class of current peptide id
    document.getElementById(peptideId).className = "text-peptide-sequence-list-item-selected";

    // Step 4. Select all corresponding wells
    let plateUniqueIdsArr = Array.from(peptideIdToWellIds[peptideId]);
    selectedWellIds = new Set();
    var plate1Count = 0;
    var plate2Count = 0;
    var plate3Count = 0;
    var plate4Count = 0;
    var plate5Count = 0;
    for (let i = 0; i < plateUniqueIdsArr.length; i++) {
        let currPlateNumber = plateUniqueIdsArr[i][0];
        let currPlateWellId = plateUniqueIdsArr[i][1];
        if (currPlateNumber == selectedPlateNumber) {
            document.getElementById(currPlateWellId).className = "btn-plate-peptide-pool";
            selectedWellIds.add(currPlateWellId);
        }
        if (currPlateNumber == 1) {
            plate1Count++;
        }
        if (currPlateNumber == 2) {
            plate2Count++;
        }
        if (currPlateNumber == 3) {
            plate3Count++;
        }
        if (currPlateNumber == 4) {
            plate4Count++;
        }
        if (currPlateNumber == 5) {
            plate5Count++;
        }
    }
    if (plate1Count > 0) {
        document.getElementById("panel-1-pools-count").className = "text-panel-pools-count-enabled";
        document.getElementById("panel-1-pools-count").innerHTML = plate1Count + " pools";
    } else {
        document.getElementById("panel-1-pools-count").className = "text-panel-pools-count-disabled";
    }
    if (plate2Count > 0) {
        document.getElementById("panel-2-pools-count").className = "text-panel-pools-count-enabled";
        document.getElementById("panel-2-pools-count").innerHTML = plate2Count + " pools";
    } else {
        document.getElementById("panel-2-pools-count").className = "text-panel-pools-count-disabled";
    }
    if (plate3Count > 0) {
        document.getElementById("panel-3-pools-count").className = "text-panel-pools-count-enabled";
        document.getElementById("panel-3-pools-count").innerHTML = plate3Count + " pools";
    } else {
        document.getElementById("panel-3-pools-count").className = "text-panel-pools-count-disabled";
    }
    if (plate4Count > 0) {
        document.getElementById("panel-4-pools-count").className = "text-panel-pools-count-enabled";
        document.getElementById("panel-4-pools-count").innerHTML = plate4Count + " pools";
    } else {
      document.getElementById("panel-4-pools-count").className = "text-panel-pools-count-disabled";
    }
    if (plate5Count > 0) {
        document.getElementById("panel-5-pools-count").className = "text-panel-pools-count-enabled";
        document.getElementById("panel-5-pools-count").innerHTML = plate5Count + " pools";
    } else {
      document.getElementById("panel-5-pools-count").className = "text-panel-pools-count-disabled";
    }

}

function onclickWell(wellId) {
    // Step 1. Unselect previously selected well
    if (previouslySelectedWellId) {
        if (selectedWellIds.has(previouslySelectedWellId)) {
            document.getElementById(previouslySelectedWellId).className = "btn-plate-peptide-pool";
        } else {
            document.getElementById(previouslySelectedWellId).className = "btn-plate-assigned";
        }
    }

    // Step 2. Update previously selected well id
    previouslySelectedWellId = wellId;

    // Step 3. Select well
    if (selectedWellIds.has(wellId)) {
        document.getElementById(wellId).className = "btn-plate-selected-peptide-pool";
    } else {
        document.getElementById(wellId).className = "btn-plate-selected-pool";
    }

    // Step 4. Update selected pool features
    document.getElementById("container-selected-pool-features").className = "container-selected-pool-active";
    var htmlValue = "Peptides in this pool:";
    let peptideIds = Array.from(wellIdToPeptideIds[selectedPlateNumber + "," + wellId]);
    for (let i = 0; i < peptideIdsSequences.length; i++) {
        if (peptideIds.includes(peptideIdsSequences[i][0])) {
            htmlValue = htmlValue + "<br/>" + peptideIdsSequences[i][2];
        }
    }
    document.getElementById("selected-pool-features").innerHTML = htmlValue;

    // Step 5. Update selected pool id
    document.getElementById("selected-pool-id").innerHTML = "Selected Pool: " + wellId;
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
