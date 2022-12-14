var hitConfigs = ''
var hitPeptides = '';           //
var hitPools = '';              //
var plateReadout = '';
var elispotConfiguration = '';  // Input 2
var peptideIdToWellIds = '';    // key = peptide ID, value = set of [plate number, plate well ID]
var wellIdToPeptideIds = '';    // key = [plate number, plate well ID], value = set of peptide IDs
var peptideIdsSequences = '';   // array of [peptide ID, sequence, peptide ID formatted pretty]
var hitPeptideIdsSequences = '';   // array of [peptide ID, sequence, peptide ID formatted pretty]
var selectedPlateNumber = 1;
var previouslySelectedPeptideId = '';
var previouslySelectedWellId = '';
var selectedWellIds = '';

window.onload = function() {
    hitConfigs = JSON.parse(localStorage["hit-peptides"]);
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    plateReadout = JSON.parse(localStorage['plate-readout'])
    renderHitPeptides(hitConfigs);
    for (var i = 0; i < localStorage.length; i++){
        console.log(localStorage.getItem(localStorage.key(i)));
    }
};

function loadHitConfig(hitConfigs) { // Load relevant info from JSON Streamed Pandas DF

    // Get peptides
    hitPeptides = Object.create(null); //Working List of Hit Peptides
    hitPools = Object.create(null);
    for (const [peptideIdIdx, peptideId] of Object.entries(hitConfigs['hit_peptide_id'])) {
        hitPeptides[peptideIdIdx] = peptideId;
        hitPools[peptideId] = new Set()
    }

    // Get Pools

    for (const [poolIdIdx, poolIds] of Object.entries(hitConfigs['pool_ids'])){
        for (const poolId of poolIds.split(',')){
            hitPools[hitPeptides[poolIdIdx]].add(poolId)
        }
    }
    return [hitPeptides, hitPools]
}

function loadConfigurationData(hitPools) {
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

    hitPeptideIdsSequences = [];
    for (const concatIdSeq of peptideIdsSequences){
        if (Object.keys(hitPools).includes(concatIdSeq[0])){
            hitPeptideIdsSequences.push(concatIdSeq);
        }
    }

}

function loadHitPeptides() {
    // Step 1. Clear the existing list
    document.getElementById("hit-peptide-list").innerHTML = "";
    // Step 2. Update the list
    var list = document.getElementById('hit-peptide-list');
    for (let i = 0; i < Object.keys(hitPeptideIdsSequences).length; i++) {
        var entry = document.createElement('li');
        entry.setAttribute("id", hitPeptideIdsSequences[i][0]);
        entry.classList.add('text-peptide-sequence-list-item');
        entry.innerHTML = hitPeptideIdsSequences[i][2];
        entry.setAttribute("onclick", "onclickPeptideSequenceListItem(this.id)");
        list.appendChild(entry);
    }
}

function get_well(pool_id){
    const alphabetUp = "ABCDEFGH".split("");
    let int_id = parseInt(pool_id.split("_")[1]);
    let plate = Math.floor(int_id / 96);
    let plate_spillover = int_id % 96;
    if (plate_spillover == 0){
        //Catch Last H12
        var row = 8
        var col = 12
    } else {
        plate = plate + 1;
        var row = Math.floor(plate_spillover / 12);
        let row_spillover = plate_spillover % 12;
        if (row_spillover == 0){

            var col = 12;
        } else {
            var col = row_spillover;
            row = row + 1;
        }
    }
    return alphabetUp[row-1] + col.toString()
}

function loadHitPools() {
    // Change the class of all wells in use
    let result_pools = new Set;
    for (const set of Object.values(hitPools))
        for (const element of set)
            result_pools.add(element);

    for (const poolId of result_pools) {
        document.getElementById(get_well(poolId)).className = "btn-plate-hit-pool";
        document.getElementById(get_well(poolId)).setAttribute("onclick", "onclickWell(this.id)");
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

function rerun_from_slider(){
    var spotCount = document.getElementById("spot-slider").value
    document.getElementById("thresholdValue").innerText = spotCount
    let hitConfigs = eel.ace_identify_helper(plateReadout, elispotConfiguration, parseInt(spotCount))(renderHitPeptides)
}

function renderHitPeptides(hitConfigs){
    let [hitPeptides, hitPools] =  loadHitConfig(hitConfigs);
    loadConfigurationData(hitPools);
    loadWells(1);
    loadHitPeptides();
    loadHitPools(hitPools);
}