var hitConfigs = '';
var hitPeptides = '';
var hitPools = '';
var elispotConfiguration = '';
var peptideIdToWellIds = '';    // key = peptide ID, value = set of [plate number, plate well ID]
var wellIdToPeptideIds = '';    // key = [plate number, plate well ID], value = set of peptide IDs
var peptideIdsSequences = '';   // array of [peptide ID, sequence, peptide ID formatted pretty]
var selectedPlateNumber = 1;

window.onload = function() {
    hitConfigs = JSON.parse(localStorage["hit-peptides"]);
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    console.log(hitConfigs);
    console.log(elispotConfiguration);
    loadConfigurationData();
    loadHitConfig();
    loadHitPeptides();
};

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
function loadHitConfig() { // Load relevant info from JSON Streamed Pandas DF

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
    console.log(hitPools)
}

function loadHitPeptides() {
    // Step 1. Clear the existing list
    document.getElementById("hit-peptide-list").innerHTML = "";
    // Step 2. Update the list
    var list = document.getElementById('hit-peptide-list');
    for (let i = 0; i < Object.keys(hitPeptides).length; i++) {
        var entry = document.createElement('li');
        entry.setAttribute("id", hitPeptides[i]);
        entry.classList.add('text-peptide-sequence-list-item');
        entry.innerHTML = hitPeptides[i];
        //entry.setAttribute("onclick", "onclickPeptideSequenceListItem(this.id)");
        list.appendChild(entry);
    }
}