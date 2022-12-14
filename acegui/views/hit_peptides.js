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

window.onload = function() {
    hitConfigs = JSON.parse(localStorage["hit-peptides"]);
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    plateReadout = JSON.parse(localStorage['plate-readout']);
    let [hitPeptides, hitPools] =  loadHitConfig(hitConfigs);
    loadConfigurationData(hitPools);
    loadWells(1);
    renderHitPeptides(hitConfigs);
    console.log(plateReadout);
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

}

function getHitPeptides(hitPools) {
    hitPeptideIdsSequences = [];
    for (const concatIdSeq of peptideIdsSequences){
        if (Object.keys(hitPools).includes(concatIdSeq[0])){
            hitPeptideIdsSequences.push(concatIdSeq);
        }
    }
    return hitPeptideIdsSequences
}

async function loadHitPeptides(hitPeptideIdsSequences) {
    // Step 1. Clear the existing list
    document.getElementById("hit-peptide-list").innerHTML = "";
    // Step 2. Update the list
    var list = document.getElementById('hit-peptide-list');
    for (let i = 0; i < Object.keys(hitPeptideIdsSequences).length; i++) {
        var entry = document.createElement('li');
        entry.setAttribute("id", hitPeptideIdsSequences[i][0]);
        entry.classList.add('text-hit-peptide-list-item');
        entry.innerHTML = hitPeptideIdsSequences[i][2];
        entry.setAttribute("onclick", "onclickPeptideSequenceListItem(this.id)");
        list.appendChild(entry);
    }
}

function getWell(pool_id) {
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

async function loadHitPools() {
    // Change the class of all wells in use
    let result_pools = new Set();
    for (const set of Object.values(hitPools))
        for (const element of set) {
            result_pools.add(element);
        }

    for (const poolId of result_pools) {
        document.getElementById(getWell(poolId)).className = "btn-plate-hit-pool";
        document.getElementById(getWell(poolId)).setAttribute("onclick", "onclickWell(this.id)");
    }
}

async function loadWells(plateNumber) {
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

async function updateHitRate(hitPeptideIdsSequences) {
    let numHits = hitPeptideIdsSequences.length;
    let numTotal = peptideIdsSequences.length;
    let hitRate = Math.round((numHits / numTotal) * 100 * 100) / 100;
    document.getElementById("hit-rate").innerHTML = "Hit rate: " + numHits + "/" + numTotal + " (" + hitRate + "%)"
}

function onupdateSlider() {
    let spotCount = document.getElementById("spot-slider").value
    document.getElementById("thresholdValue").innerHTML = spotCount
    let hitConfigs = eel.ace_identify_helper(
        plateReadout,
        elispotConfiguration,
        parseInt(spotCount)
    )(renderHitPeptides)
}

async function renderHitPeptides(hitConfigs) {
    let [hitPeptides, hitPools] =  loadHitConfig(hitConfigs);
    let hitPeptideIdsSequences = getHitPeptides(hitPools)
    loadHitPeptides(hitPeptideIdsSequences);
    loadWells(1);
    loadHitPools(hitPools);
    updateHitRate(hitPeptideIdsSequences);
}