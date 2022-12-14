var elispotConfiguration = '';

window.onload = function() {
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    console.log(elispotConfiguration);
    updatePeptideSequencesList();
    renderWells();
};

function onclickPeptideSequenceListItem(peptideId) {
    // Step 1. Figure out all 
    console.log(peptideId);
}

function onclickWell(wellId) {
    console.log(wellId);
}

function updatePeptideSequencesList() {
    // Step 1. Clear the existing list
    document.getElementById("peptide-list").innerHTML = "";

    // Step 2. Fetch the peptide IDs in order
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
    var peptideIdSequences = [];
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
        peptideIdSequences.push([currPeptideId, currPeptideIdPretty + ' (' + currPeptideSequence + ')']);
    }

    // Step 4. Update the list
    var list = document.getElementById('peptide-list');
    for (let i = 0; i < peptideIdSequences.length; i++) {
        var entry = document.createElement('li');
        entry.setAttribute("id", peptideIdSequences[i][0]);
        entry.classList.add('text-peptide-sequence-list-item');
        entry.innerHTML = peptideIdSequences[i][1];
        entry.setAttribute("onclick", "onclickPeptideSequenceListItem(this.id)");
        list.appendChild(entry);
    }
}

function renderWells() {
    // Step 1. Figure out the well IDs in use
    var wellIdsSet = new Set();
    for (const [key, value] of Object.entries(elispotConfiguration['plate_well_id'])) {
        wellIdsSet.add(value);
    }
    // Step 2. Change the class of all wells in use
    let wellIdsArr = Array.from(wellIdsSet);
    for (let i = 0; i < wellIdsArr.length; i++) {
        document.getElementById(wellIdsArr[i]).className = "btn-plate-assigned";
        document.getElementById(wellIdsArr[i]).setAttribute("onclick", "onclickWell(this.id)");
    }
}
