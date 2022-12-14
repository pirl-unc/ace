var peptideSequences = [];
var isGenerationReady = false;

function generateConfiguration() {
    document.getElementById("container-parameters").style.display = "none";
    document.getElementById("container-running").style.display = "block";

    let numPeptides = document.getElementById("number-of-peptides").value
    let numPeptidesPerPool = document.getElementById("number-of-peptides-per-pool").value
    let numCoverage = document.getElementById("number-of-coverage").value
    let numCores = document.getElementById("number-of-cores").value
    let configuration = eel.generate_configuration(
        numPeptides,
        numPeptidesPerPool,
        numCoverage,
        numCores,
        peptideSequences
    )(renderConfiguration)
}

function renderConfiguration(configuration) {
    localStorage["elispot-configuration"] = JSON.stringify(configuration);
    window.location.href="configuration.html";
}

function showSelectedFilePath(csvFilePath) {
    document.getElementById("div-file-select").style.display = "none";
    document.getElementById("div-file-selected").style.display = "flex";
    document.getElementById("peptide-sequence-file-path").innerHTML = csvFilePath;
}

function clearSelectedFilePath() {
    document.getElementById("peptide-sequence-file-path").innerHTML = "";
    document.getElementById("div-file-selected").style.display = "none";
    document.getElementById("div-file-select").style.display = "block";
    document.getElementById('number-of-peptides').readOnly = false;
    document.getElementById("number-of-peptides").value = '';
    document.getElementById("enforce-peptide-constraints").checked = false;
}

function readPeptideSequencesFile(csvFilePath) {
    sequences = eel.read_peptide_sequences_csv_file(csvFilePath)(renderPeptidesCount);
}

function renderPeptidesCount(sequences) {
    peptideSequences = sequences;
    document.getElementById("number-of-peptides").value = peptideSequences.length;
    document.getElementById('number-of-peptides').readOnly = true;
    document.getElementById("enforce-peptide-constraints").checked = true;
}

function checkInputValidity() {
    let numPeptides = document.getElementById("number-of-peptides").value;
    let numPeptidesPerPool = document.getElementById("number-of-peptides-per-pool").value;
    let numRepeats = document.getElementById("number-of-coverage").value;
    if (numPeptides > 0 && numPeptidesPerPool > 0 && numRepeats > 0) {
        isGenerationReady = true;
        document.getElementById("btn-generate").className = "btn-primary-round-active";
    } else {
        isGenerationReady = false;
        document.getElementById("btn-generate").className = "btn-primary-round-inactive";
        window.alert("This is not a valid set of Conifguration Parameters please clear and try again.")
    }
}

async function getCSV() {
	var csv = await eel.upload_csv_bttn()();
        let csvFilePath = String(csv);
        if (csvFilePath.length > 0) {
            showSelectedFilePath(csvFilePath);
            readPeptideSequencesFile(csvFilePath);
        }
	}

