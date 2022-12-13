function generateConfiguration() {
    var path_to_csv = document.getElementById("path-to-csv").value
    var numPeptides = document.getElementById("number-of-peptides").value
    var numPeptidesPerPool = document.getElementById("number-of-peptides-per-pool").value
    var numCoverage = document.getElementById("number-of-coverage").value
    var numCores = document.getElementById("number-of-cores").value
    let configuration = eel.generate_configuration(
        numPeptides,
        numPeptidesPerPool,
        numCoverage,
        numCores,
        path_to_csv
    )(renderConfiguration)
}

function renderConfiguration(configuration) {
    console.log(configuration);
}

function clearAssayParameters() {
    document.getElementById("number-of-peptides").innerHTML = ""
    document.getElementById("number-of-peptides-per-pool").innerHTML = ""
    document.getElementById("number-of-coverage").innerHTML = ""
}

async function getCSV() {
	var csv = await eel.upload_csv_bttn()();
        let csvFilePath = String(csv);
        if (csvFilePath.length > 0) {
            document.getElementById("peptide-sequence-file-path").innerHTML = csvFilePath;
        }
	}