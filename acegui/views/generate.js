function generateConfiguration() {
    var numPeptides = document.getElementById("number-of-peptides").value
    var numPeptidesPerPool = document.getElementById("number-of-peptides-per-pool").value
    var numCoverage = document.getElementById("number-of-coverage").value
    var numCores = document.getElementById("number-of-cores").value
    let configuration = eel.generate_configuration(
        numPeptides,
        numPeptidesPerPool,
        numCoverage,
        numCores
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
		if (csv) {
			document.getElementById("input_CSV").value = csv;
		}
		return csv;
	}