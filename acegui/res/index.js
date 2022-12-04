function generateConfiguration() {
    var numPeptides = document.getElementById("number-of-peptides").value
    var numPeptidesPerPool = document.getElementById("number-of-peptides-per-pool").value
    var numCoverage = document.getElementById("number-of-coverage").value
    var numCores = document.getElementById("number-of-cores").value
    var pathToSeqs = document.getElementById('input_CSV').value
    let configuration = eel.generate_configuration(
        numPeptides,
        numPeptidesPerPool,
        numCoverage,
        numCores,
        pathToSeqs
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

function clearDisallowedPeptidePairs() {
    document.getElementById("disallowed-peptide-pairs").innerHTML = ""
}

function auto_grow(element) {
    element.style.height = "10px";
    element.style.height = (element.scrollHeight)+"px";
}

async function getCSV() {
	var csv = await eel.upload_csv_bttn()();
		if (csv) {
			document.getElementById("input_CSV").value = csv;
		}
		return csv;
	}

