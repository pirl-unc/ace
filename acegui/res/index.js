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

function clearDisallowedPeptidePairs() {
    document.getElementById("disallowed-peptide-pairs").innerHTML = ""
}

function auto_grow(element) {
    element.style.height = "10px";
    element.style.height = (element.scrollHeight)+"px";
}
