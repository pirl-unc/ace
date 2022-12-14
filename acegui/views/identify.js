function identifyPositives() {
    var readoutPath = document.getElementById("counts-file-path").innerHTML
    var configPath = document.getElementById("config-file-path").innerHTML
    let joined_output = eel.identify_positives(
        readoutPath,
        configPath
    )(renderHitPeptides)
}

function showSelectedConfigurationFilePath(xlsFilePath) {
    document.getElementById("config-file-select").style.display = "none";
    document.getElementById("config-file-selected").style.display = "flex";
    document.getElementById("config-file-path").innerHTML = xlsFilePath;
}

function clearConfigurationFilePath() {
    document.getElementById("config-file-path").innerHTML = "";
    document.getElementById("config-file-selected").style.display = "none";
    document.getElementById("config-file-select").style.display = "block";
}

function showSelectedCountsFilePath(xlsFilePath) {
    document.getElementById("counts-file-select").style.display = "none";
    document.getElementById("counts-file-selected").style.display = "flex";
    document.getElementById("counts-file-path").innerHTML = xlsFilePath;
}

function clearCountsFilePath() {
    document.getElementById("counts-file-path").innerHTML = "";
    document.getElementById("counts-file-selected").style.display = "none";
    document.getElementById("counts-file-select").style.display = "block";
}

function renderHitPeptides(joined_output) {
    localStorage["hit-peptides"] = JSON.stringify(joined_output[0]);
    localStorage["elispot-configuration"] = JSON.stringify(joined_output[1]);
    localStorage["plate-readout"] = JSON.stringify(joined_output[2]);
    window.location.href="hit_peptides.html";
}

async function getXLS(type) {
	var xls = await eel.upload_xls_bttn()();
        let xlsFilePath = String(xls);
        if (xlsFilePath.length > 0 && type=="peptides") {
            document.getElementById("peptide-file-path").innerHTML = xlsFilePath;
        }
        if (xlsFilePath.length > 0 && type=="configs") {
            showSelectedConfigurationFilePath(xlsFilePath);
        }
        if (xlsFilePath.length > 0 && type=="counts") {
            showSelectedCountsFilePath(xlsFilePath);
        }
	}