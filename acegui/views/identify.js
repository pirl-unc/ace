function identifyPositives() {
    var readoutPath = document.getElementById("counts-file-path").innerHTML
    var configPath = document.getElementById("config-file-path").innerHTML
    let joined_output = eel.identify_positives(
        readoutPath,
        configPath
    )(renderHitPeptides)
}

function renderHitPeptides(joined_output) {

    localStorage["hit-peptides"] = JSON.stringify(joined_output[0]);
    localStorage["elispot-configuration"] = JSON.stringify(joined_output[1]);
    window.location.href="hit_peptides.html";
}