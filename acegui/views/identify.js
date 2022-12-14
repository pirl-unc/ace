function identifyPositives() {
    var readoutPath = document.getElementById("counts-file-path").innerHTML
    var configPath = document.getElementById("config-file-path").innerHTML
    let hit_df = eel.identify_positives(
        readoutPath,
        configPath
    )
}