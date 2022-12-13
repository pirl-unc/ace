window.onload = function() {
    let isSuccessful = eel.download_probert_data()(goToLandingPage)
};

function goToLandingPage(isSuccessful) {
    if (isSuccessful) {
        window.location.href="landing.html";
    } else {
        window.alert("There was an issue with downloading the necessary data.");
    }
}
