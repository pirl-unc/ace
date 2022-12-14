var elispotConfiguration = '';

window.onload = function() {
    elispotConfiguration = JSON.parse(localStorage["elispot-configuration"]);
    console.log(elispotConfiguration);
};
