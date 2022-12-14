
async function getCSV() {
	var csv = await eel.upload_csv_bttn()();
        let csvFilePath = String(csv);
        if (csvFilePath.length > 0) {
            document.getElementById("peptide-sequence-file-path").innerHTML = csvFilePath;
        }
	}

async function getXLS(type) {
	var xls = await eel.upload_xls_bttn()();
        let xlsFilePath = String(xls);
        if (xlsFilePath.length > 0 && type=="configs") {
            document.getElementById("config-file-path").innerHTML = xlsFilePath;
        }
        if (xlsFilePath.length > 0 && type=="counts") {
            document.getElementById("counts-file-path").innerHTML = xlsFilePath;
        }
	}