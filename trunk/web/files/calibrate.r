library("rZipHMM")

filePath = paste(path.expand("~"), "/.ziphmm.devices", sep="")

cat("Calibrating rZipHMM for multi-threaded applications.\n")
calibrate(filePath)
cat("Calibration saved in preference file: ", filePath, "\n")
