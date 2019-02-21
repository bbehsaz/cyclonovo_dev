import sys, string


def readStringMGF(peaksFile, scanID):
    protonMass = 1.00728
    peaks = {}
    peaksnIntensity = {}
    pepMasses = {}
    retention = -1
    charges = {}
    n = -1
    retentions = {}
    specLines = {}
    fileLines = {}
    # Reading the peaks and intensities.
    alllines = peaksFile.split("\n")
    pointer = 0
    allSpectraVector = {}
    while pointer < len(alllines):
        # line = peaksFile.readline()
        line = alllines[pointer]
        if line != "BEGIN IONS":
            pointer += 1
            continue
        originalLine = line
        line = line.strip()
        if line == "BEGIN IONS":
            # initiate an array with 2,000,000 entires of 0
            spectrumVector = [0] * 2000000
            intensityVector = [0] * 2000000
            specLines = [originalLine + "\n"]
            peptide = ""
            while (True):
                pointer += 1
                line = alllines[pointer]
                if not line:
                    continue
                # line = peaksFile.readline()
                if not line[0].isdigit():
                    specLines.append(line + "\n")
                    if line[0:6] == "CHARGE":
                        charge = int(''.join(c for c in line.split("=")[1] if c.isdigit()))
                    if line[0:6] == "PEPMAS":
                        pepMass = float(line.strip().split()[0][8:])
                    if line[0:5] == "DBID=":
                        line = line.strip()
                        if line[5:] != "":
                            antimartin = line[15:]
                    if "TITLE=" in line:
                        peptide = line.strip().split("=")[1]
                    if "RTINSECONDS" in line:
                        retention = line.split("=")[1].strip()
                else:
                    break
            n += 1
            # if peptide == "":
            # 	peptide=str(n)
            peptide = str(scanID)
            # peptide = "ion" + "_"+str(n)
            retentions[peptide] = retention
            peaksnIntensity[peptide] = {}
            charges[peptide] = charge
            pepMasses[peptide] = pepMass
            pointer2 = 0
            while line != "END IONS":

                specLines.append(line.strip() + "\n")
                peakLine = line.strip().split()
                peakMass = round(float(peakLine[0]), 3)
                intensity = float(peakLine[1])
                # intensity =1
                peaksnIntensity[peptide][peakMass] = intensity
                if charge == 1:
                    pointer2 = int(round(peakMass, 3) * 1000)
                    if pointer2 < 2000001:
                        spectrumVector[pointer2] = intensity
                # intensityVector[round(peakMass,3)] = 1
                elif charge == 2:
                    charge1PeakMass = int(round((2 * peakMass) - protonMass, 3) * 1000)
                    # peaksnIntensity[peptide][charge1PeakMass]=intensity
                    if charge1PeakMass < 2000001:
                        spectrumVector[charge1PeakMass] = intensity
                # spectrumVector[int(round(peakMass,3)*1000)] = intensity
                pointer += 1
                line = alllines[pointer].strip()

            if line == "END IONS":
                allSpectraVector[peptide] = spectrumVector
                specLines.append(line + "\n")
                fileLines[peptide] = specLines
                continue

    return peaksnIntensity, pepMasses, charges, retentions, fileLines, allSpectraVector
