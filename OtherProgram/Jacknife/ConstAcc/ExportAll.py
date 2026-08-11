from ConstAcc.UsefullFunction import MeasureValueType, Export, ExportZSlice

diskname = "H:\\"

foldernames1 = ["ConstAccNew\\{}\\Polyakov\\ACCQ{}__{}_polyakov.csv",
                "ConstAccNew\\{}\\Chiral\\ACCQ{}__{}_condensateZSlicepCCLightChiralKS.csv",
                "ConstAccNew\\{}\\Chiral\\ACCQ{}__{}_condensateZSlicepCCLightCMTKSGamma3.csv",
                "ConstAccNew\\{}\\Chiral\\ACCQ{}__{}_condensateZSlicepCCLightCMTKSGamma4.csv",
                "ConstAccNewMid\\{}\\Polyakov\\ACCQM{}__{}_polyakov.csv",
                "ConstAccNewMid\\{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightChiralKS.csv",
                "ConstAccNewMid\\{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightCMTKSGamma3.csv",
                "ConstAccNewMid\\{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightCMTKSGamma4.csv"]

foldernames2 = ["ConstAccNew\\{}\\Polyakov\\ACCQ{}__{}_polyakov_ZSlice.csv",
                "ConstAccNew\\{}\\Chiral\\ACCQ{}__{}_condensateZSlicepCCLightChiralKS.csv",
                "ConstAccNew\\{}\\Chiral\\ACCQ{}__{}_condensateZSlicepCCLightCMTKSGamma3.csv",
                "ConstAccNew\\{}\\Chiral\\ACCQ{}__{}_condensateZSlicepCCLightCMTKSGamma4.csv",
                "ConstAccNewMid\\{}\\Polyakov\\ACCQM{}__{}_polyakov_ZSlice.csv",
                "ConstAccNewMid\\{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightChiralKS.csv",
                "ConstAccNewMid\\{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightCMTKSGamma3.csv",
                "ConstAccNewMid\\{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightCMTKSGamma4.csv"]

saveFold = "ExportCSV\\"

measuretypes = [MeasureValueType.POLYA,
                MeasureValueType.CHIRAL,
                MeasureValueType.G3,
                MeasureValueType.G4,
                MeasureValueType.POLYA,
                MeasureValueType.CHIRAL,
                MeasureValueType.G3,
                MeasureValueType.G4]
mids = [False, False, False, False, True, True, True, True]
startIdxs = [0, 0, 0, 0, 0, 0, 0, 0]

# i = 6
# Export(diskname + foldernames1[i], diskname + saveFold, mids[i], measuretypes[i], startIdxs[i])
# i = 4
# Export(diskname + foldernames1[i], diskname + saveFold, mids[i], measuretypes[i], startIdxs[i])

for i in range(8):
    Export(diskname + foldernames1[i], diskname + saveFold, mids[i], measuretypes[i], startIdxs[i])
#     ExportZSlice(diskname + foldernames2[i], diskname + saveFold, mids[i], measuretypes[i], startIdxs[i])