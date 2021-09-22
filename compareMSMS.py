#!/usr/bin/env python
# coding: utf-8

## run in python >= 3.8
## activate conda environment on jucuda
## conda activate python3.8

# Imports
import os
import pickle
import tempfile
import tqdm
import numpy as np
from numba import jit
import numba
from collections import OrderedDict

os.environ["CUDA_VISIBLE_DEVICES"]="0"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
tf.config.experimental.set_memory_growth(tf.config.experimental.list_physical_devices('GPU')[0], True)

import sys
sys.path.append(os.path.join(".", "peakbot", "src"))
import peakbot
import peakbot.Chromatogram
import peakbot.cuda            
from peakbot.core import tic, toc, TabLog


def loadFile(path):
    tic()
    mzxml = None
    if os.path.exists(path+".pickle"):
        with open(path+".pickle", "rb") as inF:
            mzxml = pickle.load(inF)
        print("Imported chromatogram.pickle for '%s'"%(path))
    else:
        mzxml = peakbot.Chromatogram.Chromatogram()
        mzxml.parse_file(path)
        with open(path+".pickle", "wb") as outF:
            pickle.dump(mzxml, outF)
        print("Imported chromatogram for '%s'"%(path))
    return mzxml


@jit(nopython = True)
def match(msmss, peaks, msmsforPeak, walls):
    
    _msmsOfPeaks = 0
    _msmsOfWalls = 0
    _msmsOfUnspe = 0

    for i in range(msmss.shape[0]):
        usedForPeak = False
        for j in range(peaks.shape[0]):
            if peaks[j,2] <= msmss[i,0] <= peaks[j,3] and peaks[j,4] <= msmss[i,1] <= peaks[j,5]:
                msmsforPeak[j] = 1
                _msmsOfPeaks += 1
                usedForPeak = True
                break
        
        usedForWall = False
        if not usedForPeak:
            for j in range(walls.shape[0]):
                if abs(msmss[i,0] - walls[j,0]) <= 5 and abs(msmss[i,1] - walls[j,1]) * 1E6/walls[j,1] <= 10:
                    _msmsOfWalls += 1
                    usedForWall = True
                    break

        if not usedForPeak and not usedForWall:
            _msmsOfUnspe += 1

    return _msmsOfPeaks, _msmsOfWalls, _msmsOfUnspe



###############################################
### Process files 
##
if __name__ == "__main__":
    tic(label="overall")

    ###############################################
    ### data parameters
    ##
    ## Different LC-HRMS settings can be used per chromatogram. To be able to reuse them, they are summarized
    ##    in dictionaries. The keys are then used as the setting values
    ##
    ## polarities: specifies which filter lines are to be used for detecting the chromatographic peaks
    ## noiseLevel: Everything below this threshold is considered noise and removed directly after the import
    ## minRT / maxRT: Area of the chromatogram in which chromatographic peaks are expected
    ## RTpeakWidth: array of [minimum, maximum] peak-width in scans
    ## SavitzkyGolayWindowPlusMinus: specifies the degree of smoothing. A value of x results in a smoothing window of 2*x + 1
    ## intraScanMaxAdjacentSignalDifferencePPM: Maximum difference of signals belonging to the same profile mode peak
    ## interScanMaxSimilarSignalDifferencePPM: Maximum difference of signals representing the same profile mode signal
    ## minIntensity: All signals below this threshold are not considered for the local maximum detection
    expParams = {"WheatEar": {"polarities": {"negative": "Q Exactive (MS lvl: 1, pol: -)", "positive": "Q Exactive (MS lvl: 1, pol: +)"},
                              "PeakBotModel": "./temp/PBmodel_WheatEar.model.h5",
                              "minRT":150, "maxRT":2250, "RTpeakWidth":[8,120], "SavitzkyGolayWindowPlusMinus": 3,
                              "intraScanMaxAdjacentSignalDifferencePPM":15, "interScanMaxSimilarSignalDifferencePPM":3,
                              "noiseLevel":1E3, "minIntensity":1E5},
                }

    ###############################################
    ### chromatograms to process
    ##
    ## Different LC-HRMS chromatograms can be used for generating a training or validation dataset
    ##
    ## file: Path of the mzXML file
    ## params: parameter collection for the particular sample (see variable expParams)
    inFiles = OrderedDict()

    # Unreleased data
    # Wheat ear (similar to untreated samples of: Stable Isotope–Assisted Plant Metabolomics: Combination of Global and Tracer-Based Labeling for Enhanced Untargeted Profiling and Compound Annotation)
    # https://doi.org/10.3389/fpls.2019.01366
    inFiles["823_SampleLvl1_Exp670_ddMSMS"                     ] = {"file": "./peakbot_example/Data/WheatEar/823_SampleLvl1_Exp670_ddMSMS.mzXML"                       , "params": "WheatEar"}
    #inFiles["823_sampleLvl1_Exp670_ddMSMS_withExclusionList"   ] = {"file": "./peakbot_example/Data/WheatEar/823_sampleLvl1_Exp670_ddMSMS_withExclusionList.mzXML"     , "params": "WheatEar"}
    inFiles["823_sampleLvl1_Exp670_ddMSMS_withExclusionList_V2"] = {"file": "./peakbot_example/Data/WheatEar/823_sampleLvl1_Exp670_ddMSMS_withExclusionList_V2.mzXML"  , "params": "WheatEar"}

    ###############################################
    ## GPU information
    ##
    ## These values specify how the GPU is used for generating the training examples
    ## Please consult the documentation of your GPU.
    ## Values for an old Nvidia GTX 970 graphics card with 4GB GPU-memory are blockdim = 256, griddim = 64
    ## These should thus work for most newer card, however, for maximum performance these should be optimized to the GPU used
    ## The strategy specifies on which device tensorflow shall be executed.
    ## exportBatchSize: specifies how many putative areas shall be exported in one batch
    ## peakBotModelFile: specifies which model to load from the file system
    blockdim = 512
    griddim  = 256
    strategy = tf.distribute.OneDeviceStrategy(device="/gpu:0")
    exportBatchSize = 2048
    peakBotModelFile = "./peakbot_example/temp/PBmodel_WheatEar.model.h5"        
    
    
    ###############################################
    ## Finisehd with specifying LC-HRMS chromatogram files and LC-HRMS settings
    ## Nothing to change from here on
    

    ###############################################
    ### Iterate files and polarities (if FPS is used)
    for inFile, fileProps in inFiles.items():
        tic(label="sample")
    
        ###############################################
        ### data parameters for chromatograms        
        params = expParams[fileProps["params"]]
        
        for polarity, filterLine in params["polarities"].items():
            print("Processing sample '%s', polarity '%s'"%(inFile, polarity))
            walls = []
            backgrounds = []
            errors = []
            tic("instance")
            
            ###############################################
            ### Preprocess chromatogram
            tic("sample")
            mzxml = loadFile(fileProps["file"])
            print(mzxml.getFilterLines())
            mzxml.keepOnlyFilterLine(filterLine)
            print("  | .. filtered mzXML file for %s scan events only"%(polarity))
            mzxml.removeBounds(minRT = params["minRT"], maxRT = params["maxRT"])
            mzxml.removeNoise(params["noiseLevel"])
            print("  | .. removed noise and bounds")
            print("")

            with tempfile.TemporaryDirectory() as tmpdirname:
                ###############################################
                ### Detect high-quality local maxima with peak-like shapes with GD and PeakBot
                peaks, maximaProps, maximaPropsAll = peakbot.cuda.preProcessChromatogram(
                        mzxml, "'%s':'%s'"%(inFile, filterLine), 
                        intraScanMaxAdjacentSignalDifferencePPM = params["intraScanMaxAdjacentSignalDifferencePPM"],
                        interScanMaxSimilarSignalDifferencePPM = params["interScanMaxSimilarSignalDifferencePPM"],
                        RTpeakWidth = params["RTpeakWidth"],
                        SavitzkyGolayWindowPlusMinus = params["SavitzkyGolayWindowPlusMinus"], 
                        minIntensity = params["minIntensity"],
                        exportPath = tmpdirname, 
                        exportLocalMaxima = "all",
                        exportBatchSize = exportBatchSize, 
                        blockdim = blockdim,
                        griddim  = griddim, 
                        verbose = True)
                tic("PeakBotDetection")
                peaks = []
                walls = []
                backgrounds = []
                with strategy.scope():
                    peaks, walls, backgrounds, errors_ = peakbot.runPeakBot(tmpdirname, peakBotModelFile)
                print("")
                
                ###############################################
                ### Postprocessing
                tic("postProcessing")
                peaks = peakbot.cuda.postProcess(mzxml, "'%s':'%s'"%(inFile, filterLine), peaks, 
                                                 blockdim = blockdim,
                                                 griddim  = griddim, 
                                                 verbose = True)
                print("")
                
                ## Log features
                TabLog().addData("%s - %s"%(inFile, filterLine), "Features", len(peaks))
                TabLog().addData("%s - %s"%(inFile, filterLine), "Walls", len(walls))
                TabLog().addData("%s - %s"%(inFile, filterLine), "Backgrounds", len(backgrounds))
                TabLog().addData("%s - %s"%(inFile, filterLine), "Errors", len(errors))
                print("")

                mzxml = loadFile(fileProps["file"])
                msmsScans = []
                for msms in mzxml.MS2_list:
                    if (msms.polarity == "+" and polarity == "positive") or (msms.polarity == "-" and polarity == "negative"):
                        if params["minRT"] <= msms.retention_time <= params["maxRT"]:
                            msmsScans.append([msms.retention_time, msms.precursor_mz])
                TabLog().addData("%s - %s"%(inFile, filterLine), "MSMS scans", len(msmsScans))

                msmsforPeak = np.zeros((len(peaks)))
                msmsOfPeaks, msmsOfWalls, msmsOfUnspecific = match(np.array(msmsScans, np.float32), np.array(peaks, np.float32), msmsforPeak, np.array(walls, np.float32))
                
                TabLog().addData("%s - %s"%(inFile, filterLine), "MSMS 4 peaks", msmsOfPeaks)
                TabLog().addData("%s - %s"%(inFile, filterLine), "Peaks with MSMS", np.sum(msmsforPeak))
                TabLog().addData("%s - %s"%(inFile, filterLine), "MSMS 4 walls", msmsOfWalls)
                TabLog().addData("%s - %s"%(inFile, filterLine), "MSMS unspecific", msmsOfUnspecific)

                print("\n\n\n\n\n")
        
    TabLog().addData("Total time all files", "time (sec)", "%.1f"%toc("overall"))
    print("")
    TabLog().print()
