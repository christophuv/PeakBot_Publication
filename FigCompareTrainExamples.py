#!/usr/bin/env python
# coding: utf-8

## run in python >=3.9
## activate conda environment on jucuda
## conda activate python3.9

location = "JuCuda"  ## JuCuda, HomePC

# Imports
import os
import pickle
import math
import shutil
import tempfile
import uuid
import random
import numpy as np
import pandas as pd
import plotnine as p9

os.environ["CUDA_VISIBLE_DEVICES"]="0"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import tensorflow as tf
from tensorflow.python.client import device_lib 
tf.config.experimental.set_memory_growth(tf.config.experimental.list_physical_devices('GPU')[0], True)
tf.get_logger().setLevel('WARNING')




######################################################################
###
### Load and import PeakBot
##
#

import sys
sys.path.append(os.path.join(".", "peakbot", "src"))
import peakbot.train.cuda
import peakbot.Chromatogram
            
from peakbot.core import tic, toc, tocP, tocAddStat, addFunctionRuntime, timeit, printRunTimesSummary, TabLog


def loadFile(path):
    tic()
    mzxml = None
    if os.path.exists(path + ".pickle"):
        with open(path + ".pickle", "rb") as handle:
            mzxml = pickle.load(handle)
        print("Imported pickle file for '%s' (faster than parsing mzXML file, but much larger file size)"%(inFile))
    else:
        mzxml = peakbot.Chromatogram.Chromatogram()
        mzxml.parse_file(path)
        with open(path + ".pickle", "wb") as handle:
            pickle.dump(mzxml, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Loaded chromatogram for '%s' and exported as pickle"%(inFile))
    return mzxml
    

###############################################
### Process files 
##
if __name__ == "__main__":
    tic(label="overall")

    ###############################################
    ### chromatograms to process
    inFiles = {
        "670_Sequence3_LVL1_1"  : {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1_1.mzXML" , "params": "WheatEar"},
        "670_Sequence3_LVL1_2"  : {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1_2.mzXML" , "params": "WheatEar"},

        "670_Sequence3_LVL1x2_1": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x2_1.mzXML" , "params": "WheatEar"},
        "670_Sequence3_LVL1x2_2": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x2_2.mzXML" , "params": "WheatEar"},

        "670_Sequence3_LVL1x4_1": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x4_1.mzXML" , "params": "WheatEar"},
        "670_Sequence3_LVL1x4_2": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x4_2.mzXML" , "params": "WheatEar"},
    }
    exFiles = {
        "670_Sequence3_LVL1_3"  : {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1_3.mzXML"   , "params": "WheatEar"},
        "670_Sequence3_LVL1x2_3": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x2_3.mzXML" , "params": "WheatEar"},
        "670_Sequence3_LVL1x4_3": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x4_3.mzXML" , "params": "WheatEar"},
    }
    
    detFiles = {
        "670_Sequence3_LVL1_1"  : {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1_1.mzXML"  , "params": "WheatEar"},
        "670_Sequence3_LVL1_2"  : {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1_2.mzXML"  , "params": "WheatEar"},
        "670_Sequence3_LVL1_3"  : {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1_3.mzXML"  , "params": "WheatEar"},

        #"670_Sequence3_LVL1x2_1": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x2_1.mzXML", "params": "WheatEar"},
        #"670_Sequence3_LVL1x2_2": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x2_2.mzXML", "params": "WheatEar"},
        #"670_Sequence3_LVL1x2_3": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x2_3.mzXML", "params": "WheatEar"},

        #"670_Sequence3_LVL1x4_1": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x4_1.mzXML", "params": "WheatEar"},
        #"670_Sequence3_LVL1x4_2": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x4_2.mzXML", "params": "WheatEar"},
        #"670_Sequence3_LVL1x4_3": {"file": "./peakbot_example/Data/WheatEar/670_Sequence3_LVL1x4_3.mzXML", "params": "WheatEar"},

        #"05_EB3388_AOH_p_0" : {"file": "./peakbot_example/Data/PHM/05_EB3388_AOH_p_0.mzXML" , "params": "PHM"},
        #"06_EB3389_AOH_p_10": {"file": "./peakbot_example/Data/PHM/06_EB3389_AOH_p_10.mzXML", "params": "PHM"},
        #"07_EB3390_AOH_p_20": {"file": "./peakbot_example/Data/PHM/07_EB3390_AOH_p_20.mzXML", "params": "PHM"},
        #"08_EB3391_AOH_p_60": {"file": "./peakbot_example/Data/PHM/08_EB3391_AOH_p_60.mzXML", "params": "PHM"},
        #"16_EB3392_AME_p_0":  {"file": "./peakbot_example/Data/PHM/16_EB3392_AME_p_0.mzXML" , "params": "PHM"},
        #"17_EB3393_AME_p_10": {"file": "./peakbot_example/Data/PHM/17_EB3393_AME_p_10.mzXML", "params": "PHM"},
        #"18_EB3394_AME_p_20": {"file": "./peakbot_example/Data/PHM/18_EB3394_AME_p_20.mzXML", "params": "PHM"},
        #"19_EB3395_AME_p_60": {"file": "./peakbot_example/Data/PHM/19_EB3395_AME_p_60.mzXML", "params": "PHM"}
    }
        
    
    ###############################################
    ### data parameters    
    expParams = {"WheatEar" : {"polarities": {"positive": "Q Exactive (MS lvl: 1, pol: +)"},
                               "noiseLevel":1E3, "minRT":150, "maxRT":2250, "RTpeakWidth":[8,120],
                               "intraScanMaxAdjacentSignalDifferencePPM":15, "interScanMaxSimilarSignalDifferencePPM":3,
                               "minApexBorderRatio":4, "minIntensity":1E5},
                }
    
    ## GPU information
    blockdim = 256
    griddim  = 512
    exportBatchSize= 2048
    examplesDir = "/home/users/cbueschl/_ScratchFromJuCUDA/burning_scratch/cbueschl/examples/CUDA"
    peakBotModelFile = "/home/users/cbueschl/_ScratchFromJuCUDA/burning_scratch/cbueschl/examples/CUDA/PBmodel.model.h5"
    logDir = os.path.join(".", "logs")
        
    strategy = tf.distribute.OneDeviceStrategy(device="/gpu:0")
    
    maxPopulation = 4
    intensityScales = 10
    randomnessFactor = 0.1
    
    ###############################################
    ### Generate train instances
    if True:
        headers, wepeaks       = peakbot.readTSVFile(os.path.join(".", "peakbot_example", "Reference", "WheatEar_Peaks.tsv") , convertToMinIfPossible = True)
        headers, wewalls       = peakbot.readTSVFile(os.path.join(".", "peakbot_example", "Reference", "WheatEar_Walls.tsv")      , convertToMinIfPossible = True)
        headers, webackgrounds = peakbot.readTSVFile(os.path.join(".", "peakbot_example", "Reference", "WheatEar_Backgrounds.tsv"), convertToMinIfPossible = True)
        random.shuffle(wepeaks)
        a = 2400
        wepeaksTrain = wepeaks[:a]
        wepeaksVal = wepeaks[a:]

        dsProps = {
            "T"  : {"files": inFiles , "peaks": wepeaksTrain, "walls": wewalls, "backgrounds": webackgrounds, "n": math.ceil(peakbot.Config.BATCHSIZE*peakbot.Config.STEPSPEREPOCH*peakbot.Config.EPOCHS/len(inFiles)), "shuffleSteps": 1E5},
            "V"  : {"files": inFiles , "peaks": wepeaksTrain, "walls": wewalls, "backgrounds": webackgrounds, "n": math.ceil(peakbot.Config.BATCHSIZE*peakbot.Config.STEPSPEREPOCH*128/len(inFiles))                  , "shuffleSteps": 1E4},
            "iT" : {"files": exFiles , "peaks": wepeaksVal  , "walls": wewalls, "backgrounds": webackgrounds, "n": math.ceil(peakbot.Config.BATCHSIZE*peakbot.Config.STEPSPEREPOCH*128/len(exFiles))                  , "shuffleSteps": 1E4},
            "iV" : {"files": exFiles , "peaks": wepeaksVal  , "walls": wewalls, "backgrounds": webackgrounds, "n": math.ceil(peakbot.Config.BATCHSIZE*peakbot.Config.STEPSPEREPOCH*128/len(exFiles))                  , "shuffleSteps": 1E4},
        }

        print("There are %d training examples in the T dataset"%len(wepeaksTrain))

        try:
            os.remove(os.path.join(".", "FigCompareTrainExamples.pandas.pickle"))
        except FileNotFoundError:
            pass

        random.seed(0)
        histAll = None
        for trainExamples in [len(wepeaksTrain), 2000, 1500, 1000, 500, 250, 150, 100, 50, 25, 10]:
            for rep in range(3):

                try:
                    os.remove(peakBotModelFile)
                except FileNotFoundError:
                    pass

                tic("Generated training and validation instances")

                for ds in dsProps.keys():
                    try:
                        shutil.rmtree(os.path.join(examplesDir, ds))
                    except:
                        pass
                    os.mkdir(os.path.join(examplesDir, ds))
                print("removed old training instances in '%s'"%(examplesDir))

                ###############################################
                ### Iterate files and polarities (if FPS is used)

                for ds in dsProps.keys():
                    print("Processing dataset '%s'"%ds)
                    print("")

                    for inFile, fileProps in dsProps[ds]["files"].items():
                        tic(label="sample")

                        params = expParams[fileProps["params"]]

                        ###############################################
                        ### data parameters for chromatograms
                        polarities = params["polarities"]
                        intraScanMaxAdjacentSignalDifferencePPM = params["intraScanMaxAdjacentSignalDifferencePPM"]
                        interScanMaxSimilarSignalDifferencePPM = params["interScanMaxSimilarSignalDifferencePPM"]
                        RTpeakWidth = params["RTpeakWidth"]
                        minApexBorderRatio = params["minApexBorderRatio"]
                        minIntensity = params["minIntensity"]

                        for polarity, filterLine in polarities.items():
                            print("Processing sample '%s', polarity '%s'"%(inFile, polarity))
                            print("")

                            ###############################################
                            ### Load chromatogram
                            tic()
                            mzxml = loadFile(fileProps["file"])
                            mzxml.keepOnlyFilterLine(filterLine)
                            print("Filtered mzXML file for %s scan events only"%(polarity))
                            print("  | .. took %.1f seconds"%(toc()))
                            print("")

                            ###############################################
                            ### Generate train data
                            peaks = dsProps[ds]["peaks"]
                            if ds == "T":
                                print("  | .. Restricting peaks to %d examples"%trainExamples)
                                peaks = random.sample(dsProps[ds]["peaks"], trainExamples)
                            peakbot.train.cuda.generateTestInstances(mzxml, "'%s':'%s'"%(inFile, filterLine), 
                                                                     peaks, dsProps[ds]["walls"], dsProps[ds]["backgrounds"],
                                                                     nTestExamples = dsProps[ds]["n"], exportPath = os.path.join(examplesDir, ds), 
                                                                     intraScanMaxAdjacentSignalDifferencePPM=intraScanMaxAdjacentSignalDifferencePPM, 
                                                                     interScanMaxSimilarSignalDifferencePPM=interScanMaxSimilarSignalDifferencePPM, 

                                                                     updateToLocalPeakProperties = True,
                                                                     RTpeakWidth = RTpeakWidth, minApexBorderRatio = minApexBorderRatio, minIntensity = minIntensity, 

                                                                     maxPopulation = maxPopulation, intensityScales = intensityScales, randomnessFactor = randomnessFactor, 
                                                                     blockdim = blockdim, griddim = griddim,
                                                                     verbose = True)

                    print("\n\n\n\n\n")        

                peakbot.train.shuffleResultsSampleNames(os.path.join(examplesDir, ds), verbose = True)
                peakbot.train.shuffleResults(os.path.join(examplesDir, ds), steps = dsProps[ds]["shuffleSteps"], samplesToExchange = 50, verbose = True)

                tocP("Generated training and validation instances", label="Generated training and validation instances")
                print("\n\n\n\n\n")





                ###############################################
                ### Train new PeakBot Model
                tic("train new PeakBot model")
                pb = None
                with strategy.scope():

                    addValDS = []
                    for ds in dsProps.keys():
                        addValDS.append({"folder": os.path.join(examplesDir, ds), "name": ds, "numBatches": 8})

                    pb, hist = peakbot.trainPeakBotModel(trainInstancesPath = os.path.join(examplesDir, "T"), 
                                                         addValidationInstances = addValDS,
                                                         logBaseDir = logDir, 
                                                         verbose = True)

                    pb.saveModelToFile(peakBotModelFile)
                    print("Newly trained peakbot saved to file '%s'"%(peakBotModelFile))

                    hist.drop("model", axis=1, inplace=True)
                    hist.insert(0, "TrainExamples", [trainExamples for i in range(hist.shape[0])], True)
                    hist.insert(0, "Rep", [rep for i in range(hist.shape[0])], True)

                    if histAll is None:
                        histAll = hist
                    else:
                        histAll = histAll.append(hist, ignore_index=True)

                    histAll.to_pickle(os.path.join(".", "FigCompareTrainExamples.pandas.pickle"))
                    print("")
                tocP("train new PeakBot model","train new PeakBot model")
                print("\n\n\n\n\n")




                ###############################################
                ### Predict in new data
                for inFile, fileProps in detFiles.items():
                    tic(label="sample")

                    params = expParams[fileProps["params"]]

                    ###############################################
                    ### Load chromatogram
                    tic("sample")
                    tic()
                    mzxml = None
                    if os.path.exists(fileProps["file"] + ".pickle"):
                        with open(fileProps["file"] + ".pickle", "rb") as handle:
                            mzxml = pickle.load(handle)
                        print("Imported pickle file for '%s' (faster than parsing mzXML file, but much larger file size)"%(inFile))
                    else:
                        mzxml = peakbot.Chromatogram.Chromatogram()
                        mzxml.parse_file(fileProps["file"])
                        with open(fileProps["file"] + ".pickle", "wb") as handle:
                            pickle.dump(mzxml, handle, protocol=pickle.HIGHEST_PROTOCOL)
                        print("Loaded chromatogram for '%s' and exported as pickle"%(inFile))
                    print("  | .. took %.1f seconds"%(toc()))
                    print("")



                    ###############################################
                    ### data parameters for chromatograms
                    polarities = params["polarities"]
                    noiseLevel = params["noiseLevel"]
                    minRT = params["minRT"]
                    maxRT = params["maxRT"]
                    intraScanMaxAdjacentSignalDifferencePPM = params["intraScanMaxAdjacentSignalDifferencePPM"]
                    interScanMaxSimilarSignalDifferencePPM = params["interScanMaxSimilarSignalDifferencePPM"]
                    RTpeakWidth = params["RTpeakWidth"]
                    minApexBorderRatio = params["minApexBorderRatio"]
                    minIntensity = params["minIntensity"]

                    for polarity, filterLine in polarities.items():
                        print("Processing sample '%s', polarity '%s'"%(inFile, polarity))
                        with tempfile.TemporaryDirectory() as tmpdirname:
                            tic("instance")

                            ###############################################
                            ### Preprocess chromatogram
                            tic()
                            mzxml.keepOnlyFilterLine(filterLine)
                            print("Filtered mzXML file for %s scan events only\n  | .. took %.1f seconds"%(polarity, toc()))
                            print("")

                            tic()
                            mzxml.removeBounds(minRT=minRT, maxRT=maxRT, minMZ=100)
                            mzxml.removeNoise(noiseLevel)
                            print("Removed noise (%g) and bounds\n  | .. took %.1f seconds"%(noiseLevel, toc()))
                            print("")



                            ###############################################
                            ### Detect local maxima with peak-like shapes## CUDA-GPU
                            tic(label="preProcessing")
                            peaks, maximaProps, maximaPropsAll = peakbot.cuda.preProcessChromatogram(
                                    mzxml, "'%s':'%s'"%(inFile, filterLine), 
                                    intraScanMaxAdjacentSignalDifferencePPM = intraScanMaxAdjacentSignalDifferencePPM,
                                    interScanMaxSimilarSignalDifferencePPM = interScanMaxSimilarSignalDifferencePPM,
                                    RTpeakWidth = RTpeakWidth,
                                    minApexBorderRatio = minApexBorderRatio,
                                    minIntensity = minIntensity,
                                    exportPath = tmpdirname, 
                                    exportLocalMaxima = "peak-like-shape", # "all", "localMaxima-with-mzProfile", "peak-like-shape"
                                    exportBatchSize = exportBatchSize,
                                    blockdim = blockdim,
                                    griddim  = griddim, 
                                    verbose = True)
                            expPeaks = peaks # peaks # select which local maxima to use: peaks, maximaProps, maximaPropsAll


                            ###############################################
                            ### Detect peaks with PeakBot
                            tic("PeakBotDetection")
                            with strategy.scope():
                                peaks, walls, backgrounds, errors = peakbot.runPeakBot(tmpdirname, peakBotModelFile)
                                df = pd.DataFrame([[rep, trainExamples, "Predicted", "nPeaks %s"%inFile, len(peaks)]], columns=["Rep", "TrainExamples", "set", "metric", "value"])
                                histAll = histAll.append(df)
                                histAll.to_pickle(os.path.join(".", "FigCompareTrainExamples.pandas.pickle"))
                            tocP("predicted with PeakBot", label="PeakBotDetection")
                            print("")

        
        
        
    ###############################################
    ### Show training metrices
    df = pd.read_pickle(os.path.join(".", "FigCompareTrainExamples.pandas.pickle"))
    df = df[df["metric"]!="loss"]
    df["metric"] = df["metric"].str.replace("07_EB3390_AOH_p_20", "07.EB3390.AOH.p.20")
    df["metric"] = df["metric"].str.replace("08_EB3391_AOH_p_60", "08.EB3391.AOH.p.60")
    df["metric"] = df["metric"].str.replace("670_Sequence3_LVL1_1", "670.Sequence3.LVL1.1")
    df["metric"] = df["metric"].str.replace("670_Sequence3_LVL1_2", "670.Sequence3.LVL1.2")
    df["metric"] = df["metric"].str.replace("670_Sequence3_LVL1_3", "670.Sequence3.LVL1.3")
    df["metric"] = df["metric"].str.replace("_", "\n")
    df["metric"] = df["metric"].str.replace("07.EB3390.AOH.p.20", "07_EB3390_AOH_p_20")
    df["metric"] = df["metric"].str.replace("08.EB3391.AOH.p.60", "08_EB3391_AOH_p_60")
    df["metric"] = df["metric"].str.replace("670.Sequence3.LVL1.1", "670_Sequence3_LVL1_1")
    df["metric"] = df["metric"].str.replace("670.Sequence3.LVL1.2", "670_Sequence3_LVL1_2")
    df["metric"] = df["metric"].str.replace("670.Sequence3.LVL1.3", "670_Sequence3_LVL1_3")
    df["metric"] = df["metric"].str.replace("nPeaks ", "Number of features\n")

    plot = (p9.ggplot(df[-df["metric"].str.startswith("Number of features")], 
                      p9.aes("factor(TrainExamples)", "value", color="metric", group="factor(TrainExamples)"))
            + p9.facet_grid("metric~set", scales="free_y")
            + p9.geom_boxplot(p9.aes("factor(TrainExamples)", "value", group="factor(TrainExamples)"))
            + p9.geom_jitter()
            + p9.stat_summary(geom="line")
            + p9.ggtitle("Comparison of different sets for training")
            + p9.xlab("Number of training examples used for training a new PeakBot model")
            + p9.ylab("Loss/Metric value")
            + p9.theme(legend_position="none", axis_text_x=p9.element_text(angle=45)))
    p9.options.figure_size = (12, 12)
    p9.ggsave(plot=plot, filename="./FigCompareTrainExamples_Metrics.png", height=12, width=12, dpi=300)

    plot = (p9.ggplot(df[df["metric"].str.startswith("Number of features")], 
                      p9.aes("factor(TrainExamples)", "value", color="metric", group="factor(TrainExamples)"))
            + p9.facet_grid("metric~.", scales="free_y")
            + p9.geom_boxplot(p9.aes("factor(TrainExamples)", "value", group="factor(TrainExamples)"))
            + p9.geom_jitter()
            + p9.stat_summary(geom="line")
            + p9.ggtitle("Comparison of different sets for training")
            + p9.xlab("Number of training examples used for training a new PeakBot model")
            + p9.ylab("Number of detected features")
            + p9.theme(legend_position="none"))
    p9.options.figure_size = (19, 19)
    p9.ggsave(plot=plot, filename="./FigCompareTrainExamples_Detected.png", height=12, width=12, dpi=300)


    print("\n\n\n\n\n\n\n\n\n\n")

 