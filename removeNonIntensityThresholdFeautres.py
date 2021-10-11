import sys
sys.path.append(".")
import Chromatogram

for fil in ["HT_SOL1_LYS_010_pos","HT_SOL1_SUP_025_pos","HT_SOL2_LYS_014_pos","HT_SOL2_SUP_029_pos","HT_SOL3_LYS_018_pos","HT_SOL3_SUP_033_pos"]:

    file = "/home/cbueschl/Documents/_backups/PeakBot/peakbot_example/MTBLS1358/Data/%s.mzXML"%fil
    chrom = Chromatogram.Chromatogram()
    chrom.parse_file(file)

    found = []
    all = 0
    with open("./%s.txt"%fil, "r") as fin:
        headers = {}
        for linei, line in enumerate(fin):
            cells = line.split("\t")
            if linei == 0:
                for celli, cell in enumerate(cells):
                    headers[cell] = celli
                found.append(line)
            else:
                mz = float(cells[headers["Precursor m/z"]].replace(",", "."))
                rt = float(cells[headers["RT (min)"]].replace(",", "."))
                rtS = float(cells[headers["RT left(min)"]].replace(",", "."))
                rtE = float(cells[headers["RT right (min)"]].replace(",", "."))
                
                scans, times, scanIDs = chrom.getSpecificArea(rtS*60, rtE*60, mz*(1-10/1E6), mz*(1+10/1E6), filterLine="", intThreshold=1E6)
                
                f = False
                for scan in scans:
                    if len(scan)>0:
                        f = True
                if f:
                    found.append(line)
                all = all + 1
                
    print("Found %d feautres, but only %d with an intensity of >= 1E6"%(all, len(found)))

    with open("./%s_proc.txt"%fil, "w") as fout:
        fout.write("\n".join(found))
                


