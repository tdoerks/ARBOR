import os
base = "/fastscratch/tylerdoe/ARBOR/runs"
fq = base + "/fastq_20261023"
old = base + "/RVFV_inhouse_MP12/samplesheet.csv"
new = base + "/samplesheet_combined_20261026.csv"
rows = open(old).readlines()
for f in sorted(os.listdir(fq)):
    if not f.endswith("_R1_001.fastq.gz"):
        continue
    s = f.split("_S")[0]
    r1 = fq + "/" + f
    r2 = r1.replace("_R1_", "_R2_")
    rows.append(s + "_10-23-26," + r1 + "," + r2 + "\n")
open(new, "w").writelines(rows)
print("Wrote " + str(len(rows)) + " lines to " + new)
