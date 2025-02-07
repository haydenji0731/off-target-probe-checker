
### convert csv to fasta format ###

import pandas as pd

# get probes
# !wget https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv

# get rid of first four rows of the csv file
# !tail -n +6 /home/caleb/Desktop/off-target-probe-checker/data/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > /home/caleb/Desktop/off-target-probe-checker/data/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A_no_header.csv

# read csv file
vis_probes = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/data/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A_no_header.csv")
vis_probes.head()


# create fasta file
with open("/home/caleb/Desktop/off-target-probe-checker/data/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.fa", "w") as f:
    for index, row in vis_probes.iterrows():
        f.write(f">{row['probe_id']}\n{row['probe_seq']}\n")
