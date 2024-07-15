%pip install sanbomics

from sanbomics.plots import volcano
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------
# Fig 6A
# "volcano plot for Differentialy expressed genes in SSPE Frontal Cortex"
# importing the gene expression analysis result
#-----------------------------------------------------------------
import pandas as pd
df = pd.read_csv('DE_SSPEvsControl_bwaFC_coding_symbol.csv')

#list of interferone genes from HALLMARK_INTERFERON_ALPHA_RESPONSE from MSigDB database
# (https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.html)

interferon = """ADAR
B2M
BATF2
BST2
C1S
CASP1
CASP8
CCRL2
CD47
CD74
CMPK2
CNP
CSF1
CXCL10
CXCL11
DDX60
DHX58
EIF2AK2
ELF1
EPSTI1
MVB12A
TENT5A
CMTR1
GBP2
GBP4
GMPR
HERC6
HLA-C
IFI27
IFI30
IFI35
IFI44
IFI44L
IFIH1
IFIT2
IFIT3
IFITM1
IFITM2
IFITM3
IL15
IL4R
IL7
IRF1
IRF2
IRF7
IRF9
ISG15
ISG20
LAMP3
LAP3
LGALS3BP
LPAR6
LY6E
MOV10
MX1
NCOA7
NMI
NUB1
OAS1
OASL
OGFR
PARP12
PARP14
PARP9
PLSCR1
PNPT1
HELZ2
PROCR
PSMA3
PSMB8
PSMB9
PSME1
PSME2
RIPK2
RNF31
RSAD2
RTP4
SAMD9
SAMD9L
SELL
SLC25A28
SP110
STAT2
TAP1
TDRD7
TMEM140
TRAFD1
TRIM14
TRIM21
TRIM25
TRIM26
TRIM5
TXNIP
UBA7
UBE2L6
USP18
WARS1
"""
df.head()

interferon_list = interferon.split('\n')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define lists to hold symbols for upregulated and downregulated genes
upregulated_symbols = df[(df['adj.P.Val'] < 0.01) & (df['logFC'] >= 2)]['symbol'].tolist()  # Upregulated
downregulated_symbols = df[(df['adj.P.Val'] < 0.01) & (df['logFC'] <= -2)]['symbol'].tolist()  # Downregulated

# Assign colors to each point in the DataFrame based on the conditions
df['color'] = 'whitesmoke'  # Default color
df.loc[df['symbol'].isin(upregulated_symbols), 'color'] = 'lightgray'  # Upregulated genes
df.loc[df['symbol'].isin(downregulated_symbols), 'color'] = 'lightgray'  # Downregulated genes
df.loc[df['symbol'].isin(interferon_list), 'color'] = 'blue'  # Interferon genes

# Compute the negative logarithm of adjusted p-values
df['nlog10'] = -np.log10(df['adj.P.Val'])

# Generate the volcano plot with inverted y-axis
plt.figure(figsize=(12, 10))

# Plot all points except interferon genes first
plt.scatter(df.loc[~df['symbol'].isin(interferon_list), 'logFC'], 
            df.loc[~df['symbol'].isin(interferon_list), 'nlog10'], 
            c=df.loc[~df['symbol'].isin(interferon_list), 'color'])

# Plot interferon genes on top
plt.scatter(df.loc[df['symbol'].isin(interferon_list), 'logFC'], 
            df.loc[df['symbol'].isin(interferon_list), 'nlog10'], 
            c=df.loc[df['symbol'].isin(interferon_list), 'color'])

plt.xlabel('log2 Fold Change', fontsize=24)  # Increase font size of x-axis label
plt.ylabel('-log10(adj.P.Val)', fontsize=24)  # Increase font size of y-axis label
plt.grid(False)  # Remove grid lines

# Add vertical and horizontal dashed lines for thresholds
plt.axvline(x=2, color='gray', linestyle='--')
plt.axvline(x=-2, color='gray', linestyle='--')
plt.axhline(y=-np.log10(0.01), color='gray', linestyle='--')

# Customize axis spines and ticks
ax = plt.gca()
ax.spines['bottom'].set_linewidth(2)  # Increase thickness of bottom spine
ax.spines['left'].set_linewidth(2)  # Increase thickness of left spine
ax.spines['top'].set_visible(False)  # Turn off top spine
ax.spines['right'].set_visible(False)  # Turn off right spine

# Increase font size of x and y-axis tick labels
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

# Show the plot
plt.show()

#-----------------------------------------------------------
## "Heatmap for Top 20 GO Biological Processes"
#-----------------------------------------------------------

# Convert Regulation to numerical values for heatmap
regulation_map = {'Down-regulated': -1, 'Up-regulated': 1}
df['Regulation_Num'] = df['Regulation'].map(regulation_map)
df = df.sort_values(by='Regulation_Num', ascending=False).drop(columns='Regulation_Num')

# Pivot table for heatmap
heatmap_data = df.pivot_table(index='Term:Description', columns='Regulation', values='(-LogP)')

# Increase font size
sns.set_context("notebook", font_scale=2)

# Create heatmap
plt.figure(figsize=(8, 10))
ax = sns.heatmap(heatmap_data, cmap= 'mako_r', annot=False, cbar_kws={'shrink': 0.4, 'aspect': 7, 'label': '-Log10(P)'},  vmin=2, vmax=79)  # No annotations and shrink legend size
ax.set_facecolor('lightgray')  # Set background color
#plt.title('Regulation Heatmap')
#plt.xlabel('Regulation')

# Remove tick marks
ax.tick_params(length=0)

plt.ylabel('GO Biological Processes')
plt.show()

#---------------------------------------
## "Generating 10 neclotides before and after the BP and RI sites"
#---------------------------------------
import os
from Bio import SeqIO
import weblogo
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from weblogo import LogoOptions, LogoData, png_formatter

# Function to read the genome sequence from a FASTA file
def read_genome_sequence(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return str(record.seq)

# Function to extract sequences around breakpoints
def extract_sequences(genome_sequence, breakpoints, window_size=10):
    sequences = []
    for bp in breakpoints:
        start = bp - window_size
        end = bp + window_size + 1
        if start < 0 or end > len(genome_sequence):
            continue
        sequences.append(genome_sequence[start:end])
    return sequences

# List of breakpoints

Reinitiation = [15778,15797,15798,15796,15793,15791,15794,15791,15798,15785,15672,15802,15789,15793,15794,15802,15801,15797,15801,15799,15782,15785,15766,15730,15748,15797,15736,15793,15800,15794,15793,15798,15796,15789,15788,15735,15786,15799,15797,15765,15736,15593,15744,15794,15755,15742,15796,15488,15797,15735,15728,15727,15736,15751,15752,15751,15788,15753,15766,15769,15765,15738,15736,15753,15385,15764,15751,15752,15752,15722,15762,15740,15759,15802,15796,15741,15729,15742,15742,15748,15742,15741,15704,15742,15696,15747,15746,15762,15718,15734,15796,15745,15799,15705,15731,15682,15706,15705,15735,15773,15729,15156,15745,15744,15748,15765,15733,15729,15741,15734,15733,15748,15747,15722,15724,15706,15737,15706,15750,15719,15735,15735,15734,15755,15753,15795,15763,15752,15761,15738,15736,15703,15751,15741,15769,15748,15735,15735,15798,15797,15795,15743,15741,15744,15744,15747,15796,15751,15759,14805,15736,15751,15735,15714,15729,15759,15736,15735,15803,15801,15755,15793,15741,15792,15753,15751,15737,15790,15752,15749,15742,15735,15782,15683,15740,15739,15738,15735,15737,15734,14729,15729,15764,15763,15735,15763,15762,15690,15689,15529,15735,15755,15784,15734,15794,15789,15764,15763,15742,15723,15717,15753,15683,15760,15751,15798,15770,15795,15742,15724,15722,15735,15794,15721,15790,15735,15716,15714,15735,14431,15685,15763,15705,15742,15740,15751,15739,15738,15735,15718,15717,15704,15746,15751,15797,14270,15802,15799,15796,15719,15797,15803,14046,15800,15797,15746,15744,15735,15731,15734,15735,15742,15794,15729,15728,15241,15752,15788,15743,15706,15704,15758,13580,15800,13284,15765,15728,15681,15606,15604,15741,15735,15751,15751,15561,15048]

breakpoints = [15706,15681,15680,15682,15673,15675,15666,15669,15662,15675,15788,15652,15665,15655,15654,15646,15647,15645,15635,15637,15648,15645,15628,15658,15634,15573,15634,15577,15570,15570,15571,15566,15562,15569,15570,15623,15572,15547,15549,15575,15592,15735,15560,15510,15549,15562,15484,15792,15477,15503,15498,15499,15484,15463,15462,15451,15414,15449,15418,15415,15419,15428,15430,15383,15751,15372,15385,15378,15372,15402,15356,15372,15341,15298,15298,15335,15341,15310,15304,15298,15298,15299,15336,15292,15332,15275,15276,15254,15298,15276,15202,15241,15181,15275,15243,15256,15190,15191,15155,15117,15155,15728,15127,15128,15112,15095,15127,15113,15065,15066,15061,15040,15041,15060,15058,15064,14997,15022,14966,14991,14969,14957,14958,14913,14915,14867,14899,14910,14901,14918,14920,14953,14893,14897,14869,14890,14903,14897,14822,14817,14819,14859,14843,14828,14822,14819,14770,14815,14807,15749,14812,14797,14813,14828,14807,14753,14776,14777,14703,14705,14739,14701,14753,14696,14735,14737,14751,14698,14736,14739,14746,14753,14694,14787,14730,14731,14732,14729,14727,14730,15729,14729,14688,14689,14687,14647,14648,14720,14721,14869,14639,14607,14572,14622,14544,14549,14544,14545,14560,14573,14579,14531,14583,14458,14425,14372,14400,14375,14422,14440,14436,14423,14364,14437,14368,14411,14430,14432,14405,15685,14431,14341,14399,14356,14358,14341,14341,14324,14327,14332,14333,14340,14286,14269,14223,15750,14104,14101,14104,14163,14055,14043,15800,14046,14043,14094,14096,14099,14103,14088,14081,14050,13992,14003,14004,14479,13908,13752,13731,13720,13722,13578,15756,13488,15734,12701,12720,12689,12758,12760,12491,12473,12211,11083,10949,11096]

# Read the genome sequence from the FASTA file
genome_sequence = read_genome_sequence('SSPE_MeV_RNA.fa')

# Extract sequences around RI positions
sequences_RI = extract_sequences(genome_sequence, Reinitiation)

# Convert the list of sequences into FASTA format
fasta_sequences_RI = "\n".join([f">seq_{i}\n{seq}" for i, seq in enumerate(sequences_RI)])

print(fasta_sequences_RI)


# Extract sequences around breakpoints

sequences_BP = extract_sequences(genome_sequence, breakpoints)

# Convert the list of sequences into FASTA format
fasta_sequences_BP = "\n".join([f">seq_{i}\n{seq}" for i, seq in enumerate(sequences_BP)])

print(fasta_sequences_BP)



