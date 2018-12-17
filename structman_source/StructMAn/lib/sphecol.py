import os
import sys
import getopt
import traceback
import math
import numpy
import time

threeToOne = {
    "00C": "C", "01W": "X", "02K": "A", "03Y": "C", "07O": "C",
    "08P": "C", "0A0": "D", "0A1": "Y", "0A2": "K", "0A8": "C",
    "0AA": "V", "0AB": "V", "0AC": "G", "0AD": "G", "0AF": "W",
    "0AG": "L", "0AH": "S", "0AK": "D", "0AM": "A", "0AP": "C",
    "0AU": "U", "0AV": "A", "0AZ": "P", "0BN": "F", "0C ": "C",
    "0CS": "A", "0DC": "C", "0DG": "G", "0DT": "T", "0FL": "A",
    "0G ": "G", "0NC": "A", "0SP": "A", "0U ": "U", "0YG": "YG",
    "10C": "C", "125": "U", "126": "U", "127": "U", "128": "N",
    "12A": "A", "143": "C", "175": "ASG", "193": "X", "1AP": "A",
    "1MA": "A", "1MG": "G", "1PA": "F", "1PI": "A", "1PR": "N",
    "1SC": "C", "1TQ": "W", "1TY": "Y", "1X6": "S", "200": "F",
    "23F": "F", "23S": "X", "26B": "T", "2AD": "X", "2AG": "A",
    "2AO": "X", "2AR": "A", "2AS": "X", "2AT": "T", "2AU": "U",
    "2BD": "I", "2BT": "T", "2BU": "A", "2CO": "C", "2DA": "A",
    "2DF": "N", "2DM": "N", "2DO": "X", "2DT": "T", "2EG": "G",
    "2FE": "N", "2FI": "N", "2FM": "M", "2GT": "T", "2HF": "H",
    "2LU": "L", "2MA": "A", "2MG": "G", "2ML": "L", "2MR": "R",
    "2MT": "P", "2MU": "U", "2NT": "T", "2OM": "U", "2OT": "T",
    "2PI": "X", "2PR": "G", "2SA": "N", "2SI": "X", "2ST": "T",
    "2TL": "T", "2TY": "Y", "2VA": "V", "2XA": "C", "32S": "X",
    "32T": "X", "3AH": "H", "3AR": "X", "3CF": "F", "3DA": "A",
    "3DR": "N", "3GA": "A", "3MD": "D", "3ME": "U", "3NF": "Y",
    "3QN": "K", "3TY": "X", "3XH": "G", "4AC": "N", "4BF": "Y",
    "4CF": "F", "4CY": "M", "4DP": "W", "4F3": "GYG", "4FB": "P",
    "4FW": "W", "4HT": "W", "4IN": "W", "4MF": "N", "4MM": "X",
    "4OC": "C", "4PC": "C", "4PD": "C", "4PE": "C", "4PH": "F",
    "4SC": "C", "4SU": "U", "4TA": "N", "4U7": "A", "56A": "H",
    "5AA": "A", "5AB": "A", "5AT": "T", "5BU": "U", "5CG": "G",
    "5CM": "C", "5CS": "C", "5FA": "A", "5FC": "C", "5FU": "U",
    "5HP": "E", "5HT": "T", "5HU": "U", "5IC": "C", "5IT": "T",
    "5IU": "U", "5MC": "C", "5MD": "N", "5MU": "U", "5NC": "C",
    "5PC": "C", "5PY": "T", "5SE": "U", "5ZA": "TWG", "64T": "T",
    "6CL": "K", "6CT": "T", "6CW": "W", "6HA": "A", "6HC": "C",
    "6HG": "G", "6HN": "K", "6HT": "T", "6IA": "A", "6MA": "A",
    "6MC": "A", "6MI": "N", "6MT": "A", "6MZ": "N", "6OG": "G",
    "70U": "U", "7DA": "A", "7GU": "G", "7JA": "I", "7MG": "G",
    "8AN": "A", "8FG": "G", "8MG": "G", "8OG": "G", "9NE": "E",
    "9NF": "F", "9NR": "R", "9NV": "V", "A  ": "A", "A1P": "N",
    "A23": "A", "A2L": "A", "A2M": "A", "A34": "A", "A35": "A",
    "A38": "A", "A39": "A", "A3A": "A", "A3P": "A", "A40": "A",
    "A43": "A", "A44": "A", "A47": "A", "A5L": "A", "A5M": "C",
    "A5N": "N", "A5O": "A", "A66": "X", "AA3": "A", "AA4": "A",
    "AAR": "R", "AB7": "X", "ABA": "A", "ABR": "A", "ABS": "A",
    "ABT": "N", "ACB": "D", "ACL": "R", "AD2": "A", "ADD": "X",
    "ADX": "N", "AEA": "X", "AEI": "D", "AET": "A", "AFA": "N",
    "AFF": "N", "AFG": "G", "AGM": "R", "AGT": "C", "AHB": "N",
    "AHH": "X", "AHO": "A", "AHP": "A", "AHS": "X", "AHT": "X",
    "AIB": "A", "AKL": "D", "AKZ": "D", "ALA": "A", "ALC": "A",
    "ALM": "A", "ALN": "A", "ALO": "T", "ALQ": "X", "ALS": "A",
    "ALT": "A", "ALV": "A", "ALY": "K", "AN8": "A", "AP7": "A",
    "APE": "X", "APH": "A", "API": "K", "APK": "K", "APM": "X",
    "APP": "X", "AR2": "R", "AR4": "E", "AR7": "R", "ARG": "R",
    "ARM": "R", "ARO": "R", "ARV": "X", "AS ": "A", "AS2": "D",
    "AS9": "X", "ASA": "D", "ASB": "D", "ASI": "D", "ASK": "D",
    "ASL": "D", "ASM": "X", "ASN": "N", "ASP": "D", "ASQ": "D",
    "ASU": "N", "ASX": "B", "ATD": "T", "ATL": "T", "ATM": "T",
    "AVC": "A", "AVN": "X", "AYA": "A", "AYG": "AYG", "AZK": "K",
    "AZS": "S", "AZY": "Y", "B1F": "F", "B1P": "N", "B2A": "A",
    "B2F": "F", "B2I": "I", "B2V": "V", "B3A": "A", "B3D": "D",
    "B3E": "E", "B3K": "K", "B3L": "X", "B3M": "X", "B3Q": "X",
    "B3S": "S", "B3T": "X", "B3U": "H", "B3X": "N", "B3Y": "Y",
    "BB6": "C", "BB7": "C", "BB8": "F", "BB9": "C", "BBC": "C",
    "BCS": "C", "BE2": "X", "BFD": "D", "BG1": "S", "BGM": "G",
    "BH2": "D", "BHD": "D", "BIF": "F", "BIL": "X", "BIU": "I",
    "BJH": "X", "BLE": "L", "BLY": "K", "BMP": "N", "BMT": "T",
    "BNN": "F", "BNO": "X", "BOE": "T", "BOR": "R", "BPE": "C",
    "BRU": "U", "BSE": "S", "BT5": "N", "BTA": "L", "BTC": "C",
    "BTR": "W", "BUC": "C", "BUG": "V", "BVP": "U", "BZG": "N",
    "C  ": "C", "C12": "TYG", "C1X": "K", "C25": "C", "C2L": "C",
    "C2S": "C", "C31": "C", "C32": "C", "C34": "C", "C36": "C",
    "C37": "C", "C38": "C", "C3Y": "C", "C42": "C", "C43": "C",
    "C45": "C", "C46": "C", "C49": "C", "C4R": "C", "C4S": "C",
    "C5C": "C", "C66": "X", "C6C": "C", "C99": "TFG", "CAF": "C",
    "CAL": "X", "CAR": "C", "CAS": "C", "CAV": "X", "CAY": "C",
    "CB2": "C", "CBR": "C", "CBV": "C", "CCC": "C", "CCL": "K",
    "CCS": "C", "CCY": "CYG", "CDE": "X", "CDV": "X", "CDW": "C",
    "CEA": "C", "CFL": "C", "CFY": "FCYG", "CG1": "G", "CGA": "E",
    "CGU": "E", "CH ": "C", "CH6": "MYG", "CH7": "KYG", "CHF": "X",
    "CHG": "X", "CHP": "G", "CHS": "X", "CIR": "R", "CJO": "GYG",
    "CLE": "L", "CLG": "K", "CLH": "K", "CLV": "AFG", "CM0": "N",
    "CME": "C", "CMH": "C", "CML": "C", "CMR": "C", "CMT": "C",
    "CNU": "U", "CP1": "C", "CPC": "X", "CPI": "X", "CQR": "GYG",
    "CR0": "TLG", "CR2": "GYG", "CR5": "G", "CR7": "KYG", "CR8": "HYG",
    "CRF": "TWG", "CRG": "THG", "CRK": "MYG", "CRO": "GYG", "CRQ": "QYG",
    "CRU": "EYG", "CRW": "ASG", "CRX": "ASG", "CS0": "C", "CS1": "C",
    "CS3": "C", "CS4": "C", "CS8": "N", "CSA": "C", "CSB": "C",
    "CSD": "C", "CSE": "C", "CSF": "C", "CSH": "SHG", "CSI": "G",
    "CSJ": "C", "CSL": "C", "CSO": "C", "CSP": "C", "CSR": "C",
    "CSS": "C", "CSU": "C", "CSW": "C", "CSX": "C", "CSY": "SYG",
    "CSZ": "C", "CTE": "W", "CTG": "T", "CTH": "T", "CUC": "X",
    "CWR": "S", "CXM": "M", "CY0": "C", "CY1": "C", "CY3": "C",
    "CY4": "C", "CYA": "C", "CYD": "C", "CYF": "C", "CYG": "C",
    "CYJ": "X", "CYM": "C", "CYQ": "C", "CYR": "C", "CYS": "C",
    "CZ2": "C", "CZO": "GYG", "CZZ": "C", "D11": "T", "D1P": "N",
    "D3 ": "N", "D33": "N", "D3P": "G", "D3T": "T", "D4M": "T",
    "D4P": "X", "DA ": "A", "DA2": "X", "DAB": "A", "DAH": "F",
    "DAL": "A", "DAR": "R", "DAS": "D", "DBB": "T", "DBM": "N",
    "DBS": "S", "DBU": "T", "DBY": "Y", "DBZ": "A", "DC ": "C",
    "DC2": "C", "DCG": "G", "DCI": "X", "DCL": "X", "DCT": "C",
    "DCY": "C", "DDE": "H", "DDG": "G", "DDN": "U", "DDX": "N",
    "DFC": "C", "DFG": "G", "DFI": "X", "DFO": "X", "DFT": "N",
    "DG ": "G", "DGH": "G", "DGI": "G", "DGL": "E", "DGN": "Q",
    "DHA": "S", "DHI": "H", "DHL": "X", "DHN": "V", "DHP": "X",
    "DHU": "U", "DHV": "V", "DI ": "I", "DIL": "I", "DIR": "R",
    "DIV": "V", "DLE": "L", "DLS": "K", "DLY": "K", "DM0": "K",
    "DMH": "N", "DMK": "D", "DMT": "X", "DN ": "N", "DNE": "L",
    "DNG": "L", "DNL": "K", "DNM": "L", "DNP": "A", "DNR": "C",
    "DNS": "K", "DOA": "X", "DOC": "C", "DOH": "D", "DON": "L",
    "DPB": "T", "DPH": "F", "DPL": "P", "DPP": "A", "DPQ": "Y",
    "DPR": "P", "DPY": "N", "DRM": "U", "DRP": "N", "DRT": "T",
    "DRZ": "N", "DSE": "S", "DSG": "N", "DSN": "S", "DSP": "D",
    "DT ": "T", "DTH": "T", "DTR": "W", "DTY": "Y", "DU ": "U",
    "DVA": "V", "DXD": "N", "DXN": "N", "DYG": "DYG", "DYS": "C",
    "DZM": "A", "E  ": "A", "E1X": "A", "ECC": "Q", "EDA": "A",
    "EFC": "C", "EHP": "F", "EIT": "T", "ENP": "N", "ESB": "Y",
    "ESC": "M", "EXB": "X", "EXY": "L", "EY5": "N", "EYS": "X",
    "F2F": "F", "FA2": "A", "FA5": "N", "FAG": "N", "FAI": "N",
    "FB5": "A", "FB6": "A", "FCL": "F", "FFD": "N", "FGA": "E",
    "FGL": "G", "FGP": "S", "FHL": "X", "FHO": "K", "FHU": "U",
    "FLA": "A", "FLE": "L", "FLT": "Y", "FME": "M", "FMG": "G",
    "FMU": "N", "FOE": "C", "FOX": "G", "FP9": "P", "FPA": "F",
    "FRD": "X", "FT6": "W", "FTR": "W", "FTY": "Y", "FVA": "V",
    "FZN": "K", "G  ": "G", "G25": "G", "G2L": "G", "G2S": "G",
    "G31": "G", "G32": "G", "G33": "G", "G36": "G", "G38": "G",
    "G42": "G", "G46": "G", "G47": "G", "G48": "G", "G49": "G",
    "G4P": "N", "G7M": "G", "GAO": "G", "GAU": "E", "GCK": "C",
    "GCM": "X", "GDP": "G", "GDR": "G", "GFL": "G", "GGL": "E",
    "GH3": "G", "GHG": "Q", "GHP": "G", "GL3": "G", "GLH": "Q",
    "GLJ": "E", "GLK": "E", "GLM": "X", "GLN": "Q", "GLQ": "E",
    "GLU": "E", "GLX": "Z", "GLY": "G", "GLZ": "G", "GMA": "E",
    "GMS": "G", "GMU": "U", "GN7": "G", "GND": "X", "GNE": "N",
    "GOM": "G", "GPL": "K", "GS ": "G", "GSC": "G", "GSR": "G",
    "GSS": "G", "GSU": "E", "GT9": "C", "GTP": "G", "GVL": "X",
    "GYC": "CYG", "GYS": "SYG", "H2U": "U", "H5M": "P", "HAC": "A",
    "HAR": "R", "HBN": "H", "HCS": "X", "HDP": "U", "HEU": "U",
    "HFA": "X", "HGL": "X", "HHI": "H", "HHK": "AK", "HIA": "H",
    "HIC": "H", "HIP": "H", "HIQ": "H", "HIS": "H", "HL2": "L",
    "HLU": "L", "HMR": "R", "HOL": "N", "HPC": "F", "HPE": "F",
    "HPH": "F", "HPQ": "F", "HQA": "A", "HRG": "R", "HRP": "W",
    "HS8": "H", "HS9": "H", "HSE": "S", "HSL": "S", "HSO": "H",
    "HTI": "C", "HTN": "N", "HTR": "W", "HV5": "A", "HVA": "V",
    "HY3": "P", "HYP": "P", "HZP": "P", "I  ": "I", "I2M": "I",
    "I58": "K", "I5C": "C", "IAM": "A", "IAR": "R", "IAS": "D",
    "IC ": "C", "IEL": "K", "IEY": "HYG", "IG ": "G", "IGL": "G",
    "IGU": "G", "IIC": "SHG", "IIL": "I", "ILE": "I", "ILG": "E",
    "ILX": "I", "IMC": "C", "IML": "I", "IOY": "F", "IPG": "G",
    "IPN": "N", "IRN": "N", "IT1": "K", "IU ": "U", "IYR": "Y",
    "IYT": "T", "IZO": "M", "JJJ": "C", "JJK": "C", "JJL": "C",
    "JW5": "N", "K1R": "C", "KAG": "G", "KCX": "K", "KGC": "K",
    "KNB": "A", "KOR": "M", "KPI": "K", "KST": "K", "KYQ": "K",
    "L2A": "X", "LA2": "K", "LAA": "D", "LAL": "A", "LBY": "K",
    "LC ": "C", "LCA": "A", "LCC": "N", "LCG": "G", "LCH": "N",
    "LCK": "K", "LCX": "K", "LDH": "K", "LED": "L", "LEF": "L",
    "LEH": "L", "LEI": "V", "LEM": "L", "LEN": "L", "LET": "X",
    "LEU": "L", "LEX": "L", "LG ": "G", "LGP": "G", "LHC": "X",
    "LHU": "U", "LKC": "N", "LLP": "K", "LLY": "K", "LME": "E",
    "LMF": "K", "LMQ": "Q", "LMS": "N", "LP6": "K", "LPD": "P",
    "LPG": "G", "LPL": "X", "LPS": "S", "LSO": "X", "LTA": "X",
    "LTR": "W", "LVG": "G", "LVN": "V", "LYF": "K", "LYK": "K",
    "LYM": "K", "LYN": "K", "LYR": "K", "LYS": "K", "LYX": "K",
    "LYZ": "K", "M0H": "C", "M1G": "G", "M2G": "G", "M2L": "K",
    "M2S": "M", "M30": "G", "M3L": "K", "M5M": "C", "MA ": "A",
    "MA6": "A", "MA7": "A", "MAA": "A", "MAD": "A", "MAI": "R",
    "MBQ": "Y", "MBZ": "N", "MC1": "S", "MCG": "X", "MCL": "K",
    "MCS": "C", "MCY": "C", "MD3": "C", "MD6": "G", "MDH": "X",
    "MDO": "ASG", "MDR": "N", "MEA": "F", "MED": "M", "MEG": "E",
    "MEN": "N", "MEP": "U", "MEQ": "Q", "MET": "M", "MEU": "G",
    "MF3": "X", "MFC": "GYG", "MG1": "G", "MGG": "R", "MGN": "Q",
    "MGQ": "A", "MGV": "G", "MGY": "G", "MHL": "L", "MHO": "M",
    "MHS": "H", "MIA": "A", "MIS": "S", "MK8": "L", "ML3": "K",
    "MLE": "L", "MLL": "L", "MLY": "K", "MLZ": "K", "MME": "M",
    "MMO": "R", "MMT": "T", "MND": "N", "MNL": "L", "MNU": "U",
    "MNV": "V", "MOD": "X", "MP8": "P", "MPH": "X", "MPJ": "X",
    "MPQ": "G", "MRG": "G", "MSA": "G", "MSE": "M", "MSL": "M",
    "MSO": "M", "MSP": "X", "MT2": "M", "MTR": "T", "MTU": "A",
    "MTY": "Y", "MVA": "V", "N  ": "N", "N10": "S", "N2C": "X",
    "N5I": "N", "N5M": "C", "N6G": "G", "N7P": "P", "NA8": "A",
    "NAL": "A", "NAM": "A", "NB8": "N", "NBQ": "Y", "NC1": "S",
    "NCB": "A", "NCX": "N", "NCY": "X", "NDF": "F", "NDN": "U",
    "NEM": "H", "NEP": "H", "NF2": "N", "NFA": "F", "NHL": "E",
    "NIT": "X", "NIY": "Y", "NLE": "L", "NLN": "L", "NLO": "L",
    "NLP": "L", "NLQ": "Q", "NMC": "G", "NMM": "R", "NMS": "T",
    "NMT": "T", "NNH": "R", "NP3": "N", "NPH": "C", "NPI": "A",
    "NRP": "LYG", "NRQ": "MYG", "NSK": "X", "NTY": "Y", "NVA": "V",
    "NYC": "TWG", "NYG": "NYG", "NYM": "N", "NYS": "C", "NZH": "H",
    "O12": "X", "O2C": "N", "O2G": "G", "OAD": "N", "OAS": "S",
    "OBF": "X", "OBS": "X", "OCS": "C", "OCY": "C", "ODP": "N",
    "OHI": "H", "OHS": "D", "OIC": "X", "OIP": "I", "OLE": "X",
    "OLT": "T", "OLZ": "S", "OMC": "C", "OMG": "G", "OMT": "M",
    "OMU": "U", "ONE": "U", "ONH": "A", "ONL": "X", "OPR": "R",
    "ORN": "A", "ORQ": "R", "OSE": "S", "OTB": "X", "OTH": "T",
    "OTY": "Y", "OXX": "D", "P  ": "G", "P1L": "C", "P1P": "N",
    "P2T": "T", "P2U": "U", "P2Y": "P", "P5P": "A", "PAQ": "Y",
    "PAS": "D", "PAT": "W", "PAU": "A", "PBB": "C", "PBF": "F",
    "PBT": "N", "PCA": "E", "PCC": "P", "PCE": "X", "PCS": "F",
    "PDL": "X", "PDU": "U", "PEC": "C", "PF5": "F", "PFF": "F",
    "PFX": "X", "PG1": "S", "PG7": "G", "PG9": "G", "PGL": "X",
    "PGN": "G", "PGP": "G", "PGY": "G", "PHA": "F", "PHD": "D",
    "PHE": "F", "PHI": "F", "PHL": "F", "PHM": "F", "PIA": "AYG",
    "PIV": "X", "PLE": "L", "PM3": "F", "PMT": "C", "POM": "P",
    "PPN": "F", "PPU": "A", "PPW": "G", "PQ1": "N", "PR3": "C",
    "PR5": "A", "PR9": "P", "PRN": "A", "PRO": "P", "PRS": "P",
    "PSA": "F", "PSH": "H", "PST": "T", "PSU": "U", "PSW": "C",
    "PTA": "X", "PTH": "Y", "PTM": "Y", "PTR": "Y", "PU ": "A",
    "PUY": "N", "PVH": "H", "PVL": "X", "PYA": "A", "PYO": "U",
    "PYX": "C", "PYY": "N", "QLG": "QLG", "QMM": "Q", "QPA": "C",
    "QPH": "F", "QUO": "G", "R  ": "A", "R1A": "C", "R4K": "W",
    "RC7": "HYG", "RE0": "W", "RE3": "W", "RIA": "A", "RMP": "A",
    "RON": "X", "RT ": "T", "RTP": "N", "S1H": "S", "S2C": "C",
    "S2D": "A", "S2M": "T", "S2P": "A", "S4A": "A", "S4C": "C",
    "S4G": "G", "S4U": "U", "S6G": "G", "SAC": "S", "SAH": "C",
    "SAR": "G", "SBL": "S", "SC ": "C", "SCH": "C", "SCS": "C",
    "SCY": "C", "SD2": "X", "SDG": "G", "SDP": "S", "SEB": "S",
    "SEC": "A", "SEG": "A", "SEL": "S", "SEM": "S", "SEN": "S",
    "SEP": "S", "SER": "S", "SET": "S", "SGB": "S", "SHC": "C",
    "SHP": "G", "SHR": "K", "SIB": "C", "SIC": "DC", "SLA": "P",
    "SLR": "P", "SLZ": "K", "SMC": "C", "SME": "M", "SMF": "F",
    "SMP": "A", "SMT": "T", "SNC": "C", "SNN": "N", "SOC": "C",
    "SOS": "N", "SOY": "S", "SPT": "T", "SRA": "A", "SSU": "U",
    "STY": "Y", "SUB": "X", "SUI": "DG", "SUN": "S", "SUR": "U",
    "SVA": "S", "SVV": "S", "SVW": "S", "SVX": "S", "SVY": "S",
    "SVZ": "X", "SWG": "SWG", "SYS": "C", "T  ": "T", "T11": "F",
    "T23": "T", "T2S": "T", "T2T": "N", "T31": "U", "T32": "T",
    "T36": "T", "T37": "T", "T38": "T", "T39": "T", "T3P": "T",
    "T41": "T", "T48": "T", "T49": "T", "T4S": "T", "T5O": "U",
    "T5S": "T", "T66": "X", "T6A": "A", "TA3": "T", "TA4": "X",
    "TAF": "T", "TAL": "N", "TAV": "D", "TBG": "V", "TBM": "T",
    "TC1": "C", "TCP": "T", "TCQ": "Y", "TCR": "W", "TCY": "A",
    "TDD": "L", "TDY": "T", "TFE": "T", "TFO": "A", "TFQ": "F",
    "TFT": "T", "TGP": "G", "TH6": "T", "THC": "T", "THO": "X",
    "THR": "T", "THX": "N", "THZ": "R", "TIH": "A", "TLB": "N",
    "TLC": "T", "TLN": "U", "TMB": "T", "TMD": "T", "TNB": "C",
    "TNR": "S", "TOX": "W", "TP1": "T", "TPC": "C", "TPG": "G",
    "TPH": "X", "TPL": "W", "TPO": "T", "TPQ": "Y", "TQI": "W",
    "TQQ": "W", "TRF": "W", "TRG": "K", "TRN": "W", "TRO": "W",
    "TRP": "W", "TRQ": "W", "TRW": "W", "TRX": "W", "TS ": "N",
    "TST": "X", "TT ": "N", "TTD": "T", "TTI": "U", "TTM": "T",
    "TTQ": "W", "TTS": "Y", "TY1": "Y", "TY2": "Y", "TY3": "Y",
    "TY5": "Y", "TYB": "Y", "TYI": "Y", "TYJ": "Y", "TYN": "Y",
    "TYO": "Y", "TYQ": "Y", "TYR": "Y", "TYS": "Y", "TYT": "Y",
    "TYU": "N", "TYW": "Y", "TYX": "X", "TYY": "Y", "TZB": "X",
    "TZO": "X", "U  ": "U", "U25": "U", "U2L": "U", "U2N": "U",
    "U2P": "U", "U31": "U", "U33": "U", "U34": "U", "U36": "U",
    "U37": "U", "U8U": "U", "UAR": "U", "UCL": "U", "UD5": "U",
    "UDP": "N", "UFP": "N", "UFR": "U", "UFT": "U", "UMA": "A",
    "UMP": "U", "UMS": "U", "UN1": "X", "UN2": "X", "UNK": "X",
    "UR3": "U", "URD": "U", "US1": "U", "US2": "U", "US3": "T",
    "US5": "U", "USM": "U", "VAD": "V", "VAF": "V", "VAL": "V",
    "VB1": "K", "VDL": "X", "VLL": "X", "VLM": "X", "VMS": "X",
    "VOL": "X", "WCR": "GYG", "X  ": "G", "X2W": "E", "X4A": "N",
    "X9Q": "AFG", "XAD": "A", "XAE": "N", "XAL": "A", "XAR": "N",
    "XCL": "C", "XCN": "C", "XCP": "X", "XCR": "C", "XCS": "N",
    "XCT": "C", "XCY": "C", "XGA": "N", "XGL": "G", "XGR": "G",
    "XGU": "G", "XPR": "P", "XSN": "N", "XTH": "T", "XTL": "T",
    "XTR": "T", "XTS": "G", "XTY": "N", "XUA": "A", "XUG": "G",
    "XX1": "K", "XXY": "THG", "XYG": "DYG", "Y  ": "A", "YCM": "C",
    "YG ": "G", "YOF": "Y", "YRR": "N", "YYG": "G", "Z  ": "C",
    "Z01": "A", "ZAD": "A", "ZAL": "A", "ZBC": "C", "ZBU": "U",
    "ZCL": "F", "ZCY": "C", "ZDU": "U", "ZFB": "X", "ZGU": "G",
    "ZHP": "N", "ZTH": "T", "ZU0": "T", "ZZJ": "A"}

oneToThree = {'C':'CYS', 
'D':'ASP',
'S':'SER',
'V':'VAL',
'Q':'GLN',
'K':'LYS',
'P':'PRO',
'T':'THR',
'F':'PHE',
'A':'ALA',
'H':'HIS',
'G':'GLY',
'I':'ILE',
'E':'GLU',
'L':'LEU',
'R':'ARG',
'W':'TRP',
'N':'ASN',
'Y':'TYR',
'M':'MET',
'X':'UNK'}

vdw_radius = {"C":1.7,"O":1.52,"N":1.55,"F":1.47,"P":1.8,"S":1.8,"X":1.7,'A':1.85,'I':1.98,'B':1.85,'R':2.0, 'V':1.53}
#threshs = [3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0]
#threshs = [3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25,8.5,8.75,
#9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.25,11.75,12.0,12.25,12.5,12.75,
#13.0,13.25,13.4,13.75,14.0,14.25,14.5,14.75,15.0,15.25,15.0,15.25,15.5,15.75,16.0,16.25,16.5,16.75]

cone_cut_volumes = {'CYS': 542.8672105403161, 'ILE': 859.5397500221673, 'SER': 592.3734747731354, 'GLN': 859.5397500221673, 'LYS': 662.0644718052689, 'ASN': 662.0644718052689, 'PRO': 662.0644718052689, 'THR': 592.3734747731354, 'PHE': 1092.8291844899893, 'ALA': 487.836979224935, 'HIS': 696.9099703213358, 'GLY': 471.23889803846896, 'ASP': 662.0644718052689, 'SEC': 542.8672105403161, 'LEU': 814.3008158104743, 'ARG': 814.3008158104743, 'TRP': 1364.9172882296452, 'VAL': 814.3008158104743, 'GLU': 662.0644718052689, 'TYR': 1092.8291844899893, 'MET': 859.5397500221673}

threshs = numpy.linspace(2,14,25)
angle_range = numpy.linspace(-1.0,0.5,16)
#threshs = angle_range
max_sphere = 40.0
radii_map = {'CYS': 3.029584442714009, 'SEC': 3.029584442714009,'ILE': 3.3236043438778844, 'SER': 2.942862799880502, 'GLN': 3.392265173640708, 'LYS': 3.4323769248759275, 'TRP': 4.020778501571494, 'PRO': 3.1681967871564782, 'THR': 3.1209640826827605, 'PHE': 3.7193696421024685, 'ALA': 2.8009640812798096, 'GLY': 2.5743877263600927, 'HIS': 3.534660601593296, 'ASN': 3.2435252870320834, 'LEU': 3.3236043438778844, 'ARG': 3.631344304999469, 'ASP': 3.236792124221129, 'VAL': 3.1681967871564782, 'GLU': 3.3861111330098956, 'TYR': 3.8021338661518342, 'MET': 3.351107764635289}

avg_cent_side = {'ILE': (-0.90613294951389123, -1.2696649328274769, 0.75787678394555402), 'GLN': (-1.2338776755913703, -1.6271318863182545, 0.98235396198659397), 'GLY': (-0.0047299661865784631, -0.0088845691075302401, 9.8567539791531887e-05), 'GLU': (-1.2149070583116599, -1.5273224890561412, 1.0878295037481032), 'CYS': (-0.58968190337597537, -0.89958163307781258, 0.62776825003201231), 'HIS': (-0.9776890826599115, -1.6575223413455686, 0.78038430859923047), 'SER': (-0.59170713881015236, -0.65044117400421442, 0.79631448621746626), 'LYS': (-1.4201311703581181, -1.8190188812624906, 1.0949225828431293), 'PRO': (0.21186055970214385, -0.81932530473289122, 1.0393498184079177), 'SEC': (-0.83941591976390406, -0.71433578942505138, 0.85735910033472262), 'ASN': (-0.84857508137952831, -1.3074529598975138, 0.7776964830504921), 'VAL': (-0.79991957561390237, -0.94684856115170413, 0.62104228622222646), 'THR': (-0.65195583478889552, -0.90724937534582129, 0.82850637927897564), 'ASP': (-0.90583763237409776, -1.1848068294814498, 0.82225350991223289), 'TRP': (-1.528171801920162, -1.7410619994461729, 0.89554186538385805), 'UNK': (-0.3674831972563582, -0.32646621880413212, 0.49726141188550765), 'PHE': (-1.1835851942281017, -1.7644149833489573, 0.692845774001975), 'ALA': (-0.40699060927186675, -0.36885898636123338, 0.49092580988001461), 'MET': (-1.1860699223797722, -1.5396420763646954, 0.84853591931547168), 'LEU': (-1.0557131262899473, -1.4510720809764079, 0.67350078632184163), 'ARG': (-1.7179619489063336, -2.0937800414713896, 1.3885442733060744), 'TYR': (-1.2761868751958905, -1.9287097734708729, 0.69640859762760676)}

unknown_avg_cent_side = (-0.866137516594, -1.16682367833, 0.769774990962) #average vector over all residue types
unknown_thresh = 7.5
unknown_angle_thresh = -1.0
unknown_rad = 3.23342916549


#Residue distribution:  {'CYS': 0.012580215522460346, 'ASP': 0.0612362271461436, 'SER': 0.059280784598619685, 'VAL': 0.07643782540259111, 'GLN': 0.034986075796101225, 'LYS': 0.05867538442910764, 'ILE': 0.06109093110546071, 'PRO': 0.04416999636759898, 'THR': 0.05295435282721879, 'PHE': 0.04113694151834363, 'ALA': 0.08321225329943092, 'GLY': 0.07504540501271341, 'HIS': 0.022236348226177503, 'GLU': 0.06778060297856883, 'LEU': 0.09070710739799007, 'ARG': 0.05026032207289018, 'TRP': 0.011520765225814264, 'ASN': 0.042135851798038505, 'TYR': 0.03321225329943092, 'MET': 0.021340355975299673}

avg_cent = {'ILE': (-0.51939923132192789, -0.42924346469309471, 0.57478473070875846), 'GLN': (-0.80834986146513732, -0.76972077063003586, 0.68335543707029034), 'GLY': (0.003710108790125123, 0.74850253699086222, -0.011107678632392972), 'GLU': (-0.78586343234711975, -0.68596611605471958, 0.72074313020193725), 'CYS': (-0.25453284884232952, 0.025290152329503414, 0.46347585732429047), 'HIS': (-0.67516717655691016, -0.87046086894289898, 0.61423651641973831), 'SER': (-0.25186295200626596, 0.15435437507216301, 0.50821742795023583), 'LYS': (-0.93157547397808582, -0.90270366121918721, 0.7564363100362167), 'PRO': (0.19092956143389894, -0.054932017506649618, 0.7190732935459736), 'SEC': (-0.4041486429034104, 0.19276231247130715, 0.31022206884216508), 'ASN': (-0.52337772135489558, -0.46464652365549852, 0.5839222978918297), 'VAL': (-0.40625381747595429, -0.12555581649592371, 0.49308718269365215), 'THR': (-0.3458279994581282, -0.10156321384665247, 0.56845475732264672), 'ASP': (-0.5395255947207902, -0.37285940531331668, 0.57492539764983719), 'TRP': (-1.1801359800054108, -1.1546490468404469, 0.74010431975141766), 'UNK': (-0.16216122024915333, 0.48606836524146568, 0.13849611887759408), 'PHE': (-0.83581175433622168, -1.0098420145081217, 0.57161519965229146), 'ALA': (-0.11519817930399062, 0.43429992078978386, 0.21401964289355344), 'MET': (-0.70928160313770094, -0.59184549532823161, 0.58483387893500882), 'LEU': (-0.62691092485685773, -0.54401570724965409, 0.4676475183702658), 'ARG': (-1.2457704954369544, -1.2737121830387164, 1.037835446531018), 'TYR': (-0.9381954559286525, -1.1994693553734199, 0.57799294311076543)}

average_correction = {'CYS': -0.022946150704117624, 'SEC': -0.022946150704117624, 'GLN': -0.02961937816659793, 'ASP': -0.013910148202608552, 'SER': 0.013356648242644154, 'VAL': -0.03795902095389797, 'LYS': -0.03942889003072603, 'ILE': -0.05977847074136369, 'PRO': -0.004513042646024423, 'THR': -0.009849324564857875, 'PHE': -0.04798429467867554, 'ALA': 0.009567509980448752, 'GLY': 0.01743352671606324, 'HIS': -0.022794305142026405, 'GLU': -0.044020993580337175, 'LEU': -0.0641172463977473, 'ARG': -0.0848773838260547, 'TRP': -0.10482640485221231, 'ASN': 0.012512801252194339, 'TYR': -0.04262395441484643, 'MET': -0.059906004162288945}

average_correction_alpha = {'CYS': -0.08178375570565236, 'SEC': -0.08178375570565236, 'ILE': -0.0763942769641649, 'SER': -0.056276272825953266, 'GLN': -0.05432425352766346, 'LYS': -0.04347620974426458, 'TRP': -0.13120029100479658, 'PRO': -0.02893434811025425, 'THR': -0.06906454564366421, 'PHE': -0.08442494055881009, 'ALA': -0.05209043431185281, 'GLY': -0.029668853492853364, 'HIS': -0.07807760374070327, 'ASN': -0.04205405440700721, 'LEU': -0.11715899801205235, 'ARG': -0.08525007619343955, 'ASP': -0.05334347102978171, 'VAL': -0.09920832263668614, 'GLU': -0.05777287475315668, 'TYR': -0.08999231904470734, 'MET': -0.10386081073354511}


#values learned on old training set
#threshs = {'CYS' : 7.25, 'SEC': 7.25, 'GLN' : 7.25, 'ASP' : 6.75, 'SER' : 6.5, 'VAL' : 7.5, 'LYS' : 7.0, 'ASN' : 7.0, 'PRO' : 7.25, 'THR' : 7.0, 'PHE' : 8.25, 'ALA' : 6.75, 'HIS' : 7.25, 'GLY' : 6.75, 'ILE' : 7.5 , 'LEU' : 7.75, 'ARG' : 7.75, 'TRP' : 8.25, 'GLU' : 7.0, 'TYR' : 8.0, 'MET' : 7.5}

#angle_threshs = {'CYS' : -0.85, 'SEC' : -0.85, 'GLN' : -1.0, 'ASP' : -0.95, 'SER' : -0.9, 'VAL' : -0.85, 'LYS' : -1.0, 'ASN' : -0.9, 'PRO' : -0.9, 'THR' : -0.9, 'PHE' : -0.8, 'ALA' : -0.85, 'HIS' : -0.75, 'GLY' : -1.0, 'ILE' : -0.9, 'LEU' : -0.9, 'ARG' : -1.0, 'TRP' : -0.95, 'GLU' : -1.0, 'TYR' : -0.75, 'MET' : -0.95}

#threshs_alpha = {'CYS' : 7.25, 'SEC' : 7.25,'ILE' : 7.75, 'SER' : 6.75, 'GLN' : 8.0, 'LYS' : 8.0, 'PRO' : 7.25, 'ASP' : 7.0, 'THR' : 7.0, 'PHE' : 8.5, 'ALA' : 7.0, 'GLY' : 6.75, 'HIS' : 8.0, 'GLU' : 7.75, 'LEU' : 7.75, 'ARG' : 8.75, 'TRP' : 8.75, 'VAL' : 7.5, 'ASN' : 7.25, 'TYR' : 8.75, 'MET' : 8.0}

#angle_threshs_alpha = {'CYS' : -0.9, 'SEC' : -0.9, 'GLN' : -1.0, 'ASP' : -1.0, 'SER' : -1.0, 'VAL' : -1.0, 'LYS' : -1.0, 'ASN' : -1.0, 'PRO' : -1.0, 'THR' : -1.0, 'PHE' : -0.85, 'ALA' : -0.9, 'HIS' : -1.0, 'GLY' : -1.0, 'ILE' : -0.85, 'LEU' : -1.0, 'ARG' : -1.0, 'TRP' : -1.0, 'GLU' : -1.0, 'TYR' : -0.95, 'MET' : -1.0}

#learned on golden train set without included target residue
#"""
threshs_alpha ={
'CYS' : 7.25, 
'ASP' : 7.25,
'SER' : 7.0,
'VAL' : 7.5,
'GLN' : 7.75,
'LYS' : 8.0,
'PRO' : 7.5,
'THR' : 7.0,
'PHE' : 8.25,
'ALA' : 7.0,
'HIS' : 8.0,
'GLY' : 6.75,
'ILE' : 7.75,
'GLU' : 7.5,
'LEU' : 7.75,
'ARG' : 8.5,
'TRP' : 9.0,
'ASN' : 7.25,
'TYR' : 8.75,
'MET' : 8.0
}

angle_threshs_alpha = {
'CYS' : -0.85,
'ASP' : -0.90,
'SER' : -0.85,
'VAL' : -0.90,
'GLN' : -1.0,
'LYS' : -1.0,
'PRO' : -1.0,
'THR' : -1.0,
'PHE' : -0.90,
'ALA' : -0.85,
'HIS' : -0.95,
'GLY' : -1.0,
'ILE' : -0.85,
'GLU' : -1.0,
'LEU' : -1.0,
'ARG' : -1.0,
'TRP' : -0.6,
'ASN' : -1.0,
'TYR' : -0.85,
'MET' : -1.0
}
#"""
"""
#learned on golden train with including intersecting residue spheres
threshs_alpha ={
'CYS' : 8.5,
'ASP' : 9.0,
'SER' : 9.0,
'VAL' : 9.0,
'GLN' : 9.5,
'LYS' : 9.5,
'PRO' : 9.0,
'THR' : 9.0,
'PHE' : 10.0,
'ALA' : 9.0,
'HIS' : 9.5,
'GLY' : 8.5,
'ILE' : 10.0,
'GLU' : 9.5,
'LEU' : 9.5,
'ARG' : 10.0,
'TRP' : 10.5,
'ASN' : 9.0,
'TYR' : 10.5,
'MET' : 9.5
}
angle_threshs_alpha = {
'CYS' : -0.9,
'ASP' : -0.9,
'SER' : -0.8,
'VAL' : -0.9,
'GLN' : -1.0,
'LYS' : -0.9,
'PRO' : -1.0,
'THR' : -0.9,
'PHE' : -0.4,
'ALA' : -0.8,
'HIS' : -0.9,
'GLY' : -1.0,
'ILE' : -0.8,
'GLU' : -0.9,
'LEU' : -0.9,
'ARG' : -0.9,
'TRP' : -0.4,
'ASN' : -0.9,
'TYR' : -0.4,
'MET' : -1.0
}
"""
"""
#learned on golden train set with target residue included
threshs_alpha ={
'CYS' : 7.25,
'ASP' : 7.25,
'SER' : 7.0,
'VAL' : 7.5,
'GLN' : 7.75,
'LYS' : 8.0,
'PRO' : 7.5,
'THR' : 7.0,
'PHE' : 8.25,
'ALA' : 7.0,
'HIS' : 8.0,
'GLY' : 6.75,
'ILE' : 7.75,
'GLU' : 7.5,
'LEU' : 7.75,
'ARG' : 8.5,
'TRP' : 9.0,
'ASN' : 7.25,
'TYR' : 8.75,
'MET' : 8.0}

angle_threshs_alpha = {
'CYS' : -0.85,
'ASP' : -0.90,
'SER' : -0.85,
'VAL' : -0.90,
'GLN' : -1.0,
'LYS' : -1.0,
'PRO' : -1.0,
'THR' : -1.0,
'PHE' : -0.90,
'ALA' : -0.85,
'HIS' : -0.95,
'GLY' : -1.0,
'ILE' : -0.85,
'GLU' : -1.0,
'LEU' : -1.0,
'ARG' : -1.0,
'TRP' : -0.60,
'ASN' : -1.0,
'TYR' : -0.85,
'MET' : -1.0
}
"""
#Best threshs for c_alpha and unknown target residue: 7.5,-1.0

#learned values from golden train set
threshs = {'CYS' : 7.0, 'SEC': 7.0, 'GLN' : 7.25, 'ASP' : 6.75, 'ASX': 6.75, 'SER' : 6.5, 'VAL' : 7.5, 'LYS' : 7.5, 'ASN' : 6.75, 'PRO' : 7.0, 'THR' : 6.75, 'PHE' : 7.75, 'ALA' : 7.0, 'HIS' : 7.25, 'GLY' : 6.75, 'ILE' : 7.75 , 'LEU' : 7.75, 'ARG' : 7.75, 'TRP' : 8.0, 'GLU' : 7.0, 'GLX': 7.0, 'TYR' : 7.75, 'MET' : 7.5, 'UNK': unknown_thresh}

angle_threshs = {'CYS' : -0.85, 'SEC' : -0.85, 'GLN' : -0.95, 'ASP' : -0.95, 'ASX' : -0.95, 'SER' : -0.9, 'VAL' : -0.85, 'LYS' : -0.9, 'ASN' : -0.95, 'PRO' : -0.85, 'THR' : -0.85, 'PHE' : -0.85, 'ALA' : -0.85, 'HIS' : -0.8, 'GLY' : -1.0, 'ILE' : -0.9, 'LEU' : -0.9, 'ARG' : -0.95, 'TRP' : -0.95, 'GLU' : -0.95, 'GLX' : -0.95, 'TYR' : -0.85, 'MET' : -1.0, 'UNK': unknown_angle_thresh}

residue_max_acc = { 
# Miller max acc: Miller et al. 1987 http://dx.doi.org/10.1016/0022-2836(87)90038-6 
# Wilke: Tien et al. 2013 http://dx.doi.org/10.1371/journal.pone.0080635 
# Sander: Sander & Rost 1994 http://dx.doi.org/10.1002/prot.340200303 
    'Miller': { 
        'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0, 
        'CYS': 140.0, 'GLN': 189.0, 'GLU': 183.0, 'GLY': 85.0, 
        'HIS': 194.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 211.0, 
        'MET': 204.0, 'PHE': 218.0, 'PRO': 143.0, 'SER': 122.0, 
        'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0 
    }, 
    'Wilke': { 
        'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 
        'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 
        'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 
        'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 
        'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0 
    }, 
    'Sander': { 
        'ALA': 106.0, 'ARG': 248.0, 'ASN': 157.0, 'ASP': 163.0, 
        'CYS': 135.0, 'GLN': 198.0, 'GLU': 194.0, 'GLY': 84.0, 
        'HIS': 184.0, 'ILE': 169.0, 'LEU': 164.0, 'LYS': 205.0, 
        'MET': 188.0, 'PHE': 197.0, 'PRO': 136.0, 'SER': 130.0, 
        'THR': 142.0, 'TRP': 227.0, 'TYR': 222.0, 'VAL': 142.0,
        'SEC': 135.0 
    } 
} 

def calcVol(r,cos):
    vol = (2.0/3.0)*math.pi*(r**3.0)*(1.0-cos)
    return vol

def sissCorrectionVol(siss,r,cos):
    v = calcVol(r,cos)
    return (v-siss)/v

def sissCorrection(siss,res):
    v = cone_cut_volumes[res]
    max_acc = residue_max_acc['Sander'][res]
    #corrected_siss = siss/(v*max_acc)
    #corrected_siss = siss/(v)
    corrected_siss = (v-siss)/v
    #corrected_siss = siss/(max_acc)
    #corrected_siss = siss/(v+max_acc)
    return corrected_siss

def parsePDB(input_file,chain,c_alpha):
    """
    Parses a PDB-file and takes all atomic coordinates of a specified chain.

    Input:
    input_file: String ; Path to a PDB file
    chain: String or None ; Chain identifier, if None is given, the first Chain found in the file is taken
    c_alpha: Boolean ; If True, then only C alpha atoms are taken

    Output:
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.   
    """
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()
    coordinate_map = {}

    x_total = 0.0
    y_total = 0.0
    z_total = 0.0
    n = 0.0

    #print lines

    for line in lines:
        if len(line) > 5:        
            record_name = line[0:6].replace(" ","")
            if record_name == "ENDMDL":
                #print len(coordinate_map)
                break
        #ignore short lines
        if len(line) > 20:
            atom_nr = line[6:11].replace(" ","")
            atom_name = line[12:16].replace(" ","")
            res_name = line[17:20].replace(" ","")
            
            #len(res_name) == 3 is only true for amino acid chains
            if len(res_name) == 3:
                if len(line) > 21:
                    chain_id = line[21]
                    res_nr = int(line[22:27].replace(" ",""))
                    #insertion_code = line[26]
                
                #consider only lines with record name ATOM
                if record_name == "ATOM" or record_name == 'HETATM':
                    if record_name == 'HETATM':
                        if res_name not in threeToOne:
                            continue
                        else:
                            res_name = OneToThree[threeToOne[res_name][0]]
                    #if chain not given, take first one found
                    if chain == None:
                        chain = chain_id
                    if len(line) > 50:
                        #print chain,chain_id
                        #consider only lines with the correct chain id
                        if chain_id == chain:
                            #if c_alpha is true, then take only C alpha atoms
                            if res_name != 'UNK':
                                if (not c_alpha) or atom_name == 'CA':
                                    if atom_name[0] != 'H' and atom_name[0] != 'D':
                                        x = float(line[30:38].replace(" ",""))
                                        y = float(line[38:46].replace(" ",""))
                                        z = float(line[46:54].replace(" ",""))
                                        if res_nr not in coordinate_map:
                                            coordinate_map[res_nr]=[res_name,{}]
                                        coordinate_map[res_nr][1][atom_nr] = (atom_name,x,y,z)
                                        x_total += x
                                        y_total += y
                                        z_total += z
                                        n += 1.0
    if n > 0.0:                                        
        protein_centroid = numpy.array([x_total/n,y_total/n,z_total/n])
    else:
        protein_centroid = numpy.array([0.0,0.0,0.0])
    #print coordinate_map
    return coordinate_map,protein_centroid

def calcCentroidMap(coordinate_map,target_residues,c_alpha,double_unknown_mode = False):
    #print coordinate_map
    centroid_map = {}
    if c_alpha:
        c1 = None
        c2 = None
        res_2 = None
        res_1 = None
        res_name_1 = None
        res_name_2 = None
        for res in coordinate_map:
            c0 = c1
            c1 = c2
            res_0 = res_1
            res_1 = res_2
            res_2 = res
            res_name_1 = res_name_2
            res_name_2 = coordinate_map[res_2][0]
            atomlist = coordinate_map[res_2][1]
            #if len(atomlist) > 1:
            #    print atomlist
            #    raise NameError('More than one atom in C alpha only case.')
            (atomname,x,y,z) = atomlist.values()[0]
            if atomname != ('CA'):
                raise NameError('Only atom is not C alpha')
            c2 = numpy.array([x,y,z])           
            if res_0 != None and res_1 != None:
                if res_2 - res_1 == 1 and res_1 - res_0 == 1:
                    avg_cent_side_vec = avg_cent_side[res_name_1]
                    avg_cent_vec = avg_cent[res_name_1]
                    if double_unknown_mode: #in double unknown mode, the residue types are not needed!
                        avg_cent_side_vec = unknown_avg_cent_side
                    side_centroid = predict_centroid(c0,c1,c2,avg_cent_side_vec)
                    centroid = predict_centroid(c0,c1,c2,avg_cent_vec)
                    
                    centroid_map[res_1] = (side_centroid,centroid)

                else:
                    centroid_map[res_1] = (c1,c1)
            elif res_1 != None:
                centroid_map[res_1] = (c1,c1)
        if len(coordinate_map) > 0:
            centroid_map[res_2] = (c2,c2)
    else:
        for res in target_residues:
            atomlist = coordinate_map[res][1]
            centroid = getCentroid(atomlist)
            centroid_map[res] = centroid

    return centroid_map

"""
def calcDistMatrix(coordinate_map,target_residues):
    
    Calculates distances matrix for all atoms in a set of target residues and all atoms in a coordinate map. Distances between atoms of the same residue are not computed
    
    Input:
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.   
    target_residues: [String] or None; contains the residue-ids for all residues, for which the distances are needed. If None is given, compute all distances inbetween the coordinate_map.
    
    Output:
    dist_matrix: {String:{String:Float}} ; Matrix of [Atom-ID][Atom-ID] -> euclidean distance. IDs are not ordered!
    
    dist_matrix= {}

    #If None is given, take all residues from the coordinate_map
    if target_residues == None:
        target_residues = coordinate_map.keys()
    for res_id_1 in target_residues:
        for atom_id_1 in coordinate_map[res_id_1][1]:
            (atom_name,x,y,z) = coordinate_map[res_id_1][1][atom_id_1]
            vec1 = numpy.array([x,y,z])
            if not atom_id_1 in dist_matrix:
                dist_matrix[atom_id_1 = {}
            for res_id_2 in coordinate_map:
                #ignore if same residue
                if res_id_2 != res_id_1:
                    for atom_id_2 in coordinate_map[res_id_2][1]:
                        #check if distance already computed
                        if not atom_id_2 in dist_matrix[atom_id_1]:
                            if not atom_id_2 in dist_matrix:
                                (atom_name,x,y,z) = coordinate_map[res_id_2][1][atom_id_2]
                                vec2 = numpy.array([x,y,z])
                                diff = vec2 - vec1
                                d = numpy.sqrt(numpy.dot(diff, diff))
                                dist_matrix[atom_id_1][atom_id_2] = d
                            elif not atom_id_1 in dist_matrix[atom_id_2]:
                                (atom_name,x,y,z) = coordinate_map[res_id_2][1][atom_id_2]
                                vec2 = numpy.array([x,y,z])
                                diff = vec2 - vec1
                                d = numpy.sqrt(numpy.dot(diff, diff))
                                dist_matrix[atom_id_1][atom_id_2] = d

    return dist_matrix 
"""

def calcDistMatrix(coordinate_map,centroid_map,target_residues,c_alpha):
    dist_matrix = {}
    if c_alpha:
        for res_1 in centroid_map:
            #dist_matrix[res_1] = {}
            centroid_1 = centroid_map[res_1][0]
            for res_2 in centroid_map:
                if not (res_1,res_2) in dist_matrix:
                    centroid_2 = centroid_map[res_2][0]
                    #if centroid_2 == None:
                    #    print centroid_map
                    diff = centroid_2 - centroid_1
                    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                    dist_matrix[(res_1,res_2)] = d
                    dist_matrix[(res_2,res_1)] = d
    else:
        for res_1 in target_residues:
            centroid = centroid_map[res_1]
            res_name = coordinate_map[res_1][0]
            thresh = threshs[res_name]
            dist_matrix[res_1] = {}
            for res_2 in coordinate_map:
                #if res_2 != res_1: #if this is commented, the atoms of the target residue are included in the measure calculation
                dist_matrix[res_1][res_2] = {}
                atomlist = coordinate_map[res_2][1]
                for atom in atomlist:
                    (atomname,x,y,z) = atomlist[atom]
                    diff = centroid - numpy.array([x,y,z])
                    #d = numpy.sqrt(numpy.dot(diff, diff))
                    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                    #If one atom is further than 10A+threshhold, then ignore the whole residue
                    if d - 10.0 > thresh:
                        break

                    dist_matrix[res_1][res_2][atom] = (d,atomname,x,y,z)
    return dist_matrix

def get_gly_cb_vector(residue): 
    """ 
    Return a pseudo CB vector for a Gly residue. 
    The pseudoCB vector is centered at the origin. 

    CB coord=N coord rotated over -120 degrees 
    along the CA-C axis. 
    """ 
    try: 
      n_v = residue["N"].get_vector() 
      c_v = residue["C"].get_vector() 
      ca_v = residue["CA"].get_vector() 
    except: 
      return None 
    # center at origin 
    n_v = n_v - ca_v 
    c_v = c_v - ca_v 
    # rotation around c-ca over -120 deg 
    rot = rotaxis(-math.pi * 120.0 / 180.0, c_v) 
    cb_at_origin_v = n_v.left_multiply(rot)
    vec = numpy.array([cb_at_origin_v[0],cb_at_origin_v[1],cb_at_origin_v[2]])
    return vec


def sphere_intersection(R,r,d):
    #print R,r,d,calcVol(R,-1.0),calcVol(r,-1.0)
    if R < r:
        a = r
        r = R
        R = a
    if r+R <= d:
        return 0.0
    if d+r <= R:
        return (4.0/3.0)*math.pi*r**3.0
    sum1 = (R+r-d)**2
    sum2 = (d**2.0+2.0*d*r-3.0*r**2.0+2.0*d*R+6.0*r*R-3.0*R**2.0)
    si = math.pi*sum1*sum2/(12.0*d)
    #print si
    return si

def createRotMatrix(axis,cos):
    sin = (1.0-cos**2.0)**(1.0/2.0)
    if axis == "x":
        Rot = [[1.0,0.0,0.0],[0.0,cos,-sin],[0.0,sin,cos]]
    if axis == "y":
        Rot = [[cos,0.0,sin],[0.0,1.0,0.0],[-sin,0.0,cos]]
    if axis == "z":
        Rot = [[cos,-sin,0.0],[sin,cos,0.0],[0.0,0.0,1.0]]
    return numpy.matrix(Rot)

def createPlaneRotMatrix(plane,vec):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    if plane == 'xy':
        if z == 0.0:
            cos = 1.0
        elif y == 0.0:
            cos = 0.0
        else:
            cos = (((z**2.0/y**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('x',cos)
        if (vec*rot).A1[2]**2 > 0.0000001:
            rot = rot.T
    if plane == 'xz':
        if y == 0.0:
            cos = 1.0
        elif x == 0.0:
            cos = 0.0
        else:
            cos = (((y**2.0/x**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('z',cos)
        if (vec*rot).A1[1]**2 > 0.0000001:
            rot = rot.T
    if plane == 'yz':
        if x == 0.0:
            cos = 1.0
        elif z == 0.0:
            cos = 0.0
        else:
            cos = (((x**2.0/z**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('y',cos)
        if (vec*rot).A1[0]**2 > 0.0000001:
            rot = rot.T
    return rot

def createAxisRotMatrix(axis,vec):
    if axis == 'x':
        rot1 = createPlaneRotMatrix('xz',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('y',getCos('x',ivec))
        if (ivec*rot2).A1[2]**2 > 0.0000001:
            rot2 = rot2.T
    if axis == 'y':
        rot1 = createPlaneRotMatrix('xy',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('z',getCos('y',ivec))
        if (ivec*rot2).A1[0]**2 > 0.0000001:
            rot2 = rot2.T
    if axis == 'z':
        rot1 = createPlaneRotMatrix('yz',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('x',getCos('z',ivec))
        if (ivec*rot2).A1[1]**2 > 0.0000001:
            rot2 = rot2.T
    #print vec,'vec'
    #print (vec*rot1.T).A1
    #print ivec,'ivec'
    return rot1*rot2

def createRotAxisMatrix(axis,cos):
    sin = (1.0-cos**2.0)**(1.0/2.0)
    axis = axis/numpy.linalg.norm(axis)
    x = axis[0]
    y = axis[1]
    z = axis[2]
    m1 = [cos+x**2.0*(1.0-cos),x*y*(1.0-cos)-z*sin,x*z*(1.0-cos)+y*sin]
    m2 = [y*x*(1.0-cos)+z*sin,cos+y**2.0*(1.0-cos),y*z*(1.0-cos)-x*sin]
    m3 = [z*x*(1.0-cos)-y*sin,z*y*(1.0-cos)+x*sin,cos+z**2.0*(1.0-cos)]
    Rot = [m1,m2,m3]
    return numpy.matrix(Rot)

def gly_vector(n_v,c_v,ca_v):
    n_v = n_v - ca_v 
    c_v = c_v - ca_v
    rot = createRotAxisMatrix(c_v,-0.5)
    vec = (n_v*rot).A1
    return vec
    
def getCosAngle(vec1,vec2):
    n1 = numpy.linalg.norm(vec1)
    n2 = numpy.linalg.norm(vec2)
    norm = n1*n2
    dot = numpy.dot(vec1,vec2)
    if norm != 0.0:
        c = dot / norm
    else:
        return None
    #print c 
    # Take care of roundoff errors 
    c = min(c, 1.0) 
    c = max(-1.0, c)
    return c

def getCos(axis,vec):
    if axis == "x":
        other = numpy.array([1.0,0.0,0.0])
    if axis == "y":
        other = numpy.array([0.0,1.0,0.0])
    if axis == "z":
        other = numpy.array([0.0,0.0,1.0])
    n1 = numpy.linalg.norm(vec)
    n2 = numpy.linalg.norm(other)
    #print n1,n2
    c = ((numpy.dot(vec,other)) / (n1 * n2))
    #print c 
    # Take care of roundoff errors 
    c = min(c, 1.0) 
    c = max(-1.0, c)
    return c

def nullTest(vec):
    if vec[0] == 0.0 and vec[1] == 0.0 and vec[2] == 0.0:
        return False
    return True

def predict_centroid(c0,c1,c2,avg_cent_vec):
    if nullTest(c0-c1) and nullTest(c2-c1):
        A = (c0-c1)/numpy.linalg.norm(c0-c1)
        B = (c2-c1)/numpy.linalg.norm(c2-c1)
        rx = createAxisRotMatrix("x",A)
        B_prime = (B*rx).A1
        r2 = createPlaneRotMatrix('xy',B_prime)
        ROT = rx*r2
        support = (B_prime*r2).A1
        if support[1] < 0.0:
            flip = createRotMatrix('x',-1.0)
            ROT = ROT*flip

        return ((avg_cent_vec*ROT.T).A1)+c1
    else:
        return c1

def getCentroid(atomlist):
    n = 0.0
    t_x = 0.0
    t_y = 0.0
    t_z = 0.0
    c_a = None
    c_b = None
    #print atomlist
    """
    #This for-loop can be used, if C_alpha should be used as Sphere Center
    for atom in atomlist:
        (atomname,x,y,z) = atomlist[atom]
        if atomname == 'CA':
            return numpy.array([x,y,z])
    """
    #"""
    #If this part is commented out, the centroid is calculated by including the backbone
    for atom in atomlist:
        (atomname,x,y,z) = atomlist[atom]
        if len(atomname) > 1:
            t_x += x
            t_y += y
            t_z += z
            n += 1.0
            
    if n > 0.0:
        centroid = numpy.array([t_x/n,t_y/n,t_z/n])
        
    else:
    #"""
        for atom in atomlist:
            (atomname,x,y,z) = atomlist[atom]        
            t_x += x
            t_y += y
            t_z += z
            n += 1.0
            
        if n > 0.0:
            centroid = numpy.array([t_x/n,t_y/n,t_z/n])
        else:
            #Does this happen?
            print "This does happen"
            print atomlist
            centroid = None
    return centroid

def produceOutput(siss_map,coordinate_map,output_file):
    """
    Input:
    siss_map: {String:float} ; 
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.
    """ 

    lines = []    

    for res in siss_map:
        res_name = coordinate_map[res][0]
        siss = siss_map[res]
        lines.append("%s %s\t%s" % (res,res_name,str(siss)))

    lines = sorted(lines,key=lambda x:int(x.split()[0]))

    lines = ["Residue\tSISS value"] + lines

    f = open(output_file,'w')
    f.write("\n".join(lines))
    f.close()

def calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha,manu_thresh=None,manu_angle_thresh=None,double_unknown_mode=False,manu_parameters=None,unknown_parameters=None):
    global threshs
    global angle_threshs
    global threshs_alpha
    global angle_threshs_alpha
    global unknown_thresh
    global unknown_angle_thresh
    siss_map = {}

    if unknown_parameters != None:
        (unknown_thresh,unknown_angle_thresh) = unknown_parameters
    if manu_parameters != None:
        [threshs,angle_threshs] = manu_parameters
        [threshs_alpha,angle_threshs_alpha] = manu_parameters
    if c_alpha:
        for res in target_residues:
            res_name = coordinate_map[res][0]

            if manu_thresh == None:
                thresh = threshs_alpha[res_name]
            else:
                thresh = manu_thresh
            if manu_angle_thresh == None:
                angle_thresh = angle_threshs_alpha[res_name]
            else:
                angle_thresh = manu_angle_thresh

            if thresh == None or double_unknown_mode:
                thresh = unknown_thresh
            if angle_thresh == None or double_unknown_mode:
                angle_thresh = unknown_angle_thresh

            (atomname,x,y,z) = coordinate_map[res][1].values()[0]
            c_alpha_1 = numpy.array([x,y,z])
            centroid_1 = centroid_map[res][0]

            
            #test: substract ssi of all residue spheres from the total sum
            #res_list = []

            siss = 0.0
            for res_2 in coordinate_map:
                #if res == res_2:
                #    siss += sphere_intersection(thresh,radii_map[res_name],0.0)
                #elif res != res_2:
                if res != res_2:
                    res_name_2 = coordinate_map[res_2][0]
                    if (res,res_2) in dist_matrix:
                        d = dist_matrix[(res,res_2)]
                    else:
                        raise NameError('Did not find residue pair in distance matrix for: %s and %s' % (res,res_2))

                    rad = radii_map[res_name_2]
                    if double_unknown_mode:
                        rad = unknown_rad
                    if d - rad <= thresh:
                        if angle_thresh > -1.0: 
                            centroid_2 = centroid_map[res_2][0]
                            
                            angle = getCosAngle(centroid_1-c_alpha_1,centroid_2-c_alpha_1)
                            #print centroid_1,centroid_2,c_alpha_1,angle
                            if angle == None:    
                                siss += sphere_intersection(thresh,rad,d)
                                #res_list.append(res_2)
                            elif angle_thresh <= angle:
                                siss += sphere_intersection(thresh,rad,d)
                                #res_list.append(res_2)
                        else:
                            siss += sphere_intersection(thresh,rad,d)
                            #res_list.append(res_2)

            """
            #test: substract ssi of all residue spheres from the total sum
            done = set([])
            for r1 in res_list:
                done.add(r1)
                for r2 in res_list:
                    if r2 in done:
                        continue
                    if (r1,r2) in dist_matrix:
                        d = dist_matrix[(r1,r2)]
                    else:
                        centroid_1 = centroid_map[r1][0]
                        centroid_2 = centroid_map[r2][0]
                        diff = centroid_2 - centroid_1
                        d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                    siss -= sphere_intersection(radii_map[coordinate_map[r1][0]],radii_map[coordinate_map[r2][0]],d)
            """

            siss = sissCorrectionVol(siss,thresh,angle_thresh)
            #siss_map[res] = siss + average_correction_alpha[res_name]
            siss_map[res] = siss
    else:
        for res in target_residues:
            res_name = coordinate_map[res][0]
            atomlist = coordinate_map[res][1]
            if manu_thresh == None:
                thresh = threshs[res_name]
            else:
                thresh = manu_thresh
            if manu_angle_thresh == None:
                angle_thresh = angle_threshs[res_name]
            else:
                angle_thresh = manu_angle_thresh

            if thresh == None or double_unknown_mode:
                thresh = unknown_thresh
            if angle_thresh == None or double_unknown_mode:
                angle_thresh = unknown_angle_thresh
            centroid_1 = centroid_map[res]
            #centroid_1 = c_alpha_1

            gly_c = [None]
            gly_n = [None]
            c_alpha_1 = [None]
            c_beta = [None]
            for atom in atomlist:
                (atomname,x,y,z) = atomlist[atom]
                if atomname == 'CA':
                    c_alpha_1 = numpy.array([x,y,z])
                if atomname == 'CB':
                    c_beta = numpy.array([x,y,z])
                if atomname == 'N':
                    gly_n = numpy.array([x,y,z])
                if atomname == 'C':
                    gly_c = numpy.array([x,y,z])

            siss = 0.0
            error_flag = False

            for res_2 in dist_matrix[res]:
                for atom in dist_matrix[res][res_2]:
                    (d,atomname,x,y,z) = dist_matrix[res][res_2][atom]
                    #print d,res,res_2
                    #if res == res_2:
                        #print 'this should happen alot'
                    #    siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)
                    if d - vdw_radius[atomname[0]] <= thresh:
                        if angle_thresh > -1.0:
                            coord = numpy.array([x,y,z])

                            if res_name != 'GLY':
                                if c_beta[0] == None or c_alpha_1[0] == None:
                                    angle = None
                                else:
                                    angle = getCosAngle(c_beta-c_alpha_1,coord-c_alpha_1)

                                #angle = None #Calculate everything with gly-vector

                                if angle == None:
                                    error_flag = True
                                    #print c_beta,res,res_2,atomlist
                                    #This happens for resiudes (beside Glycin), where only the C-Alpha atom is given. (Note: This is the full atom case)
                                    #The solution is to handle it as Glycin
                                    
                                    if gly_c[0] == None or gly_n[0] == None or c_alpha_1[0] == None:
                                        angle = 1.0
                                    else:
                                        gly_vec = gly_vector(gly_n,gly_c,c_alpha_1)
                                        angle = getCosAngle(gly_vec,coord-centroid_1)
                                
                            else:
                                if gly_c[0] == None or gly_n[0] == None or c_alpha_1[0] == None:
                                    error_flag = True
                                    #print res,res_2
                                    angle = 1.0
                                else:    
                                    gly_vec = gly_vector(gly_n,gly_c,c_alpha_1)
                                    angle = getCosAngle(gly_vec,coord-centroid_1)
                            
                            if angle_thresh <= angle or angle == None:
                                siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)
                            elif angle_thresh == -1.0:
                                print "WADDAFAQ"
                        else:
                            #If the full sphere is taken, there is no need for calculating the angle
                            siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)
                        
            #siss = sissCorrection(siss,res_name)
            #if not error_flag: #comment this out for final version
            siss = sissCorrectionVol(siss,thresh,angle_thresh)
            siss_map[res] = siss #+ average_correction[res_name]
            #else:
            #    print "Error_flage was hit"

    return siss_map

def calcAASiss(coordinate_map,target_residues,manu_thresh=None):
    siss_map = {}    
    
    for res in target_residues:
        siss = 0.0
        res_name = coordinate_map[res][0]
        atomlist = coordinate_map[res][1]
        
        if manu_thresh == None:
            thresh = threshs[res_name]
        else:
            thresh = manu_thresh
        
        for atom in atomlist:
            atom_siss = 0.0
            (atomname,x,y,z) = atomlist[atom]    
            for res_2 in coordinate_map:
                if res != res_2:
                    atomlist_2 = coordinate_map[res_2][1]
                    for atom_2 in atomlist_2:
                        (atomname_2,x_2,y_2,z_2) = atomlist_2[atom_2]
                        diff = numpy.array([x_2,y_2,z_2]) - numpy.array([x,y,z])
                        d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                        #If one atom is further than 10A+threshhold, then ignore the whole residue
                        if d - 15.0 > thresh:
                            break
                        if d - vdw_radius[atomname_2[0]] <= thresh:                    
                            atom_siss += sphere_intersection(thresh,vdw_radius[atomname_2[0]],d)
                        
        #siss = sissCorrection(siss,res_name)
            siss += sissCorrectionVol(atom_siss,thresh,-1.0)
        siss_map[res] = siss/float(len(atomlist))       

    return siss_map

def siss(input_file=None,output_file=None,target_residues=None,chain=None,c_alpha=False,coordinate_map_path=None,unknown_mode=False,double_unknown_mode=False):
    
    cwd = os.getcwd()
    if output_file == None:
        if input_file == None:
            output_file = "%s/siss_output.tsv" % cwd
        else:
            output_file = "%s_siss.tsv" % input_file

    if input_file == None:
        if coordinate_map_path == None:
            raise NameError("No input structure given")
        else:
            coordinate_map = parseCM(coordinate_map_path)
    else:
        #TODO remove centroid calculation
        coordinate_map,protein_centroid = parsePDB(input_file,chain,c_alpha)


    if target_residues == None:
        target_residues = coordinate_map.keys()

    centroid_map = calcCentroidMap(coordinate_map,target_residues,c_alpha,double_unknown_mode=double_unknown_mode)

    dist_matrix = calcDistMatrix(coordinate_map,centroid_map,target_residues,c_alpha)

    if not unkown_mode:
        siss_map = calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha)
    else:
        siss_map = calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha,manu_thresh=unknown_thresh,manu_angle_thresh=unknown_angle_thresh,double_unknown_mode=double_unknown_mode)
    
    produceOutput(siss_map,coordinate_map,output_file)
    
if __name__ == "__main__":
    
    argv = sys.argv[1:]
    helptext = """
\tSphere InterSection Sum

This tool calculates a measure for relative solvent accessible area for single residues or whole amino acid chains.
Usage:
siss.py -i /Path/To/Input/File [-o /Path/To/Output/File] [-c chain] [-r residues] [--ca] [--cm coordinate_map]

-i:\tPath to an input file in PDB (Protein Data Bank) file format. If --cm are given, this is not needed.

-o:\tPath to the output file produced as tab separated text file.
\tDefault: *InputFile*_siss.tsv

-c:\tChain identifier of the amino acid chain, which should be analysed denoted as the chain identifiers of the ATOM records in the PDB file.
\tDefault: The first chain found in the input file

-r:\tList of residue identifiers for all residues for which the SISS value should be computed. If not given, the SISS values for all residues are computed.
\tIf given with -i the residue identifiers denote as the residues identifiers of the ATOM records in the PDB file. If given with --cm the residues should denote as the residue map keys of --cm.
\tExamples: 234,78,368 | 17 | 34,35,36,37

--ca:\tC alpha version of SISS. Needs only coordinates of the C alpha atoms, but is more inaccurate as the full atom coordinates version.

--cm:\tPath to a text file containing an atomic cooridante map as a python map: {Residue-ID:(Residue-Name,{Atom-ID:(Atom-Name,x,y,z)})}.
"""
    if len(argv) == 0:
        print helptext
        sys.exit()

    input_file = None
    output_file = None
    target_residues = None
    chain = None
    c_alpha = False
    coordinate_map_path = None

    try:
        opts,args = getopt.getopt(argv,"hr:o:i:c:",['ca','cm='])
    except getopt.GetoptError:
        print "siss.py -h"
        sys.exit(2)
    for opt,arg in opts:
        if opt == '-h':
            print(helptext)
            sys.exit()
        elif opt == "-i":
            input_file = arg
        elif opt == '-o':
            output_file = arg
        elif opt == '-r':
            target_residues = arg.split(',')
        elif opt == '-c':
            chain = arg
        elif opt == '--ca':
            c_alpha = True
        elif opt == '--cm':
            coordinate_map_path = arg

    siss(input_file,output_file,target_residues,chain,c_alpha,coordinate_map_path)        
