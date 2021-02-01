#!/usr/bin/python3
import argparse
import gzip
import sys, os
import igraph
from collections import defaultdict

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

boring_ligands = set(["HOH",
		            "MPT",
		            "IPA",
		            "SO4",
		            "GOL",
		            "PO4",
		            "WAT",
		            "SOL",
		            "ACY",
		            "FMT",
		            "EDO",
		            "EPE",
		            "DMS",
		            "ACT",
		            "BME",
		            "DOD",
		            "CME",
		            "NTB",
		            "ARF",
		            "PGR",
		            "MPD",
		            "NAG",
		            "GAL",
		            "TAM",
		            "D1D",
		            "IUM",
		            "CO3",
		            "MES",
		            "NCO",
		            "MES",
		            "BMA",
		            "MAN",
		            "FUC",
		            "LBT",
		            "DTV",
		            "DTU",
		            "ACE",
		            "ANL",
		            "BMT",
		            "PG6",
		            "LNK",
		            "C8E",
		            "PGW",
		            "PG0",
		            "2PE",
		            "UPL",
		            "LMG",
		            "NDG",
		            "PG6",
		            "D12",
		            "DTT",
		            "LMT",
		            "I42",
		            "TRD",
		            "NAG",
		            "2HP",
		            "12P",
		            "CXE",
		            "2PE",
		            "M2M",
		            "BNG",
		            "BNG",
		            "TRD",
		            "BNG",
		            "BCR",
		            "PGW",
		            "12P",
		            "PE4",
		            "SPM",
		            "MES",
		            "CLA",
		            "TAR",
		            "PG6",
		            "12P",
		            "SBY",
		            "PE4",
		            "BCR",
		            "POP",
		            "PLX",
		            "SBY",
		            "BEZ",
		            "MYR",
		            "P6G",
		            "PE5",
		            "LMT",
		            "SIA",
		            "1PE",
		            "LDA",
		            "PGV",
		            "C5P",
		            "1PE",
		            "LMT",
		            "AE3",
		            "15P",
		            "UNK",
		            "PE4",
		            "P4C",
		            "TAM",
		            "CHL",
		            "SBY",
		            "PEE",
		            "D12",
		            "PGW",
		            "HEX",
		            "CXE",
		            "P6G",
		            "DAO",
		            "6PL",
		            "7PE",
		            "LDA",
		            "1PE",
		            "ETE",
		            "PE4",
		            "GDL",
		            "PTL",
		            "MPO",
		            "DTT",
		            "OCT",
		            "BOG",
		            "BOG",
		            "PE3",
		            "PG4",
		            "BCR",
		            "PEE",
		            "PLM",
		            "P33",
		            "PG6",
		            "UNK",
		            "HP6",
		            "DPF",
		            "UNL",
		            "D12",
		            "5AX",
		            "PGE",
		            "WO4",
		            "15P",
		            "HTG",
		            "12P",
		            "CLA",
		            "NAG",
		            "UND",
		            "CHD",
		            "BCR",
		            "MYR",
		            "DAO",
		            "LMG",
		            "2PE",
		            "MOE",
		            "C8E",
		            "BCR",
		            "HF3",
		            "LMT",
		            "2PE",
		            "12P",
		            "UNK",
		            "1PE",
		            "CE9",
		            "HEZ",
		            "P33",
		            "PEE",
		            "LDA",
		            "GOL",
		            "PGW",
		            "PG4",
		            "CXE",
		            "PE3",
		            "UNK",
		            "LI1",
		            "CAC",
		            "PO4",
		            "CXE",
		            "C10",
		            "D12",
		            "DIO",
		            "BOG",
		            "PGE",
		            "UPL",
		            "1PE",
		            "P33",
		            "PGV",
		            "CL1",
		            "P6G",
		            "C8E",
		            "ETE",
		            "P6G",
		            "P33",
		            "C8E",
		            "ACD",
		            "LMG",
		            "D1D",
		            "PE3",
		            "BNG",
		            "12P",
		            "KMB",
		            "PG4",
		            "EGC",
		            "P4C",
		            "MYS",
		            "UPL",
		            "TBU",
		            "CXE",
		            "LMG",
		            "12P",
		            "GOL",
		            "F3S",
		            "UPL",
		            "PEE",
		            "PG4",
		            "15P",
		            "5AX",
		            "NAG",
		            "BCR",
		            "LMG",
		            "15P",
		            "D9G",
		            "PG0",
		            "PGM",
		            "2HA",
		            "NAG",
		            "BCR",
		            "DTU",
		            "NAG",
		            "NDG",
		            "PE3",
		            "NBU",
		            "LHG",
		            "D12",
		            "PG0",
		            "CXE",
		            "BCR",
		            "MRD",
		            "PT5",
		            "DKA",
		            "P4C",
		            "SO4",
		            "OCT",
		            "LI1",
		            "15P",
		            "DD9",
		            "FTT",
		            "PGW",
		            "OC9",
		            "1PE",
		            "UNK",
		            "PG4",
		            "5GP",
		            "WO6",
		            "LI1",
		            "PG6",
		            "TRD",
		            "UNL",
		            "NAG",
		            "1PG",
		            "15P",
		            "1PE",
		            "SGN",
		            "UPL",
		            "GOL",
		            "PE5",
		            "P3G",
		            "UNK",
		            "PEU",
		            "CDL",
		            "CDL",
		            "UNL",
		            "BTB",
		            "LHG",
		            "GOL",
		            "NDG",
		            "U5P",
		            "SOH",
		            "BCR",
		            "CDL",
		            "PEU",
		            "PIO",
		            "BCR",
		            "2PE",
		            "PE4",
		            "C10",
		            "P33",
		            "PAM",
		            "15P",
		            "ZRC",
		            "P6G",
		            "MYR",
		            "MRD",
		            "D10",
		            "PE3",
		            "DTD",
		            "A2G",
		            "UNL",
		            "LMG",
		            "CXS",
		            "LHG",
		            "PEF",
		            "2PE",
		            "TOE",
		            "PEF",
		            "DEP",
		            "D12",
		            "PLM",
		            "PSC",
		            "PE3",
		            "LMG",
		            "BOG",
		            "PE4",
		            "LMG",
		            "PGE",
		            "PE8",
		            "PEG",
		            "2PE",
		            "15P",
		            "ETE",
		            "HTG",
		            "BCR",
		            "D9G",
		            "PHF",
		            "ETE",
		            "PG4",
		            "NDG",
		            "UPL",
		            "MYR",
		            "3GR",
		            "BCR",
		            "PGR",
		            "PG6",
		            "UNL",
		            "12P",
		            "BNG",
		            "YBT",
		            "7PE",
		            "TAR",
		            "PG6",
		            "VO4",
		            "B3P",
		            "PGR",
		            "TLA",
		            "NDG",
		            "PO3",
		            "LI1",
		            "PG6",
		            "ACD",
		            "HTO",
		            "PLM",
		            "PLM",
		            "TRS",
		            "TRD",
		            "TAR",
		            "P6G",
		            "PGE",
		            "TTP",
		            "D10",
		            "1PE",
		            "PPI",
		            "UNL",
		            "ETE",
		            "PLM",
		            "2PE",
		            "C10",
		            "SQD",
		            "1PG",
		            "C8E",
		            "ETX",
		            "2PE",
		            "PE3",
		            "PGW",
		            "STE",
		            "SRT",
		            "UNL",
		            "BCR",
		            "LI1",
		            "PGO",
		            "PE4",
		            "DTV",
		            "BCR",
		            "12P",
		            "PEG",
		            "CXE",
		            "D12",
		            "PSE",
		            "BU1",
		            "PYR",
		            "TRD",
		            "PG6",
		            "SF4",
		            "P6G",
		            "2PE",
		            "P6G",
		            "PEE",
		            "CXE",
		            "AE3",
		            "PLM",
		            "UNK",
		            "7PE",
		            "BCR",
		            "1PE",
		            "BCN",
		            "PG4",
		            "MYR",
		            "EPE",
		            "PEG",
		            "MGE",
		            "1PG",
		            "CBS",
		            "LDA",
		            "MPD",
		            "P6G",
		            "C8E",
		            "P33",
		            "1PG",
		            "PE5",
		            "PLM",
		            "CL1",
		            "UPL",
		            "UNL",
		            "BCR",
		            "PGO",
		            "GOL",
		            "PE4",
		            "SPD",
		            "MPD",
		            "BCR",
		            "MAG",
		            "LDA",
		            "ETE",
		            "BCR",
		            "PG4",
		            "MRD",
		            "2PE",
		            "WO3",
		            "1PE",
		            "PGW",
		            "6JZ",
		            "LMU",
		            "UNL",
		            "PGE",
		            "BMA",
		            "P6G",
		            "I3P",
		            "15P",
		            "1PE",
		            "SF4",
		            "OCT",
		            "15P",
		            "BCR",
		            "1PG",
		            "PE3",
		            "1PG",
		            "D10",
		            "NGA",
		            "12P",
		            "LI1",
		            "PE5",
		            "UNL",
		            "LHG",
		            "EGC",
		            "NAG",
		            "PEF",
		            "1PG",
		            "LAP",
		            "LHG",
		            "BOG",
		            "PE5",
		            "HEM",
		            "UNL",
		            "2PO",
		            "TRD",
		            "SO4",
		            "LMG",
		            "GOL",
		            "PEE",
		            "GOL",
		            "PG5",
		            "BCR",
		            "MC3",
		            "BCR",
		            "PG4",
		            "MG8",
		            "1PG",
		            "SO3",
		            "MXE",
		            "DTT",
		            "12P",
		            "LMT",
		            "2PE",
		            "NGA",
		            "GLC",
		            "NAG",
		            "UPL",
		            "P33",
		            "UNK",
		            "LHG",
		            "MAN",
		            "LMG",
		            "HTG",
		            "1PG",
		            "LMT",
		            "MPD",
		            "C10",
		            "C8E",
		            "1PE",
		            "UPL",
		            "P33",
		            "PG4",
		            "BCR",
		            "12P",
		            "DTT",
		            "PEU",
		            "PI ",
		            "D10",
		            "PHS",
		            'D9G',
		            "P33",
		            "P4C",
		            "OCT",
		            "PLM",
		            "LDA",
		            "DKA",
		            "15P",
		            "1PE",
                    'CSD',
		            'HYP',
		            'BMT',
		            '5HP',
		            'ABA',
		            'AIB',
		            'CSW',
		            'OMT',
		            'OCS',
		            'DAL',
		            'DAR',
		            'DSG',
		            'DSP',
		            'DCY',
		            'DGL',
		            'DGN',
		            'DHI',
		            'DIL',
		            'DIV',
		            'DLE',
		            'DLY',
		            'DPN',
		            'DPR',
		            'DSN',
		            'DTH',
		            'DTR',
		            'DTY',
		            'DVA',
		            'CGU',
		            'KCX',
		            'LLP',
		            'CXM',
		            'FME',
		            'MLE',
		            'MVA',
		            'NLE',
		            'PTR',
		            'ORN',
		            'SEP',
		            'TPO',
		            'PCA',
		            'SAR',
		            'CEA',
		            'CSO',
		            'CZ2',
		            'CSS',
		            'CSX',
		            'CME',
		            'TYS',
		            'TPQ',
		            'STY',
		            'MSE',
		            'SAC',
		            '5UA',])

class ChainGraph():
    def __init__(self):
        self.res_to_vertex = defaultdict(lambda: len(self.res_to_vertex))
        self.edges = []
        self.weights = []
        self.weights2 = []

def normalize_centrality(btw_cent, min_btw, max_btw, num_residues):
    # Betweenness centrality only makes sense if there is at least one vertex connected to two edges, otherwise everything will be 0
    if max_btw == 0:
        return (None, None, None)

    btw_cent_norm_len = btw_cent / ((num_residues - 1) * (num_residues - 2) / 2.0)

    if min_btw == max_btw:
        btw_cent_norm_max = 1.0
    else:
        btw_cent_norm_max = (btw_cent - min_btw) / (max_btw - min_btw)

    return (btw_cent, btw_cent_norm_len, btw_cent_norm_max)


def main(path):
    dir_name = os.path.dirname(path)
    pdb_id = os.path.basename(path)[0:4]

    scores_path = dir_name + '/' + pdb_id + '_intsc.ea'
    if os.path.exists(scores_path):
        scores_file = open(scores_path, 'r')
    else:
        scores_path = scores_path + ".gz"
        scores_file = gzip.open(scores_path, 'rt')

    out_path = dir_name + '/' + pdb_id + '_btw_cent.txt.gz'
    out_file = gzip.open(out_path, 'wt')


    chains = defaultdict(ChainGraph)
    all_chains = ChainGraph()

    #Skip header
    next(scores_file)

    for line in scores_file:
        res1, int_type, res2, score = line.strip().split()

        if int_type != "(combi:all_all)":
            continue

        # Don't include clashing residues in our network
        score = float(score)
        score2 = score
        if score < 0:
            score = 0

        score2 = abs(score2)

        chain1 = res1.split(":")[0]
        chain2 = res2.split(":")[0]

        res_type1 = res1.split(":")[3]
        res_type2 = res2.split(":")[3]
        if res_type1 in boring_ligands and res_type1 not in threeToOne:
            continue
        if res_type2 in boring_ligands and res_type2 not in threeToOne:
            continue

        vertex_id1 = chains[chain1].res_to_vertex[res1]
        vertex_id2 = chains[chain2].res_to_vertex[res2]

        # Ignore contacts between chains unless we are calculating complex centralities
        if chain1 == chain2:
            chains[chain1].edges.append((vertex_id1, vertex_id2))
            chains[chain1].weights.append(float('inf') if score == 0 else 1.0 / score)
            chains[chain1].weights2.append(float('inf') if score2 == 0 else 1.0 / score2)

        complex_vertex_id1 = all_chains.res_to_vertex[res1]
        complex_vertex_id2 = all_chains.res_to_vertex[res2]

        all_chains.edges.append((complex_vertex_id1, complex_vertex_id2))
        all_chains.weights.append(float('inf') if score == 0 else 1.0 / score)
        all_chains.weights2.append(float('inf') if score2 == 0 else 1.0 / score2)

    scores_file.close()

    out_file.write('#Residue\t'
                   'AbsoluteCentrality\tLengthNormalizedCentrality\tMinMaxNormalizedCentrality\t'
                   'AbsoluteCentralityWithNegative\tLengthNormalizedCentralityWithNegative\tMinMaxNormalizedCentralityWithNegative\t'
                   'AbsoluteComplexCentrality\tLengthNormalizedComplexCentrality\tMinMaxNormalizedComplexCentrality\t'
                   'AbsoluteComplexCentralityWithNegative\tLengthNormalizedComplexCentralityWithNegative\tMinMaxNormalizedComplexCentralityWithNegative\n')

    num_complex_residues = len(all_chains.res_to_vertex)
    complex_graph = igraph.Graph(n=len(all_chains.res_to_vertex), edges=all_chains.edges, directed=False)
    if num_complex_residues < 2 or len(all_chains.edges) < 1:
        return

    try:
        complex_betweenness = complex_graph.betweenness(directed=False, weights=all_chains.weights)
        complex_betweenness_with_negative = complex_graph.betweenness(directed=False, weights=all_chains.weights2)
    except:
        print('Centrality error:',pdb_id,num_complex_residues,len(all_chains.edges))

    min_complex_btw_cent = min(complex_betweenness)
    max_complex_btw_cent = max(complex_betweenness)
    min_complex_btw_cent_with_negative = min(complex_betweenness_with_negative)
    max_complex_btw_cent_with_negative = max(complex_betweenness_with_negative)

    for chain, chain_graph in sorted(chains.items()):
        residues = [ r for r, v in sorted(chain_graph.res_to_vertex.items(), key=lambda x: x[1]) ]
        num_residues = len(residues)
        if num_residues < 2 or len(chain_graph.edges) < 1:
            continue
        graph = igraph.Graph(n=num_residues, edges=chain_graph.edges, directed=False)

        try:
            betweenness = graph.betweenness(directed=False, weights=chain_graph.weights)
            betweenness_with_negative = graph.betweenness(directed=False, weights=chain_graph.weights2)
        except:
            print('Centrality error:',pdb_id,chain,num_residues,len(chain_graph.edges))

        min_btw_cent = min(betweenness)
        max_btw_cent = max(betweenness)
        min_btw_cent_with_negative = min(betweenness_with_negative)
        max_btw_cent_with_negative = max(betweenness_with_negative)

        res_names = []
        res_centralities = []

        for res, btw_cent, btw_cent_with_negative in sorted(zip(residues, betweenness, betweenness_with_negative), key=lambda x: x[1], reverse=True):
            centralities = []
            centralities.extend(normalize_centrality(btw_cent, min_btw_cent, max_btw_cent, num_residues))
            centralities.extend(normalize_centrality(btw_cent_with_negative, min_btw_cent_with_negative, max_btw_cent_with_negative, num_residues))

            complex_vertex_id = all_chains.res_to_vertex[res]
            centralities.extend(normalize_centrality(complex_betweenness[complex_vertex_id],
                                                     min_btw_cent,
                                                     max_btw_cent,
                                                     num_complex_residues))
            centralities.extend(normalize_centrality(complex_betweenness_with_negative[complex_vertex_id],
                                                     min_btw_cent_with_negative,
                                                     max_btw_cent_with_negative,
                                                     num_complex_residues))
            if any(centralities) != None:
                res_names.append(res)
                res_centralities.append(centralities)

        for entries in zip(res_names, res_centralities):
            out_file.write(entries[0] + '\t' + '\t'.join(map(str, entries[1])) + '\n')

    out_file.close()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate betweenness centrality for a given RIN')
    parser.add_argument('in_file', metavar='input.sif', help='input .sif or .sif.gz file')
    args = parser.parse_args()


    main(args.in_file)
