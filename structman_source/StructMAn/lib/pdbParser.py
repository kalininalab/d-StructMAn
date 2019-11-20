import urllib.request, urllib.error, urllib.parse
import os
import sys
import subprocess
import time
import gzip
import database

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

def getCurrentPDBID(pdb_id,pdb_path):
    obsolete_path = '%s/data/structures/obsolete/pdb/%s/pdb%s.ent.gz' % (pdb_path,pdb_id[1:-1].lower(),pdb_id.lower())

    if os.path.isfile(obsolete_path):
        buf = gzip.open(obsolete_path, 'rb')
        # The replacement PDB ID is always given in the second line
        next(buf)
        pdb_id = buf.readline().split()[3].decode('ascii')
        buf.close()

    return pdb_id

def getPDBBuffer(pdb_id,pdb_path,AU=False,obsolete_check=False):
    if obsolete_check:
        pdb_id = getCurrentPDBID(pdb_id, pdb_path)

    if AU:
        path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path,pdb_id[1:-1].lower(),pdb_id.lower())
        if not os.path.isfile(path):
            url = 'https://files.rcsb.org/view/%s.pdb' %pdb_id
            if pdb_path != '':
                print("Did not find asymetric unit entry in local pdb, searching online: ",url)
            try:
                request = urllib.request.Request(url)
                return urllib.request.urlopen(request)
            except:
                print("Did not find the PDB-file: %s" % pdb_id)
                return None
        else:
            return gzip.open(path, 'rb')
    else:
        path = '%s/data/biounit/PDB/divided/%s/%s.pdb1.gz' % (pdb_path,pdb_id[1:-1].lower(),pdb_id.lower())
        if not os.path.isfile(path):
            url = 'https://files.rcsb.org/view/%s.pdb1' %pdb_id
            if pdb_path != '':
                print("Did not find entry in local pdb, searching online: ",url)
            try:
                request = urllib.request.Request(url)
                return urllib.request.urlopen(request)
            except:
                print("Did not find the PDB-file: %s" % pdb_id)
                return None
        else:
            return gzip.open(path, 'rb')

def parsePDBSequence(buf):
    seq = ""
    res_pos_map = {}
    used_res = set()
    i = 1
    firstAltLoc = None
    for line in buf:
        if len(line) > 27:
            record_name = line[0:6].strip()
            if not record_name == b'ATOM' \
            and not record_name == b'HETATM':
                continue

            altLoc = chr(line[16])
            if firstAltLoc == None and altLoc != ' ':
                firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue

            res_name = line[17:20].strip().decode('ascii')
            if record_name == b'HETATM':
                if res_name not in threeToOne:
                    continue

            chain_id = line[21].decode('ascii')
            if chain_id != chain:
                continue

            atom_nr = line[6:11].strip().decode('ascii')
            res_nr = line[22:27].strip().decode('ascii') # this includes the insertion code

            if res_nr not in used_res:
                aa = threeToOne[res_name][0]
                seq = seq + aa
                used_res.add(res_nr)
                res_pos_map[res_nr] = i
                i += 1

    return seq,res_pos_map

def getSequenceAU(pdb_id,chain,pdb_path):
    buf = getPDBBuffer(pdb_id,pdb_path,AU=True,obsolete_check=False)
    if buf == None:
        return "", {}
    seq,res_pos_map = parsePDBSequence(buf,chain)
    buf.close()
    return seq,res_pos_map

def getSequencePlain(pdb_id,chain,pdb_path):
    buf = getPDBBuffer(pdb_id,pdb_path,AU=False,obsolete_check=True)
    if buf == None:
        return "", {}
    seq,res_pos_map = parsePDBSequence(buf,chain)
    buf.close()
    return seq,res_pos_map

def getSequence(pdb_id,chain,pdb_path,AU=False):
    if AU:
        seq,res_pos_map = getSequenceAU(pdb_id,chain,pdb_path)
    else:
        seq,res_pos_map = getSequencePlain(pdb_id,chain,pdb_path)

    #Catch empty sequences and try to find the sequence in the asymetric unit
    if not AU and seq == "":
        return getSequenceAU(pdb_id,chain,pdb_path)

    return seq,res_pos_map

#called by serializedPipeline
def getSequences(pdb_dict,pdb_path,AU=False):
    pdb_sequence_map = {}
    pdb_pos_map = {}
    for pdb_chain_tuple in pdb_dict:
        [pdb,chain] = pdb_chain_tuple.split(':')
        seq,res_pos_map = getSequence(pdb,chain,pdb_path,AU=AU)
        pdb_sequence_map[pdb_chain_tuple] = seq,None,None
        pdb_pos_map[pdb_chain_tuple] = res_pos_map
    #print pdb_sequence_map#,pdb_pos_map
    return pdb_sequence_map,pdb_pos_map

#called by serializedPipeline
def getSequencesPlain(pdb_dict,pdb_path):
    pdb_sequence_map = {}
    for pdb_chain_tuple in pdb_dict:
        [pdb,chain] = pdb_chain_tuple.split(':')
        seq,res_pos_map = getSequencePlain(pdb,chain,pdb_path)
        pdb_sequence_map[pdb_chain_tuple] = seq
    #print "Plain",pdb_sequence_map
    return pdb_sequence_map

def parsePDBSequence(buf,chain):
    seq = ""
    res_pos_map = {}
    used_res = set()
    i = 1
    firstAltLoc = None
    for line in buf:
        line = line.decode('ascii')
        #print(line)
        if len(line) > 26:
            record_name = line[0:6].strip()
            if not record_name == 'ATOM' \
            and not record_name == 'HETATM':
                continue

            altLoc = line[16]
            if firstAltLoc == None and altLoc != ' ':
                firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue

            res_name = line[17:20].strip()
            if record_name == 'HETATM':
                if res_name not in threeToOne:
                    continue

            chain_id = line[21]
            if chain_id != chain:
                continue

            atom_nr = line[6:11].strip()
            res_nr = line[22:27].strip() # this includes the insertion code

            if res_nr not in used_res:
                aa = threeToOne[res_name][0]
                seq = seq + aa
                used_res.add(res_nr)
                res_pos_map[res_nr] = i
                i += 1

    return seq,res_pos_map


#called by serializedPipeline
#called by templateFiltering
def standardParsePDB(pdb_id,pdb_path,obsolete_check=False):
    AU = False
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getPDBBuffer(pdb_id,pdb_path,AU=AU,obsolete_check=obsolete_check)

    chain_ids = set()
    newlines = []

    current_serial = 1
    firstAltLoc = None
    for line in buf:
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
            line = line.rstrip(b'\n')
        else:
            continue

        if record_name == b'ENDMDL':
            newlines.append(line)
            break
        elif record_name == b'TER':
            # This can happen, if the whole chain consists of boring ligands,
            # which is very boring
            if chain_id not in chain_ids:
                continue
        elif record_name == b'MODRES':
            pass
        elif record_name == b'ATOM' or record_name == b'HETATM':
            altLoc = chr(line[16])
            if firstAltLoc == None and altLoc != ' ':
                firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue
            res_name = line[17:20].strip().decode('ascii')
            chain_id = line[21:22].decode('ascii')
            res_nr = line[22:27].strip().decode('ascii')

            if record_name == b'HETATM':
                if res_name not in threeToOne and boring(res_name):
                    continue

            if not chain_id in chain_ids:
                chain_ids.add(chain_id)

        else:
            # We are not interested in other record types
            continue

        # Adjust serials after skipping records
        if current_serial > 99999:
            newline = record_name[0:5].ljust(5, b' ') + str(current_serial).encode('ascii').rjust(6, b' ')
        else:
            newline = record_name.ljust(6, b' ') + str(current_serial).encode('ascii').rjust(5, b' ')
        newline += line[11:]
        newlines.append(newline)
        current_serial += 1

    template_page = b'\n'.join(newlines).decode('ascii')

    return template_page

#called by serializedPipeline
def getStandardizedPdbFile(pdb_id,chain,structure,pdb_path):
    times = [0.0,0.0,0.0]
    t0 = time.time()
    AU = False

    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getPDBBuffer(pdb_id,pdb_path,AU=AU,obsolete_check=False)

    if buf == None:
        return "",'',"Did not find the PDB-file: %s" % (pdb_id),times

    structure['IAP'] = []
    t1 = time.time()

    chain_type_map = {}
    modres_map = {}

    lig_set = {}

    newlines = []

    current_serial = 1
    firstAltLoc = None
    for line in buf:
        line = line.decode('ascii')
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
            line = line.rstrip('\n')
        else:
            continue

        if record_name == 'ENDMDL':
            newlines.append(line)
            break
        elif record_name == 'TER':
            # This can happen, if the whole chain consists of boring ligands,
            # which is very boring
            if chain_id not in chain_type_map:
                continue
        elif record_name == 'MODRES':
            chain_id = line[16]
            res_nr = line[18:23].strip()
            res_name = line[24:27].strip()
            if (chain_id,res_nr) not in modres_map:
                modres_map[(chain_id,res_nr)] = res_name
        elif record_name == 'ATOM' or record_name == 'HETATM':
            altLoc = line[16]
            if firstAltLoc == None and altLoc != ' ':
                firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue
            res_name = line[17:20].strip()
            chain_id = line[21:22]
            res_nr = line[22:27].strip()

            #occupancy = float(line[54:60].replace(" ",""))

            if chain_id not in chain_type_map:
                if record_name == 'ATOM':
                    if len(res_name) == 1:
                        chain_type = 'RNA'
                    elif len(res_name) == 2:
                        chain_type = 'DNA'
                    elif len(res_name) == 3:
                        chain_type = 'Protein'
                    chain_type_map[chain_id] = chain_type
                elif record_name == 'HETATM':
                    # For a hetero peptide 'boring' hetero amino acids are
                    # allowed as well as other non boring molecules not in
                    # threeToOne, which are hopefully some kind of anormal amino
                    # acids
                    if res_name in threeToOne or not boring(res_name):
                        chain_type = 'Peptide'
                        chain_type_map[chain_id] = chain_type
            elif record_name == 'ATOM' and chain_type_map[chain_id] == 'Peptide':
                if len(res_name) == 1:
                    chain_type = "RNA"
                elif len(res_name) == 2:
                    chain_type = "DNA"
                elif len(res_name) == 3:
                    chain_type = "Protein"
                chain_type_map[chain_id] = chain_type

            if record_name == 'HETATM':
                # modified residue are not parsed as ligands
                if res_name not in threeToOne and not boring(res_name):
                    if chain_id not in lig_set:
                        lig_set[chain_id] = set()
                    if not res_nr in lig_set[chain_id]:
                        lig_set[chain_id].add(res_nr)

                        if (chain_id,res_nr) in modres_map:
                            continue

                        structure['IAP'].append(["Ligand",res_name,res_nr,chain_id])

        else:
            # We are not interested in other record types
            continue

        # Adjust serials after skipping records
        if current_serial > 99999:
            newline = record_name[0:5].ljust(5, ' ') + str(current_serial).rjust(6, ' ')
        else:
            newline = record_name.ljust(6, ' ') + str(current_serial).rjust(5, ' ')
        newline += line[11:]
        newlines.append(newline)
        current_serial += 1

    pep_not_lig = []
    for chain_id in chain_type_map:
        chain_type = chain_type_map[chain_id]
        structure['IAP'].append([chain_type,chain_id])

    for chain_id in chain_type_map:
        chain_type = chain_type_map[chain_id]
        if chain_type == 'Peptide':
            for i,iap in enumerate(structure['IAP']):
                if iap[0] == 'Ligand':
                    if iap[3] == chain_id:
                        pep_not_lig.append(i)

    pep_not_lig = sorted(pep_not_lig,reverse=True)

    for i in pep_not_lig:
        del structure['IAP'][i]

    # Happens for repeated chains in asymetric units, which do not occur in the
    # biological assembly
    irregular_homo_chains = []
    for i,homo_chain in enumerate(structure['Oligo']):
        if not homo_chain in chain_type_map:
            irregular_homo_chains.append(i)

    irregular_homo_chains = sorted(irregular_homo_chains,reverse=True)
    for i in irregular_homo_chains:
        del structure['Oligo'][i]

    t2 = time.time()

    template_page = '\n'.join(newlines)

    t3 = time.time()
    times = [t1-t0,t2-t1,t3-t2]

    return template_page,structure,None,times

#called by database
def parseLigandDB(smiles_path,inchi_path):
    ligand_db = {}

    f = open(smiles_path,'rt')
    for line in f:
        words = line.rstrip('\n').split('\t')
        if len(words) < 2:
            continue
        smiles = words[0]
        name = words[1]
        ligand_db[name] = smiles
    f.close()

    f = open(inchi_path,'rt')
    for line in f:
        words = line.rstrip('\n').split('\t')
        if len(words) < 2:
            continue
        inchi = words[0]
        name = words[1]
        smiles = ligand_db[name]
        ligand_db[name] = (smiles,inchi)
    f.close()

    return ligand_db

#called by database
def updateLigandDB(new_ligands,smiles_path,inchi_path):
    smiles_lines = []
    inchi_lines = []
    for (name,smiles,inchi) in new_ligands:
        smiles_lines.append('%s\t%s' % (smiles,name))
        inchi_lines.append('%s\t%s' % (inchi,name))

    f = open(smiles_path,'a')
    f.write('\n'.join(smiles_lines)+'\n')
    f.close()

    f = open(inchi_path,'a')
    f.write('\n'.join(inchi_lines)+'\n')
    f.close()

    return

#get called by database
def getSI(pdb_id,name,res,chain,pdb_path):
    cwd = os.getcwd()

    AU = False
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getPDBBuffer(pdb_id,pdb_path,AU=AU,obsolete_check=False)

    newlines = []

    for line in buf:
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
        else:
            continue
        if record_name == b'ENDMDL':
            break

        if record_name == b'HETATM':
            res_name = line[17:20].strip().decode('ascii')
            chain_id = line[21:22].decode('ascii')
            res_nr = line[22:27].strip().decode('ascii')
            if chain == chain_id and name == res_name and res_nr == res:
                newlines.append(line)

    buf.close()

    page = (b''.join(newlines)).decode('ascii')

    FNULL = open(os.devnull, 'w')
    try:
        inchiproc = subprocess.Popen(["babel","-i","pdb","-o","inchi"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        inchi,i_err = inchiproc.communicate(page)
        smilesproc = subprocess.Popen(["babel","-i","pdb","-o","smi"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        smiles,s_err = smilesproc.communicate(page)
    except:
        print("[WARNING] Babel Package seems not to be installed")
        return ("","")

    #print(i_err)
    #print(s_err)
    smiles = smiles.split()
    if len(smiles) > 0:
        smiles = smiles[0]
    else:
        print("Ligand not found in database and Babel failed: ",name,res,chain,pdb_id,i_err,s_err)

    inchi = inchi.replace('\n','')

    return (smiles,inchi)

def getPDBHeaderBuffer(pdb_id,pdb_path):
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]

    path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path,pdb_id[1:-1].lower(),pdb_id.lower())
    if not os.path.isfile(path):
        url = 'https://files.rcsb.org/header/%s.pdb' % pdb_id
        if pdb_path != '':
            print("Did not find asymetric unit entry in local pdb, searching online: ",url)
        try:
            request = urllib.request.Request(url)
            return urllib.request.urlopen(request)
        except:
            print("Did not find the PDB-header: %s" % pdb_id)
            return None
    else:
        return gzip.open(path, 'rb')

def boring(abr):
    if abr in database.metal_atoms or abr in database.ion_atoms:
        return False
    if abr in boring_ligands:
        return True
    if len(abr) < 3:
        return True
    return False

#called by templateSelection
#called by templateFiltering
def getInfo(pdb_id,pdb_path):
    buf = getPDBHeaderBuffer(pdb_id,pdb_path)

    if buf == None:
        return None,{}
    #print(pdb_id)
    abort = False
    resolution = None

    not_crystal = False

    homomer_dict = {}

    multiple_chain_line = False

    for line in buf:
        if abort:
            break

        line = line.rstrip(b'\n')
        words = line.split()

        if len(words) > 0:
            if words[0] == b'COMPND':
                if line[11:16] == b'CHAIN' or multiple_chain_line:
                    if line.count(b';') == 0:
                        multiple_chain_line = True
                    else:
                        multiple_chain_line = False

                    chains = line[10:].replace(b'CHAIN:',b'').replace(b' ',b'').replace(b';',b'').decode('ascii').split(',')

                    for chain in chains:
                        homomer_dict[chain] = chains

            #Filter out NMR-structures, if you decide to use them in future, remember Recoord
            elif words[0] == b'EXPDTA':
                if words[2] == b'NMR':
                    resolution = 2.5
                    not_crystal = True

                if words[1] == b'POWDER':
                    resolution = 100.0
                    abort = True

                if words[1] == b'SOLUTION' and (words[2] == b'SCATTERING' or words[2] == b'SCATTERING;' or words[2] == b'NMR;' or words[2] == b'NMR'):
                    resolution = 2.5
                    not_crystal = True

                #if line.count(b'ELECTRON MICROSCOPY'):
                #    resolution = 2.5
                #    not_crystal = True

                if line.count(b'NMR') > 0:
                    resolution = 2.5
                    not_crystal = True


            elif words[0] == b'REMARK':
                if not_crystal:
                    continue
                if len(words) > 2:
                    if words[1] == b'2':
                        if words[2] == b'RESOLUTION.':
                            if words[3] == b'NULL' or words[3] == b'NOT':
                                resolution = None
                            else:
                                resolution = float(words[3])
                            abort = True

            #Filter out C-alpha only structures
            elif words[0] == b'MDLTYP':
                if len(words) > 2:
                    if b' '.join(words[1:3]) == b'CA ATOMS':
                        resolution = 100.0
                        abort = True

            elif words[0] == b'SEQRES':
                break


            #the authors of 1L9U were not able to make a proper MDLTYP entry
            if words[0] == b'REMARK':
                if b'COORDINATES' in words \
                and b'CONTAIN' in words \
                and b'ONLY' in words \
                and b'CA' in words:
                    resolution = 100.0
                    abort = True

    if resolution == None:
        resolution = 100.0
    try:
        resolution = float(resolution)
    except:
        raise NameError("Resolution buggy for: %s" % pdb_id)
    return resolution,homomer_dict

def colorStructure(pdb_id,chain,outfile,db,cursor,pdb_path):
    colors = {'Surface' : b' 10.00',
              'Core' : b' 90.00',
              'Contact Ligand' : b' 70.00',
              'Contact Protein' : b' 50.00',
              'Contact DNA' : b' 30.00',
              'Contact RNA' : b' 40.00'}

    res_class_map = database.getStructureClassification(pdb_id,chain,db,cursor)

    AU = False
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getPDBBuffer(pdb_id,pdb_path,AU=AU,obsolete_check=False)

    if buf == None:
        return ""#,'',"Did not find the PDB-file: %s" % pdb_id

    newlines = []

    for line in buf:
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
            line = line.rstrip(b'\n')

            if record_name == b'ENDMDL':
                newlines.append(line)
                break
            elif record_name == b'ATOM' or record_name == b'MODRES' or record_name == b'HETATM':
                chain_id = line[21:22].decode('ascii')
                res_nr = line[22:27].strip().decode('ascii')
                if chain_id == chain and res_nr in res_class_map:
                    res_class = res_class_map[res_nr]
                    color_code = colors[res_class]
                    line = line[:60] + color_code + line[66:]

        newlines.append(line)

    f = open(outfile,'wb')
    f.write(b'\n'.join(newlines))
    f.close()

chain_id_list = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]

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
