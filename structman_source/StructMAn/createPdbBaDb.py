import os
import gzip
import sys

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
'N':'ASP',
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
'Q':'GLU',
'L':'LEU',
'R':'ARG',
'W':'TRP',
'N':'ASN',
'Y':'TYR',
'M':'MET',
'X':'UNK'}

def parseByAtom(path,pdb_id,great_seq_map):
    global threeToOne
    global oneToThree
    f = gzip.open(path, 'rb')
    page = f.read()
    f.close()

    page = page.decode('ascii')

    lines = page.split('\n')

    seqs = {}
    used_res = {}
    firstAltLoc = None
    for line in lines:
        record_name = line[0:6].replace(" ","")
        atom_nr = line[6:11].replace(" ","")
        atom_name = line[12:16].replace(" ","")
        res_name = line[17:20].replace(" ","")
        if record_name == 'ENDMDL':
            break
        if len(res_name) < 3:
            continue
        if len(line) > 21:
            
            chain_id = line[21]
            res_nr = line[22:27].replace(" ","")


            if record_name == "ATOM" or record_name == 'HETATM':
                altLoc = line[16]

                if firstAltLoc == None and altLoc != ' ':
                    firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                    continue
                if not chain_id in seqs:
                    seqs[chain_id] = ""
                    used_res[chain_id] = set()
                if record_name == 'HETATM' and res_name not in threeToOne:
                    continue
                if not res_nr in used_res[chain_id]:
                    if res_name in threeToOne:
                        aa = threeToOne[res_name][0]
                    else:
                        aa = 'X'
                    seqs[chain_id] += aa
                    used_res[chain_id].add(res_nr)


    for chain in seqs:        
        seq = seqs[chain]
        if len(seq) > 10:
            if not seq in great_seq_map:
                great_seq_map[seq] = [(pdb_id,chain)]
            else:
                great_seq_map[seq].append((pdb_id,chain))

def parseBySeqres(path,pdb_id,great_seq_map):
    f = gzip.open(path, 'rb')
    page = f.read()
    f.close()

    page = page.decode('ascii')

    lines = page.split('\n')

    modres = False

    seqs = {}
    #modres_map = {}
    for line in lines:
        record_name = line[0:6].replace(" ","")
        
        if record_name == "MODRES":
            print(line)

        if record_name != 'SEQRES':
            if record_name == 'ATOM':
                break
            elif modres and record_name == "HETNAM":
                #mod = line[12:15]
                #tlc = line[24:27]
                #print(mod,tlc)
                #modres_map[mod] = threeToOne[tlc]"""
                print(line)
            continue


        chain = line[11]
        if not chain in seqs:
            seqs[chain] = ""
        seq = line[19:].split()
        for aa in seq:
            if len(aa) == 3:
                if aa in threeToOne:
                    if threeToOne[aa][0] in oneToThree:
                        seqs[chain] += threeToOne[aa][0]
                    else:
                        seqs[chain] += 'X'
                else:                    
                    print(pdb_id)
                    print(aa)
                    modres = True
                    seqs[chain] += 'X'
            else:
                break

    """
    if modres:
        for chain in seqs:
            seq = seqs[chain]
            for mod in modres_map:
                seq = seq.replace("X"+mod,modres_map[mod])
            seqs[chain] = seq
    """

    for chain in seqs:
        seq = seqs[chain]
        if len(seq) > 10:
            if not seq in great_seq_map:
                great_seq_map[seq] = [(pdb_id,chain)]
            else:
                great_seq_map[seq].append((pdb_id,chain))

def parse_seq_file(seq_file,recently_modified_structures):
    f = open(seq_file,'rb')
    lines = f.read().decode('ascii').split('\n')
    f.close()

    pdb_ids = set()

    great_seq_map = {}

    seq = ""

    chains = []

    for line in lines:
        if line[0] != '>':
            seq += line
        else:
            if not seq == '':
                great_seq_map[seq] = chains
                seq = ''
                chains = []
            parts = line[1:].split(',')
            for part in parts:
                if part.count('-') < 1:
                    continue
                pdb_id,chain = part.split('-')
                if pdb_id.lower() in recently_modified_structures:
                    continue
                pdb_ids.add(pdb_id)
                if not os.path.isfile('%s/data/structures/obsolete/pdb/%s/pdb%s.ent.gz' % (pdb_path,pdb_id[1:-1].lower(),pdb_id.lower())):
                    chains.append((pdb_id,chain))

    great_seq_map[seq] = chains

    return pdb_ids,great_seq_map



pdb_path = ''

divide_folder = ''
sub_folders = []
AU_folder = ''
AU_sub_folders = []

def main(seq_file,recently_modified_structures,fromScratch=False,pdb_p = ''):
    global pdb_path
    global sub_folders
    global AU_sub_folders
    global divide_folder

    pdb_path = pdb_p

    divide_folder = "%s/data/biounit/PDB/divided/" % pdb_path

    sub_folders = os.listdir(divide_folder)

    AU_folder = "%s/data/structures/divided/pdb/" % pdb_path

    AU_sub_folders = os.listdir(AU_folder)

    if not fromScratch:
        safe_list,great_seq_map = parse_seq_file(seq_file,recently_modified_structures)
    else:
        safe_list = set()
        great_seq_map = {}

    sub_folders = sorted(sub_folders)
    bio_entries = set()

    for sub_folder in sub_folders:
        files = os.listdir("%s%s" % (divide_folder,sub_folder))

        print(sub_folder)
        for filename in files:
            if filename.count('.pdb1.gz') > 0:
                pdb_id = filename[:4].upper()
                if pdb_id in safe_list:
                    continue
                bio_entries.add(pdb_id)
                path  = "%s%s/%s" % (divide_folder,sub_folder,filename)
                parseByAtom(path,pdb_id,great_seq_map)
                #parseBySeqres(path,pdb_id,great_seq_map)

    sub_folders = sorted(AU_sub_folders)
    for sub_folder in sub_folders:
        files = os.listdir("%s%s" % (AU_folder,sub_folder))
        print(sub_folder,' AU')
        for filename in files:
            if filename.count('.ent.gz') > 0:
                pdb_id = filename[3:7].upper()
                if pdb_id in bio_entries or pdb_id in safe_list:
                    continue
                pdb_id = '%s_AU' % pdb_id
                if pdb_id in safe_list:
                    continue
                path  = "%s%s/%s" % (AU_folder,sub_folder,filename)
                parseByAtom(path,pdb_id,great_seq_map)

    entries = []

    #seed_id_map_lines = []

    id_list = []

    for seq in great_seq_map:
        if seq == '':
            continue

        for pdb,chain in  great_seq_map[seq]:
            id_list.append('%s-%s' % (pdb,chain))
            #seed_id_map_lines.append('%s\t%s' % (seed_id,','.join(id_list)))

        block_seqs = []
        while len(seq) > 0:
            m = min((len(seq),80))
            block_seqs.append(seq[:m])
            seq = seq[m:]

        if id_list == []:
            continue

        entry = ">%s\n%s" % (','.join(id_list),'\n'.join(block_seqs))
        entries.append(entry)
        id_list = []

    page = '\n'.join(entries)

    f = open(seq_file,'w')
    f.write(page)
    f.close()

    #f = open(id_map_file,'w')
    #f.write('\n'.join(seed_id_map_lines))
    #f.close()

if __name__ == "__main__":
    #seq_file = 'pdbba'
    #seq_file = '/TL/sin/work/agress/StructMAn_git/lib/base/blast_db/pdbba_mmseqs2'
    main(seq_file,set(),fromScratch=True)

