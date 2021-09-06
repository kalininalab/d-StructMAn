NON_BORING_SHORT_LIGANDS = set(['PC', 'NO', 'VX', 'VR'])

BORING_LIGANDS = set([
    "HOH",
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
    '5UA',
])
