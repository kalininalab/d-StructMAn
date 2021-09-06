from structman.lib.output.output import OutputGenerator


def init_feature_table(feature_file):
    feature_output = OutputGenerator()
    headers = [
        'Input Protein ID', 'Primary Protein ID', 'Uniprot-Ac', 'WT Amino Acid', 'Position', 'Mut Amino Acid', 'AA change', 'Tags',
        'Distance-based classification', 'Distance-based simple classification',
        'RIN-based classification', 'RIN-based simple classification',
        'Classification confidence', 'Structure Location', 'Mainchain Location', 'Sidechain Location',
        'RSA', 'Mainchain RSA', 'Sidechain RSA', 'Amount of mapped structures',
        'Secondary structure assignment', 'IUPred value', 'Region structure type', 'Modres score',
        'Phi', 'Psi', 'KD mean',
        'Volume mean', 'Chemical distance', 'Blosum62',
        'Aliphatic change', 'Hydrophobic change', 'Aromatic change', 'Positive charged change',
        'Polar change', 'Negative charge change', 'Charged change', 'Small change', 'Tiny change', 'Total change',
        'B Factor',
        'AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality',
        'AbsoluteCentralityWithNegative', 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative',
        'AbsoluteComplexCentrality', 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality',
        'AbsoluteComplexCentralityWithNegative', 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative',
        'Intra_SSBOND_Propensity', 'Inter_SSBOND_Propensity', 'Intra_Link_Propensity', 'Inter_Link_Propensity',
        'CIS_Conformation_Propensity', 'CIS_Follower_Propensity',
        'Inter Chain Median KD', 'Inter Chain Distance Weighted KD', 'Inter Chain Median RSA', 'Inter Chain Distance Weighted RSA',
        'Intra Chain Median KD', 'Intra Chain Distance Weighted KD', 'Intra Chain Median RSA', 'Intra Chain Distance Weighted RSA',
        'Inter Chain Interactions Median', 'Inter Chain Interactions Distance Weighted',
        'Intra Chain Interactions Median', 'Intra Chain Interactions Distance Weighted'
    ]

    for chaintype in ['mc', 'sc']:
        for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
            feature_name = '%s %s score' % (chaintype, interaction_type)
            headers.append(feature_name)
            feature_name = '%s %s degree' % (chaintype, interaction_type)
            headers.append(feature_name)
            feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
            headers.append(feature_name)

    feature_output.add_headers(headers)

    feat_f = open(feature_file, 'a')
    feat_f.write(feature_output.get_header())
    return feature_output, feat_f
