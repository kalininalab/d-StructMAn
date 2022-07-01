import math
import os
import psutil
from zlib import adler32

from structman import settings

import structman.lib.sdsc as sdsc

import zstd
import msgpack
import pickle
import pickletools
import sys

structure_id_types = [
    'PDB',
    'Asymmetric unit PDB',
    'Alphafold model'
]

def identify_structure_id_type(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    if structe_id_type_key is None:
        return 'Unknown structure identifier'
    return structure_id_types[structe_id_type_key]

def identify_structure_id_type_key(structure_id):
    l = len(structure_id)
    if l == 7:
        if structure_id[-3:] == '_AU':
            return 1

    if len(structure_id) == 4:
        return 0

    if structure_id[:3] == 'AF-':
        return 2

    return None

def is_alphafold_model(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    return (structe_id_type_key == 2)

def alphafold_model_id_to_file_path(model_id, config):
    uniprot_ac = model_id.split('-')[1]
    topfolder_id = uniprot_ac[-2:]
    subfolder_id = uniprot_ac[-4:]

    folder_path = f'{config.path_to_model_db}/{topfolder_id}/{subfolder_id}'

    for fn in os.listdir(folder_path):
        if fn.count(model_id) == 1 and fn.count('.pdb.gz') == 1:
            return f'{config.path_to_model_db}/{topfolder_id}/{subfolder_id}/{fn}'

    config.errorlog.add_warning(f'Did not find file {model_id} in alphafold model db: {config.path_to_model_db}')
    return None

def is_alphafold_db_valid(config):
    if not os.path.exists(config.path_to_model_db):
        return False
    if not os.path.exists(f'{config.path_to_model_db}/00'):
        return False
    return True


def resolve_path(path):
    return os.path.abspath(os.path.realpath(path))

def median(l):
    n = len(l)
    l = sorted(l)
    if n == 1:
        return l[0]
    if n == 0:
        return None
    if n % 2 == 0:
        med = (l[(n // 2) - 1] + l[n // 2]) / 2.0
    else:
        med = l[(n - 1) // 2]
    return med


def distance(coord1, coord2):
    diff = [coord1[0] - coord2[0], coord1[1] - coord2[1], coord1[2] - coord2[2]]
    return math.sqrt(diff[0]**2.0 + diff[1]**2.0 + diff[2]**2.0)


def calc_checksum(filename):
    try:
        with open(filename, 'rb') as file:
            return adler32(file.read())
    except FileNotFoundError:
        return adler32(bytes(filename, 'utf-8'))

def calculate_chunksizes(n_of_chunks, n_of_items):
    small_chunksize = n_of_items // n_of_chunks
    big_chunksize = small_chunksize + 1
    n_of_small_chunks = n_of_chunks * big_chunksize - n_of_items
    n_of_big_chunks = n_of_chunks - n_of_small_chunks
    if n_of_big_chunks == 0:
        big_chunksize = 0
    return small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks


def custom_encoder(obj):

    if isinstance(obj, set):
        return {'__set__': True, 'as_list': list(obj)}

    if 'Residue' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_residue = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'interaction_profile':
                serialized_residue.append(None)
            elif attribute_name == 'centralities':
                serialized_residue.append(None)
            else:
                serialized_residue.append(obj.__getattribute__(attribute_name))
        return {'__residue__': True, 'as_list': serialized_residue}

    if 'StructureAnnotation' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_annotation = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'alignment':
                serialized_annotation.append(None)
            else:
                serialized_annotation.append(obj.__getattribute__(attribute_name))
        return {'__structureannotation__': True, 'as_list': serialized_annotation}

    if 'Structure' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_structure = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'sequence':
                serialized_structure.append(None)
            else:
                serialized_structure.append(obj.__getattribute__(attribute_name))
        return {'__structure__': True, 'as_list': serialized_structure}

    if 'Complex' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_complex = []
        for attribute_name in obj.__slots__:
            #Complex objects get slimmed down, when packed !!!!!!!!!
            if attribute_name == 'chains' or attribute_name == 'resolution' or attribute_name == 'interfaces':
                serialized_complex.append(obj.__getattribute__(attribute_name))
            else:
                serialized_complex.append(None)
        return {'__complex__': True, 'as_list': serialized_complex}

    if 'Protein' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_protein = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'sequence':
                serialized_protein.append(None)
            else:
                serialized_protein.append(obj.__getattribute__(attribute_name))
        return {'__protein__': True, 'as_list': serialized_protein}

    if 'Position_Position_Interaction' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_pos_pos_i = []
        for attribute_name in obj.__slots__:
            serialized_pos_pos_i.append(obj.__getattribute__(attribute_name))
        return {'__pos_pos_i__': True, 'as_list': serialized_pos_pos_i}

    if 'Position' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_position = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'mappings':
                serialized_position.append(None)
            elif attribute_name == 'mut_aas':
                serialized_position.append(None)
            else:
                serialized_position.append(obj.__getattribute__(attribute_name))
        return {'__position__': True, 'as_list': serialized_position}

    if 'Aggregated_interface' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_agg_interface = []
        for attribute_name in obj.__slots__:
            serialized_agg_interface.append(obj.__getattribute__(attribute_name))
        return {'__aggregated_interface__': True, 'as_list': serialized_agg_interface}

    if 'Interface' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_interface = []
        for attribute_name in obj.__slots__:
            serialized_interface.append(obj.__getattribute__(attribute_name))
        return {'__interface__': True, 'as_list': serialized_interface}


    return obj

def custom_decoder(obj):

    if '__set__' in obj:
        return set(obj['as_list'])

    if '__residue__' in obj:
        serialized_residue = obj['as_list']
        res = sdsc.residue.Residue(0)
        for i, attribute_name in enumerate(res.__slots__):
            res.__setattr__(attribute_name, serialized_residue[i])
        return res

    if '__structureannotation__' in obj:
        serialized_annotation = obj['as_list']
        u_ac = serialized_annotation[0]
        pdb_id = serialized_annotation[1]
        chain = serialized_annotation[2]
        anno = sdsc.structure.StructureAnnotation(u_ac, pdb_id, chain)
        for i, attribute_name in enumerate(anno.__slots__[3:]):
            anno.__setattr__(attribute_name, serialized_annotation[i+3])
        return anno

    if '__structure__' in obj:
        serialized_structure = obj['as_list']
        pdb_id = serialized_structure[0]
        chain = serialized_structure[1]
        struct = sdsc.structure.Structure(pdb_id, chain)
        for i, attribute_name in enumerate(struct.__slots__[2:]):
            struct.__setattr__(attribute_name, serialized_structure[i+2])
        return struct

    if '__complex__' in obj:
        serialized_complex = obj['as_list']
        pdb_id = serialized_complex[0]
        compl = sdsc.complex.Complex(pdb_id)
        for i, attribute_name in enumerate(compl.__slots__[1:]):
            compl.__setattr__(attribute_name, serialized_complex[i+1])
        return compl

    if '__protein__' in obj:
        serialized_protein = obj['as_list']
        prot = sdsc.protein.Protein(None)
        for i, attribute_name in enumerate(prot.__slots__):
            prot.__setattr__(attribute_name, serialized_protein[i])
        return prot

    if '__pos_pos_i__' in obj:
        serialized_pos_pos_i = obj['as_list']
        pos_pos_i = sdsc.interface.Position_Position_Interaction(None, None, None, None)
        for i, attribute_name in enumerate(pos_pos_i.__slots__):
            pos_pos_i.__setattr__(attribute_name, serialized_pos_pos_i[i])
        return pos_pos_i

    if '__position__' in obj:
        serialized_position = obj['as_list']
        posi = sdsc.position.Position()
        for i, attribute_name in enumerate(posi.__slots__):
            try:
                posi.__setattr__(attribute_name, serialized_position[i])
            except:
                print('Warning: undefined slot in position object in custom_decoder: ', attribute_name, i, serialized_position)
                posi.__setattr__(attribute_name, None)
        return posi

    if '__aggregated_interface__' in obj:
        serialized_agg_interface = obj['as_list']
        agg_inter = sdsc.interface.Aggregated_interface(None)
        for i, attribute_name in enumerate(agg_inter.__slots__):
            agg_inter.__setattr__(attribute_name, serialized_agg_interface[i])
        return agg_inter

    if '__interface__' in obj:
        serialized_interface = obj['as_list']
        inter = sdsc.interface.Interface(None, None)
        for i, attribute_name in enumerate(inter.__slots__):
            inter.__setattr__(attribute_name, serialized_interface[i])
        return inter


    return obj

def serialize(item):
    """Serialization using MessagePack

    Args:
        item: Any item (int, str, list, dict) to serialize

    Returns: Serialized object

    """
    return msgpack.dumps(item)

def deserialize(item):
    """Deserialization using MessagePack

    Args:
        item: Any MessagePack serialized object

    Returns: Deserialized object

    """
    return msgpack.loads(item)


def compress(value):
    """Serialize and Zstandard compress the value

    Args:
        value: Any value, could be list, dict, string, int

    Returns: Compressed serialized bytes

    """
    return zstd.compress(value, 9, 2)

def decompress(compressed_value):
    """Zstandard decompress and deserialize the compressed value

    Args:
        compressed_value: Any bytes that was compressed using 'compress' function

    Returns: Decompressed and deserialized value if not None, else empty list as default value is returned

    """
    return zstd.decompress(compressed_value)

def pack(some_object):
    #packed_object = compress(pickletools.optimize(pickle.dumps(some_object, protocol = pickle.HIGHEST_PROTOCOL)))
    packed_object = compress(msgpack.packb(some_object, default = custom_encoder))
    return packed_object

def unpack(packed_object):
    #some_object = pickle.loads(decompress(packed_object))
    some_object = msgpack.unpackb(decompress(packed_object), object_hook = custom_decoder, strict_map_key = False, use_list = False)
    return some_object
