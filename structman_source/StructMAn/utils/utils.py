import math
import os
import psutil
from zlib import adler32

import ray
from structman import settings
from structman.lib.sdsc.residue import Residue
from structman.lib.sdsc.structure import StructureAnnotation, Structure
from structman.lib.sdsc.complex import Complex
from structman.lib.sdsc.protein import Protein
from structman.lib.sdsc.position import Position

import subprocess

import zstd
import msgpack
import pickle
import pickletools

def ray_init(config):
    if ray.is_initialized():
        return
    os.environ["PYTHONPATH"] = f'{settings.ROOT_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.LIB_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.RINERATOR_DIR}:{os.environ.get("PYTHONPATH", "")}'
    if config.iupred_path != '':
        os.environ["PYTHONPATH"] = f'{os.path.abspath(os.path.realpath(config.iupred_path))}:{os.environ.get("PYTHONPATH", "")}'
    logging_level = 20
    if config.verbosity <= 1:
        logging_level = 0
    ray.init(num_cpus=config.proc_n, include_dashboard=False, ignore_reinit_error=True, logging_level = logging_level)

def ray_hack():
    # hack proposed by the devs of ray to prevent too many processes being spawned
    return
    #resources = ray.ray.get_resource_ids()
    #cpus = [v[0] for v in resources['CPU']]
    #psutil.Process().cpu_affinity(cpus)


def monitor_ray_store():
    p = subprocess.Popen(['ray','memory'])

    page, err = p.communicate()
    print(page)

    return

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
            if attribute_name == 'chains' or attribute_name == 'resolution':
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

    return obj

def custom_decoder(obj):
    if '__set__' in obj:
        return set(obj['as_list'])

    if '__residue__' in obj:
        serialized_residue = obj['as_list']
        res = Residue(0)
        for i, attribute_name in enumerate(res.__slots__):
            res.__setattr__(attribute_name, serialized_residue[i])
        return res

    if '__structureannotation__' in obj:
        serialized_annotation = obj['as_list']
        u_ac = serialized_annotation[0]
        pdb_id = serialized_annotation[1]
        chain = serialized_annotation[2]
        anno = StructureAnnotation(u_ac, pdb_id, chain)
        for i, attribute_name in enumerate(anno.__slots__[3:]):
            anno.__setattr__(attribute_name, serialized_annotation[i+3])
        return anno

    if '__structure__' in obj:
        serialized_structure = obj['as_list']
        pdb_id = serialized_structure[0]
        chain = serialized_structure[1]
        struct = Structure(pdb_id, chain)
        for i, attribute_name in enumerate(struct.__slots__[2:]):
            struct.__setattr__(attribute_name, serialized_structure[i+2])
        return struct

    if '__complex__' in obj:
        serialized_complex = obj['as_list']
        pdb_id = serialized_complex[0]
        compl = Complex(pdb_id)
        for i, attribute_name in enumerate(compl.__slots__[1:]):
            compl.__setattr__(attribute_name, serialized_complex[i+1])
        return compl

    if '__protein__' in obj:
        serialized_protein = obj['as_list']
        prot = Protein(None)
        for i, attribute_name in enumerate(prot.__slots__):
            prot.__setattr__(attribute_name, serialized_protein[i])
        return prot

    if '__position__' in obj:
        serialized_position = obj['as_list']
        posi = Position()
        for i, attribute_name in enumerate(posi.__slots__):
            posi.__setattr__(attribute_name, serialized_position[i])
        return posi

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
