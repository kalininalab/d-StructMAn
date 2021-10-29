from pathlib import Path

import structman
import structman.lib as lib
from structman.lib import rinerator
from structman.lib import output
from structman import scripts


ROOT_DIR = Path(structman.__file__).parent.absolute()
LIB_DIR = Path(lib.__file__).parent.absolute()
RINERATOR_DIR = Path(rinerator.__file__).parent.absolute()
OUTPUT_DIR = Path(output.__file__).parent.absolute()
SCRIPTS_DIR = Path(scripts.__file__).parent.absolute()

TOOLS_DIR = f'{ROOT_DIR}/tools'

RINERATOR_PATH = f'{RINERATOR_DIR}/get_chains.py'

UNIPROT_DB_SQL = f'{SCRIPTS_DIR}/struct_man_db_uniprot.sql'
STRUCTMAN_DB_SQL = f'{SCRIPTS_DIR}/database_structure.sql'
REDUCE_HET_DICT = f'{RINERATOR_DIR}/reduce_wwPDB_het_dict.txt'
PDB_SYNC_SCRIPT = f'{SCRIPTS_DIR}/pdb-rsync.sh'
