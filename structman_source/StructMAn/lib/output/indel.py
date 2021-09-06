import os

from structman.lib import database
from structman.lib.output.output import OutputGenerator


def create_indel_results_table(config, output_path, session_name, session_id):
    table = 'RS_Indel_Session'
    rows = ['Indel', 'Tags']
    eq_rows = {'Session': session_id}

    results = database.select(config, rows, table, equals_rows=eq_rows)

    tag_map = {}
    ids = []
    for row in results:
        indel_id = row[0]
        tags = row[1]
        tag_map[indel_id] = tags
        ids.append(indel_id)

    cols = ['Indel_Id', 'Indel_Notation', 'Delta_Delta_Classification']
    results = database.binningSelect(ids, cols, 'Indel', config)

    indel_output = OutputGenerator()

    headers = ['Indel', 'Tags', 'Delta delta classification']

    indel_output.add_headers(headers)

    indel_file = '%s/%s_indel_analysis.tsv' % (output_path, session_name)

    if os.path.exists(indel_file):
        os.remove(indel_file)

    f = open(indel_file, 'a')
    f.write(indel_output.get_header())

    for row in results:
        indel_id = row[0]
        indel_output.add_value('Indel', row[1])
        indel_output.add_value('Tags', tag_map[indel_id])
        indel_output.add_value('Delta delta classification', row[2])

        f.write(indel_output.pop_line())

    f.close()
    return
