import os
import shutil

import ray

from structman.lib import database, modelling
from structman.lib.output.output import OutputGenerator


# called by structman_main
def mass_model(session_id, config, outfolder, include_snvs=True):
    proteins = database.proteinsFromDb(session_id, config, with_residues=True,
                                       with_snvs=True, mutate_snvs=include_snvs, with_alignments=True,
                                       with_complexes=True, keep_ray_alive=True)
    prot_ids = proteins.get_protein_ids()

    modelling_result_ids = []

    conf_dump = ray.put(config)

    for prot_id in prot_ids:
        rec_struct = proteins[prot_id].get_major_recommend_structure()
        if rec_struct is None:
            wildtype_protein = proteins[prot_id].wildtype_protein
            if wildtype_protein is None:
                print('Error wt prot is None:', prot_id)
            rec_struct = proteins[wildtype_protein].get_major_recommend_structure()
            if rec_struct is None:
                print('Error with rec struct:', prot_id, wildtype_protein)
        pdb_id, tchain = rec_struct.split(':')
        compl_obj = proteins.complexes[pdb_id]
        structures = proteins.get_complex_structures(pdb_id)
        print('Modelling of:', prot_id, pdb_id, tchain)

        alignment_tuple = proteins.get_alignment(prot_id, pdb_id, tchain)

        seq_id = proteins[prot_id].structure_annotations[(pdb_id, tchain)].sequence_identity
        cov = proteins[prot_id].structure_annotations[(pdb_id, tchain)].coverage
        modelling_result_ids.append(modelling.model.remote(conf_dump, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, skip_analysis=True))

    models = ray.get(modelling_result_ids)

    if len(models) == 0:
        print('No models got produced :(')
        return

    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    model_database = '%s/models' % outfolder
    summary_file = '%s/model_summary.tsv' % outfolder

    if not os.path.exists(model_database):
        os.mkdir(model_database)

    if os.path.exists(summary_file):
        os.remove(summary_file)

    model_output = OutputGenerator()

    headers = ['Input ID', 'Wildtype protein', 'Protein ID', 'Template PDB ID', 'Template chain',
               'Model chain', 'Template resolution', 'Sequence identity', 'Coverage', 'File location']
    model_output.add_headers(headers)
    f = open(summary_file, 'a')
    f.write(model_output.get_header())

    for model in models:
        if isinstance(model, str):
            config.errorlog.add_warning(model)
        else:
            target_path = '%s/%s' % (model_database, model.path.split('/')[-1])
            shutil.copy(model.path, target_path)
            # model.clear_tmp()
            pdb_id, tchain = model.template_structure
            model_chain = model.chain_id_map[tchain]

            input_id = proteins[model.target_protein].input_id
            if input_id is None:
                input_id = proteins[proteins[model.target_protein].wildtype_protein].input_id

            model_output.add_value('Input ID', input_id)
            model_output.add_value('Wildtype protein', proteins[model.target_protein].wildtype_protein)
            model_output.add_value('Protein ID', model.target_protein)
            model_output.add_value('Template PDB ID', pdb_id)
            model_output.add_value('Template chain', tchain)
            model_output.add_value('Model chain', model_chain)
            model_output.add_value('Template resolution', model.template_resolution)
            model_output.add_value('Sequence identity', model.sequence_identity)
            model_output.add_value('Coverage', model.coverage)
            model_output.add_value('File location', target_path)

            f.write(model_output.pop_line())

    f.close()
    ray.shutdown()
