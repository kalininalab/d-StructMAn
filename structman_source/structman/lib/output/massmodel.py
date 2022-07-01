import os
import shutil

import ray

from structman.lib.database import database
try:
    from structman.lib import modelling
except:
    pass

from structman.lib.output import out_generator

def get_optimal_templates(prot_id, proteins, config):

    best_interface_structures = {}
    for interface_number, (positions, interface_structure_map) in enumerate(proteins[prot_id].aggregated_interface_map):
        best_score = 0.
        best_interface_structure = None
        for structure_id in interface_structure_map:
            chain, i_chain, interface_size = interface_structure_map[structure_id]
            seq_id = proteins[prot_id].structure_annotations[(structure_id, chain)].sequence_identity
            cov = proteins[prot_id].structure_annotations[(structure_id, chain)].coverage
            score = interface_size * seq_id * cov
            if score > best_score:
                best_score = score
                best_interface_structure = (structure_id, chain)
        if not best_interface_structure in best_interface_structures:
            best_interface_structures[best_interface_structure] = [interface_number]
        else:
            best_interface_structures[best_interface_structure].append(interface_number)

    for position in proteins[prot_id].positions:
        for (structure_id, tchain) in proteins[prot_id].structure_annotations:
            sub_info = protein_obj.structure_annotations[(pdb_id, chain)].get_sub_info(pos)
            res_nr = sub_info[0]
    #TODO

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
                print('No rec struct:', prot_id, wildtype_protein)
                continue
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

    model_output = out_generator.OutputGenerator()

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
