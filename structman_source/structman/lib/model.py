import shutil
import os
from structman.lib import templateFiltering


class Model:
    __slots__ = ['path', 'template_structure', 'target_protein', 'model_id', 'structural_analysis_dict', 'ligand_profiles',
                 'metal_profiles', 'ion_profiles', 'chain_chain_profiles', 'chain_type_map', 'chainlist', 'chain_id_map',
                 'template_resolution', 'sequence_identity', 'coverage', 'tmp_folder']

    def __init__(self, path='', template_structure=None, target_protein=None, model_id=None, chain_id_map=None,
                 template_resolution=None, sequence_identity=None, coverage=None, tmp_folder=None):
        self.path = path
        self.template_structure = template_structure  # tuple of template pdb id and target chain
        self.target_protein = target_protein
        self.model_id = model_id
        self.chain_id_map = chain_id_map  # maps the template chain ids to the model chain ids
        self.template_resolution = template_resolution
        self.sequence_identity = sequence_identity
        self.coverage = coverage
        self.tmp_folder = tmp_folder

    def analyze(self, config):
        config.n_of_chain_thresh = 1000  # Hinder structuralAnalysis to spawn processes, since this function is already called by a remote
        model_target_chain = self.chain_id_map[self.template_structure[1]]

        if config.verbosity >= 3:
            print('Start model self analysis:', self.model_id, self.path)

        if self.path[-12:] != '_refined.pdb':
            self.refine_model()

        (structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, nested_processes) = templateFiltering.structuralAnalysis(self.model_id, config, model_path=self.path, target_dict=[model_target_chain], keep_rin_files=True)

        self.structural_analysis_dict = structural_analysis_dict
        self.ligand_profiles = ligand_profiles
        self.metal_profiles = metal_profiles
        self.ion_profiles = ion_profiles
        self.chain_chain_profiles = chain_chain_profiles
        self.chain_type_map = chain_type_map
        self.chainlist = chainlist

        for error_text in errorlist:
            config.errorlog.add_warning(error_text)

    def clear_tmp(self):
        if self.tmp_folder is None:
            return
        try:
            shutil.rmtree(self.tmp_folder)
        except:
            print('Could not remove tmp folder:', self.tmp_folder)
        return

    # modeller increases the residue ID continuously even when a new chain begun, this can lead to residue IDs > 9999
    # modeller solves this by using letters, for example: A000 = 10000
    def refine_model(self):
        f = open(self.path, 'r')
        lines = f.readlines()
        f.close()

        new_lines = []
        current_chain = None
        current_res = None
        for line in lines:
            if len(line) >= 21:
                record_name = line[0:6].rstrip()
                if record_name == "ATOM" or record_name == 'HETATM':
                    chain_id = line[21]
                    res_nr = line[22:26].strip()  # without insertion code
                    # reset res_id counter for every new chain
                    if current_chain != chain_id:
                        current_new_res_nr = 0
                        current_chain = chain_id
                    if current_res != res_nr:
                        current_res = res_nr
                        current_new_res_nr += 1
                        digit_res_str = str(current_new_res_nr)
                        current_res_str = '%s%s' % (' ' * (4 - len(digit_res_str)), digit_res_str)
                    new_lines.append('%s%s%s' % (line[:22], current_res_str, line[26:]))

                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)

        refined_path = '%s_refined.pdb' % self.path[:-4].replace('.', '_')

        f = open(refined_path, 'w')
        f.write(''.join(new_lines))
        f.close()

        os.remove(self.path)
        self.path = refined_path
        return
