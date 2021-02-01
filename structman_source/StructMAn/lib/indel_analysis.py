import multiprocessing

def separate_structure_annotations(wt_structure_annotations,mut_structure_annotations,config):
    sep_wt = {}
    sep_mut = {}
    for struct_tuple in wt_structure_annotations:
        if struct_tuple in mut_structure_annotations:
            if wt_structure_annotations[struct_tuple].sequence_identity == None:
                config.errorlog.add_error('sequence identity was None for: %s' % str(struct_tuple))
                continue
            #This is especially important to seperate structures for substitution-type indels
            if wt_structure_annotations[struct_tuple].sequence_identity >= mut_structure_annotations[struct_tuple].sequence_identity:
                sep_wt[struct_tuple] = wt_structure_annotations[struct_tuple]
            else:
                sep_mut[struct_tuple] = mut_structure_annotations[struct_tuple]
            del mut_structure_annotations[struct_tuple]
        else:
            sep_wt[struct_tuple] = wt_structure_annotations[struct_tuple]

    for struct_tuple in mut_structure_annotations:
        sep_mut[struct_tuple] = mut_structure_annotations[struct_tuple]

    return sep_wt,sep_mut

def indel_analysis(config,input_queue,out_queue,lock):
    with lock:
        input_queue.put(None)
    while True:
        inp = input_queue.get()
        if inp == None:
            return
        indel_obj,wt_structure_annotations,mut_structure_annotations = inp

        #The structures has to be checked that they include the flanks
        wt_structure_annotations,mut_structure_annotations = indel_obj.inclusion_check(wt_structure_annotations,mut_structure_annotations)

        wt_structure_annotations,mut_structure_annotations = separate_structure_annotations(wt_structure_annotations,mut_structure_annotations,config)

        analysis_list = indel_obj.get_scenario(len(wt_structure_annotations),len(mut_structure_annotations))

        print(indel_obj.wt_prot,indel_obj.get_notation(),analysis_list,len(wt_structure_annotations),len(mut_structure_annotations))

def para_indel_analysis(proteins,config,manager,lock):
    input_queue = manager.Queue()
    out_queue = manager.Queue()

    for indel_notation in proteins.indels:
        indel_obj = proteins.indels[indel_notation]
        wt_structure_annotations = proteins.get_protein_structure_annotations(indel_obj.wt_prot)
        mut_structure_annotations = proteins.get_protein_structure_annotations(indel_obj.mut_prot)
        input_queue.put((indel_obj,wt_structure_annotations,mut_structure_annotations))

    indel_processes = config.annotation_processes
    if indel_processes >= len(proteins.indels):
        indel_processes = len(proteins.indels)

    processes = {}

    for i in range(1,indel_processes + 1):
        p = multiprocessing.Process(target=indel_analysis, args=(config,input_queue,out_queue,lock))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    return
