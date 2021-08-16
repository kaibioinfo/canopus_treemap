#!/usr/bin/env python3

import os
import canopus
import time
from sys import argv
from collections import defaultdict


if __name__ == "__main__":
    start = time.time()
    print("Start")
    if len(argv) < 3:
        raise FileNotFoundError("\nUsage: python classification_to_gnps.py "+
                                "sirius_folder gnps_folder")
    sirius_file = argv[1]  # "C:\Users\joris\Documents\iOMEGA_no_backup\NPLinker_results\sirius_test_spec"
    gnps_file = argv[2]  # "C:\Users\joris\Documents\iOMEGA_no_backup\NPLinker_results\ProteoSAFe-METABOLOMICS-SNETS-V2-eea136cb-download_clustered_spectra"

    # make canopus object
    C = canopus.Canopus(sirius_file, gnps_file)

    # loop through the nodes in molecular network
    class_p_cutoff = 0.5
    results = []
    hierarchy = ["kingdom", "superclass", "class", "subclass", "level 5"]
    for node_id, node in C.gnps.nodes.items():
        compound = C.sirius.compounds.get(node_id)  # get canopus compound obj
        if compound:
            # assignments with assign_most_specific_classes()
            # does not work by finding the most specific class, it finds the
            # class (regardless of level) with highest priority.
            # here: find deepest/most specific class with highest priority.
            # 1. find all classification trees above 0.5 (dict)
            classifications = [c for c in C.sirius.statistics.categoriesFor(
                compound, class_p_cutoff)]

            # 2. take deepest classification - find the most specific classes
            # classyFireGenus() gives a dict of the class trees
            max_tree = max(len(c.classyFireGenus()) for c in classifications)
            deepest_classifications = [c for c in classifications if
                             len(c.classyFireGenus()) == max_tree]
            deepest_names = set(c.name for c in deepest_classifications)

            # 3. choose the most specific class with top priority
            priority_names = [c.name for c in C.sirius.statistics.priority]
            chosen_classes_sorted = []
            for pr in priority_names:  # find deepest with highest priority
                if pr in deepest_names:
                    chosen_class = [c for c in deepest_classifications
                                    if c.name == pr][0]
                    chosen_classes_sorted.append(chosen_class)

            # 4. unpack all classes
            classes_dict = defaultdict(list)
            classes_set_dict = defaultdict(set)
            for h_lvl in hierarchy:
                for c in chosen_classes_sorted:
                    c_h_lvl = c.classyFireGenus().get(h_lvl, '')
                    c_h_set = classes_set_dict[h_lvl]
                    if c_h_lvl not in c_h_set:
                        classes_dict[h_lvl].append(c_h_lvl)
                        c_h_set.add(c_h_lvl)
                        if not c_h_lvl:
                            # there is no current level - go to next
                            break

            # 5. get all scores and convert ontology.Category to strings
            classes_dict_sc = {}
            for h_lvl in hierarchy:
                h_lvl_classes = classes_dict[h_lvl]
                h_lvl_result = []
                for h_lvl_class in h_lvl_classes:
                    h_lvl_str = ''
                    if h_lvl_class:
                        sc = compound.canopusfp[
                            C.sirius.statistics.workspace.revmap[h_lvl_class]]
                        h_lvl_str = f"{h_lvl_class}:{sc:.3f}"
                    h_lvl_result.append(h_lvl_str)
                classes_dict_sc[h_lvl] = '; '.join(h_lvl_result)

            formula = compound.formula
            results.append(
                [node.componentId, node_id, formula, classes_dict_sc])

    results.sort(key=lambda x: (int(x[0]), int(x[1])))
    header = ["GNPS_compnentindex", "Canopus_id", "Formula",
              "Classes_sorted"]
    print('\t'.join(header))
    for res in results:
        print(res)

    end = time.time()
    print(f"\nTime elapsed: {end-start:.2f}s")