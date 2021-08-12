#!/usr/bin/env python3

import os
import canopus
from sys import argv


if __name__ == "__main__":
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
    hierarchy = ["kingdom", "superclass", "class", "subclass", "level 5",
                 "level 6", "level 7", "level 8"]
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
            # 4. save scores of most specific classes - todo: scores of all classes
            scores = []
            for c in chosen_classes_sorted:
                sc = compound.canopusfp[
                    C.sirius.statistics.workspace.revmap[c]]
                scores.append(sc)
            chosen_classes_sorted_sc = list(zip(chosen_classes_sorted, scores))
            # alternative: rank based on highest score instead of priority
            # chosen_class = [c for c, sc in zip(deepest_classifications, scores)
            #                 if sc == highest_sc][0]

            formula = compound.formula
            results.append(
                [node.componentId, node_id, formula, hierarchy[max_tree-1],
                 chosen_classes_sorted_sc])
    results.sort(key=lambda x: int(x[0]))
    header = ["GNPS_compnentindex", "Canopus_id", "Formula", "Deepest_lvl",
              "Classes_sorted"]
    print('\t'.join(header))
    for res in results:
        print(res)
