#!/usr/bin/env python3

import os
import canopus
import time
from sys import argv
from collections import defaultdict
from typing import List, Dict, Tuple, Union, DefaultDict


def analyse_canopus(sirius_folder: str, gnps_folder: str,
                    output_folder: str = './'):
    """Wrapper for analysing and combining canopus output with gnps mol network

    :param sirius_folder: directory containing sirius/canopus output
    :param gnps_folder: directory containing gnps molecular networking output
    :param output_folder: directory to write (preliminary) output to
    :return:
    """
    # make canopus object
    C = canopus.Canopus(sirius_file, gnps_file)

    # loop through the nodes in molecular network
    class_p_cutoff = 0.5
    # it will find all trees reaching at least the Xth lvl and trace back
    max_class_depth = None
    hierarchy = ["kingdom", "superclass", "class", "subclass"] + \
                [f"level {i}" for i in range(5, 12)]
    results = get_classes_for_mol_network(C, hierarchy, class_p_cutoff,
                                          max_class_depth)
    write_classes_cluster_index(results, hierarchy, output_folder)

    # group classes per MF if fraction (count class/total spec in MF) above X
    mf_fraction_cutoff = 0.3
    #todo: treat -1 differently
    for comp_ind, cluster_ind_results in results.items():
        print('\n',comp_ind)
        num_cluster_inds = len(cluster_ind_results)
        print(num_cluster_inds)

        # 1. count the instances of each class in the componentindex (MF)
        h_counters = {h:defaultdict(int) for h in hierarchy}
        for cluster_ind_res in cluster_ind_results:
            ci_classes = cluster_ind_res[2]
            for h in hierarchy:
                for class_tup in ci_classes[h]:
                    c_name = class_tup[0]
                    if c_name:
                        h_counters[h][c_name] += 1
        print(h_counters)
        # 2. calculate fraction
        h_fraction = {}
        for h_lvl, h_dict in h_counters.items():
            cls_frac_list = []
            cls_names_set = set()
            for cls, count in h_dict.items():
                if cls:
                    frac = count / num_cluster_inds
                    if frac >= mf_fraction_cutoff:
                        cls_frac_list.append((cls, frac))
                        cls_names_set.add(cls.name)
            if not cls_frac_list:
                cls_frac_list = [('', None)]

            # 3. order each level on highest priority (above cutoff)
            priority_names = [c.name for c in C.sirius.statistics.priority]
            cls_frac_sorted = []
            for pr in priority_names:  # highest priority
                if pr in cls_names_set:
                    chosen_class_tup = [c_tup for c_tup in cls_frac_list
                                    if c_tup[0].name == pr][0]
                    cls_frac_sorted.append(chosen_class_tup)
            h_fraction[h_lvl] = cls_frac_sorted
        print(h_fraction)


def get_classes_for_mol_network(C: canopus.Canopus,
                                hierarchy: List[str],
                                class_p_cutoff: float,
                                max_class_depth: int)\
        -> DefaultDict[str, List[
            Union[str, Dict[str, List[Tuple[Union[str, float, None]]]]]]]:
    """Loop through mol network and gather CF and NPC classes

    :param C: Canopus object of canopus results with gnps mol network data
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param class_p_cutoff: probability cutoff for including a class
    :param max_class_depth: max class depth for finding CF class
    :return: classes output - dict of lists of {componentindex: [cluster index,
        formula, {CF_level: [(class, prob)]}, {NPC_level: [(class, prob)]}]}

    CF classes are found by looking for the class at deepest depth (or
    max_class_depth) and then ordering these deepest classes based on priority.
    Then, the classes are traced back to higher hierarchy and sorted in output,
    again based on priority of deepest classes.
    """
    results = defaultdict(list)
    for node_id, node in C.gnps.nodes.items():
        compound = C.sirius.compounds.get(node_id)  # get canopus compound obj
        if compound:
            cf_classes_dict = get_CF_classes(C, compound, hierarchy,
                                             class_p_cutoff, max_class_depth)

            formula = compound.formula
            results[node.componentId].append(
                [node_id, formula, cf_classes_dict])
    return results


def get_CF_classes(C: canopus.Canopus,
                   compound: canopus.ontology.Category,
                   hierarchy: List[str],
                   class_p_cutoff: float = 0.5,
                   max_class_depth: int = 5) ->\
        Dict[str, List[Tuple[Union[str, float, None]]]]:
    """Get the ClassyFire classes for each compound

    :param C: Canopus object of canopus results with gnps mol network data
    :param compound: object for the current compound
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param class_p_cutoff: probability cutoff for including a class
    :param max_class_depth: max class depth for finding CF class
    :return: dict with hierarchy as keys and list of tuples as values where
        each tuple contains a class (ontology.Category) and its score. List of
        classes is sorted on priority

    CF classes are found by looking for the class at deepest depth (or
    max_class_depth) and then ordering these deepest classes based on priority.
    Then, the classes are traced back to higher hierarchy and sorted in output,
    again based on priority of deepest classes.
    """
    # assignments with assign_most_specific_classes()
    # does not work by finding the most specific class, it finds the
    # class (regardless of level) with highest priority.
    # here: find deepest/most specific class with highest priority.
    # 1. find all classification trees above 0.5 (dict)
    classifications = [c for c in C.sirius.statistics.categoriesFor(
        compound, class_p_cutoff)]

    # 2.a. take all classifications above cutoff
    if not max_class_depth:
        deepest_classifications = classifications
        deepest_names = set(c.name for c in deepest_classifications)
    else:
    # 2.b. take deepest classification - find the most specific classes
    # classyFireGenus() gives a dict of the class trees
        max_tree = max(len(c.classyFireGenus()) for c in classifications)
        max_tree = min([max_tree, max_class_depth])
        deepest_classifications = [c for c in classifications if
                                   len(c.classyFireGenus()) >= max_tree]
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
            # h_lvl_str = ''
            h_lvl_elem = ('', None)
            if h_lvl_class:
                sc = compound.canopusfp[
                    C.sirius.statistics.workspace.revmap[h_lvl_class]]
                # h_lvl_str = f"{h_lvl_class}:{sc:.3f}"
                h_lvl_elem = (h_lvl_class, sc)
            # h_lvl_result.append(h_lvl_str)
            h_lvl_result.append(h_lvl_elem)
        # classes_dict_sc[h_lvl] = '; '.join(h_lvl_result)
        classes_dict_sc[h_lvl] = h_lvl_result
    return classes_dict_sc


def write_classes_cluster_index(results: DefaultDict[str, List[
            Union[str, Dict[str, List[Tuple[Union[str, float, None]]]]]]],
                                hierarchy: List[str],
                                output_folder: str = './'):
    """Write class results for each cluster index to file grouped by components

    :param results: dict of lists of {componentindex: [cluster index,
        formula, {CF_level: [(class, prob)]}, {NPC_level: [(class, prob)]}]}
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param output_folder: directory to write (preliminary) output to
    :return: None
    """
    output_file = os.path.join(output_folder,
                               "cluster_index_classifications.txt")
    header = ["compnentindex", "cluster index", "Formula"] + hierarchy
    with open(output_file, 'w') as outf:
        outf.write("{}\n".format('\t'.join(header)))
        for comp_ind, clust_ind_results in sorted(
                results.items(), key=lambda x: int(x[0])):
            for clust_ind_res in clust_ind_results:
                # turn CF classifications into strings
                cf_list = []
                for h in hierarchy:
                    h_str = []
                    for cl_tup in clust_ind_res[2][h]:
                        if cl_tup[1]:
                            h_str.append(f"{cl_tup[0]}:{cl_tup[1]:.3f}")
                        else:
                            h_str.append('')
                    cf_list.append('; '.join(h_str))
                # add everything to a list
                res_list = [comp_ind, clust_ind_res[0], clust_ind_res[1]] +\
                           cf_list
                res_str = '\t'.join(res_list)
                outf.write(f"{res_str}\n")

if __name__ == "__main__":
    start = time.time()
    print("Start")
    if len(argv) < 3:
        raise FileNotFoundError("\nUsage: python classification_to_gnps.py "+
            "sirius_folder gnps_folder output_folder(default: ./)")
    sirius_file = argv[1]  # "C:\Users\joris\Documents\iOMEGA_no_backup\NPLinker_results\sirius_test_spec"
    gnps_file = argv[2]  # "C:\Users\joris\Documents\iOMEGA_no_backup\NPLinker_results\ProteoSAFe-METABOLOMICS-SNETS-V2-eea136cb-download_clustered_spectra"
    output_dir = './'
    if len(argv) == 4:
        output_dir = argv[3]

    analyse_canopus(sirius_file, gnps_file, output_dir)

    end = time.time()
    print(f"\nTime elapsed: {end-start:.2f}s")