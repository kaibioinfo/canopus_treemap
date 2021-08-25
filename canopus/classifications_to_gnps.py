#!/usr/bin/env python3

import os
import canopus
import time
from sys import argv
from collections import defaultdict
from typing import List, Dict, Tuple, Union, DefaultDict, Any


def analyse_canopus(sirius_folder: str, gnps_folder: str,
                    output_folder: str = './'):
    """Wrapper for analysing and combining canopus output with gnps mol network

    :param sirius_folder: directory containing sirius/canopus output
    :param gnps_folder: directory containing gnps molecular networking output
    :param output_folder: directory to write (preliminary) output to
    :return:
    """
    # make canopus object
    can = canopus.Canopus(sirius_folder, gnps_folder)

    # find classes per cluster index (spectra): loop through the nodes in
    # molecular network
    class_p_cutoff = 0.5
    # it will use all classes (set max_class_depth to a number for only keeping
    # the deepest classes at the Xth level)
    max_class_depth = None
    hierarchy = ["kingdom", "superclass", "class", "subclass"] + \
                [f"level {i}" for i in range(5, 12)]
    npc_hierarchy = ['pathway', 'superclass', 'class']
    class_results = get_classes_for_mol_network(
        can, hierarchy, npc_hierarchy, class_p_cutoff, max_class_depth)

    # write cluster index (spectra) results
    write_classes_cluster_index(class_results, hierarchy, npc_hierarchy,
                                output_folder)

    # group classes per componentindex (molecular family)
    mf_fraction_cutoff = 0.2
    comp_ind_classes = get_classes_for_componentindices(
        class_results, can, hierarchy, npc_hierarchy, mf_fraction_cutoff,
        mf_fraction_cutoff)

    # write componentindex results
    write_classes_componentindex(
        comp_ind_classes, hierarchy, npc_hierarchy, output_folder)


def get_classes_for_mol_network(can: canopus.Canopus,
                                hierarchy: List[str],
                                npc_hierarchy: List[str],
                                class_p_cutoff: float,
                                max_class_depth: Union[int, None]) -> \
        DefaultDict[str, List[Union[str, Dict[str, List[Tuple[
            Union[str, float]]]]]]]:
    """Loop through mol network and gather CF and NPC classes

    :param can: Canopus object of canopus results with gnps mol network data
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param npc_hierarchy: the NPC class level names to be included in output in
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
    for node_id, node in can.gnps.nodes.items():
        # get canopus compound obj
        compound = can.sirius.compounds.get(node_id)
        if compound:
            cf_classes_dict = get_cf_classes(can, compound, hierarchy,
                                             class_p_cutoff, max_class_depth)
            npc_classes_dict = get_npc_classes(can, compound, npc_hierarchy)
            formula = compound.formula
            comp_id = node.componentId
            if comp_id == '-1':  # handling of singleton -1 components
                comp_id += f"_{node_id}"
            results[comp_id].append(
                [node_id, formula, cf_classes_dict, npc_classes_dict])
    return results


def get_cf_classes(can: canopus.Canopus,
                   compound: canopus.ontology.Category,
                   hierarchy: List[str],
                   class_p_cutoff: float = 0.5,
                   max_class_depth: Union[int, None] = None
                   ) -> Dict[str, List[Tuple[Any, Any]]]:
    """Get the ClassyFire classes for a compound

    :param can: Canopus object of canopus results with gnps mol network data
    :param compound: object for the current compound
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param class_p_cutoff: probability cutoff for including a class
    :param max_class_depth: None for using all classes, set this to a number
        for only taking classes into account at a certain max depth
    :return: dict with hierarchy as keys and list of tuples as values where
        each tuple contains a class (ontology.Category) and its score. List of
        classes is sorted on priority

    This is kind of an elaboration of assign_most_specific_classes() but here
    it finds all classes and orders them on priority, and then it fills the
    result lists by starting with the highest priority class and tracing that
    class back to all parent classes, then it moves to the second priority
    class, etc.
    There is also an option to only take into account the class(es) at deepest
    depth (or max_class_depth) and then ordering these deepest classes based on
    priority, etc. For this max_class_depth needs to be set to a number.
    """
    # 1. find all classification trees above 0.5 (dict)
    classifications = [c for c in can.sirius.statistics.categoriesFor(
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
    priority_names = [c.name for c in can.sirius.statistics.priority]
    chosen_classes_sorted = []
    for pr in priority_names:  # find deepest with highest priority
        if pr in deepest_names:
            chosen_class = [c for c in deepest_classifications
                            if c.name == pr][0]
            chosen_classes_sorted.append(chosen_class)

    # 4. unpack/trace back all classes to parents even if parents have low prob
    classes_dict = defaultdict(list)
    classes_set_dict = defaultdict(set)
    for c in chosen_classes_sorted:
        for h_lvl in hierarchy:
            c_h_lvl = c.classyFireGenus().get(h_lvl)
            if c_h_lvl:
                c_h_set = classes_set_dict[h_lvl]
                if c_h_lvl not in c_h_set:
                    classes_dict[h_lvl].append(c_h_lvl)
                    c_h_set.add(c_h_lvl)

    # 5. get all scores to tuple(cls, sc) and fill empty levels with empty list
    classes_dict_sc = {}
    for h_lvl in hierarchy:
        h_lvl_classes = classes_dict[h_lvl]
        h_lvl_result = []
        for h_lvl_class in h_lvl_classes:
            sc = compound.canopusfp[
                can.sirius.statistics.workspace.revmap[h_lvl_class]]
            h_lvl_elem = (h_lvl_class, sc)
            h_lvl_result.append(h_lvl_elem)
        # classes_dict is default dict so _sc will be populated with empty
        # list if there are no classes at a certain level
        classes_dict_sc[h_lvl] = h_lvl_result
    return classes_dict_sc


def get_npc_classes(can: canopus.Canopus,
                    compound: canopus.ontology.Category,
                    npc_hierarchy: List[str]):
    """Get NPClassifier classes for a compound

    :param can: Canopus object of canopus results with gnps mol network data
    :param compound: object for the current compound
    :param npc_hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :return: dict
    """
    npc_array = can.npcSummary().loc[compound.name]
    npc_dict = {}
    for h_lvl in npc_hierarchy:
        cls = npc_array[h_lvl]
        npc_dict[h_lvl] = []  # init level in dict
        if cls != "N/A":  # add class if there is one
            prob_key = h_lvl + 'Probability'
            npc_dict[h_lvl].append((cls, npc_array[prob_key]))
    return npc_dict


def write_classes_cluster_index(results: DefaultDict[str, List[
    Union[str, Dict[str, List[Tuple[Union[str, float]]]]]]],
                                hierarchy: List[str],
                                npc_hierarchy: List[str],
                                output_folder: str = './'):
    """Write class results for each cluster index to file grouped by components

    :param results: dict of lists of {componentindex: [cluster index,
        formula, {CF_level: [(class, prob)]}, {NPC_level: [(class, prob)]}]}
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param npc_hierarchy: the NPC class level names to be included in output in
        order of hierarchy
    :param output_folder: directory to write (preliminary) output to
    :return: None
    """
    output_file = os.path.join(output_folder,
                               "cluster_index_classifications.txt")
    header = ["componentindex", "cluster index", "formula"] + hierarchy + \
        npc_hierarchy
    with open(output_file, 'w') as outf:
        outf.write("{}\n".format('\t'.join(header)))
        # sort on MF/componentindex and then on clusterindex taking care of -1
        for comp_ind, clust_ind_results in sorted(
                results.items(), key=lambda x: (
                        int(x[0].split("_")[0]))):
            for clust_ind_res in sorted(clust_ind_results,
                                        key=lambda x: int(x[0])):
                # turn CF classifications into strings
                cf_list = []
                cf_res = clust_ind_res[2]
                for h in hierarchy:
                    h_str = []
                    for cl_tup in cf_res[h]:
                        h_str.append(f"{cl_tup[0]}:{cl_tup[1]:.3f}")
                    cf_list.append('; '.join(h_str))
                # turn NPC classifications into strings - todo: make subroutine
                npc_list = []
                npc_res = clust_ind_res[3]
                for h in npc_hierarchy:
                    h_str = []
                    for cl_tup in npc_res[h]:
                        h_str.append(f"{cl_tup[0]}:{cl_tup[1]:.3f}")
                    npc_list.append('; '.join(h_str))
                # add everything to a list
                res_list = [comp_ind, clust_ind_res[0], clust_ind_res[1]] + \
                    cf_list + npc_list
                res_str = '\t'.join(res_list)
                outf.write(f"{res_str}\n")


def get_classes_for_componentindices(clusterindex_results,
                                     can,
                                     hierarchy: List[str],
                                     npc_hierarchy: List[str],
                                     cf_fraction_cutoff: float = 0.3,
                                     npc_fraction_cutoff: float = 0.2):
    result_list = []
    for comp_ind, cluster_ind_results in clusterindex_results.items():
        num_cluster_inds = len(cluster_ind_results)
        comp_ind_cf_classes_dict = get_cf_classes_for_componentindex(
            cluster_ind_results, num_cluster_inds, can, hierarchy,
            cf_fraction_cutoff)
        comp_ind_npc_classes_dict = get_classes_npc_for_componentindex(
            cluster_ind_results, num_cluster_inds, npc_hierarchy,
            npc_fraction_cutoff)
        result_list.append((
            comp_ind, num_cluster_inds, comp_ind_cf_classes_dict,
            comp_ind_npc_classes_dict))
    return result_list


def get_cf_classes_for_componentindex(cluster_ind_results,
                                      num_cluster_inds,
                                      can,
                                      hierarchy: List[str],
                                      fraction_cutoff: float = 0.3):
    comp_ind_scores = {}  # dict{class: fraction_score}

    # 1. count the instances of each class in the componentindex (MF)
    h_counters = {h: defaultdict(int) for h in hierarchy}
    for cluster_ind_res in cluster_ind_results:
        ci_classes = cluster_ind_res[2]
        for h in hierarchy:
            for class_tup in ci_classes[h]:
                c_name = class_tup[0]
                if c_name:
                    h_counters[h][c_name] += 1

    # 2. calculate fraction and save to dict, irregardless of hierarchy lvl
    for h_lvl, h_dict in h_counters.items():
        for cls, count in h_dict.items():
            if cls:
                frac = count / num_cluster_inds
                if frac >= fraction_cutoff:
                    comp_ind_scores[cls] = frac

    # 3. order all classes on priority
    priority_classes = [c for c in can.sirius.statistics.priority]
    comp_ind_sorted_pr = []
    for pr in priority_classes:  # highest priority
        if pr in comp_ind_scores:
            comp_ind_sorted_pr.append(pr)

    # 4. fill hierarchy based on priority
    comp_ind_classes_dict = defaultdict(list)
    comp_ind_classes_set_dict = defaultdict(set)
    for c in comp_ind_sorted_pr:
        for h_lvl in hierarchy:
            c_h_lvl = c.classyFireGenus().get(h_lvl)
            if c_h_lvl:
                c_h_set = comp_ind_classes_set_dict[h_lvl]
                # make sure to add a classification only once to lvl
                if c_h_lvl not in c_h_set:
                    comp_ind_sc = comp_ind_scores[c_h_lvl]
                    comp_ind_classes_dict[h_lvl].append(
                        (c_h_lvl, comp_ind_sc))
                    c_h_set.add(c_h_lvl)
    return comp_ind_classes_dict


def get_classes_npc_for_componentindices(clusterindex_results,
                                         npc_hierarchy: List[str],
                                         fraction_cutoff: float = 0.3):
    result_list = []
    for comp_ind, cluster_ind_results in clusterindex_results.items():
        num_cluster_inds = len(cluster_ind_results)
        scores_dict = get_classes_npc_for_componentindex(
            cluster_ind_results, num_cluster_inds, npc_hierarchy,
            fraction_cutoff)
        result_list.append((comp_ind, num_cluster_inds, scores_dict))
    return result_list


def get_classes_npc_for_componentindex(cluster_ind_results,
                                       num_cluster_inds,
                                       npc_hierarchy,
                                       fraction_cutoff: float = 0.2):
    # 1. count the instances of each class in the componentindex (MF)
    h_counters = {h: defaultdict(int) for h in npc_hierarchy}
    for cluster_ind_res in cluster_ind_results:
        npc_classes = cluster_ind_res[3]
        for h in npc_hierarchy:
            for class_tup in npc_classes[h]:
                c_name = class_tup[0]
                if c_name:
                    h_counters[h][c_name] += 1

    # 2. calculate fraction and add to dict
    scores_dict = defaultdict(list)
    for h_lvl, h_dict in h_counters.items():
        for cls, count in h_dict.items():
            if cls:
                frac = count / num_cluster_inds
                if frac >= fraction_cutoff:
                    # comp_ind_scores[cls] = frac
                    scores_dict[h_lvl].append((cls, frac))
    return scores_dict


def write_classes_componentindex(comp_ind_cf_classes,
                                 hierarchy, npc_hierarchy, output_folder):
    """
    Write component indices class results - component_index_classifications.txt

    :param comp_ind_cf_classes:
    :param hierarchy: the CF class level names to be included in output in
        order of hierarchy
    :param npc_hierarchy: the NPC class level names to be included in output in
        order of hierarchy
    :param output_folder: directory to write (preliminary) output to
    :return: None
    """
    output_file = os.path.join(output_folder,
                               "component_index_classifications.txt")
    header = ["componentindex", "size"] + hierarchy + \
        npc_hierarchy
    with open(output_file, 'w') as outf:
        outf.write("{}\n".format('\t'.join(header)))
        for comp_ind, comp_ind_len, cf_res, npc_res in sorted(
                comp_ind_cf_classes, key=lambda x: int(x[0].split("_")[0])):
            # turn CF classifications into strings todo: make subroutine
            cf_list = []
            for h in hierarchy:
                h_str = []
                for cl_tup in cf_res[h]:
                    h_str.append(f"{cl_tup[0]}:{cl_tup[1]:.3f}")
                cf_list.append('; '.join(h_str))
            # turn NPC classifications into strings
            npc_list = []
            for h in npc_hierarchy:
                h_str = []
                for cl_tup in npc_res[h]:
                    h_str.append(f"{cl_tup[0]}:{cl_tup[1]:.3f}")
                npc_list.append('; '.join(h_str))
            # add everything to a list
            res_list = [comp_ind, str(comp_ind_len)] + \
                cf_list + npc_list
            res_str = '\t'.join(res_list)
            outf.write(f"{res_str}\n")


if __name__ == "__main__":
    start = time.time()
    print("Start")
    if len(argv) < 3:
        raise FileNotFoundError(
            "\nUsage: python classification_to_gnps.py " +
            "sirius_folder gnps_folder output_folder(default: ./)")
    sirius_file = argv[
        1]  # "C:\Users\joris\Documents\iOMEGA_no_backup\NPLinker_results\subset_canopus_crus_mgf_maxmz600"
    gnps_file = argv[
        2]  # "C:\Users\joris\Documents\iOMEGA_no_backup\NPLinker_results\ProteoSAFe-METABOLOMICS-SNETS-V2-eea136cb-download_clustered_spectra"
    output_dir = './'
    if len(argv) == 4:
        output_dir = argv[3]

    analyse_canopus(sirius_file, gnps_file, output_dir)

    end = time.time()
    print(f"\nTime elapsed: {end - start:.2f}s")
