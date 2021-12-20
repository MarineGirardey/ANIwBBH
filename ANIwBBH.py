#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 10:25:53 2021
@author: Marine Girardey
"""

# Library imports
import csv


def get_feature_intervals_length(feature_t):
    """
    Function which extract the feature interval length in the feature table file to compute the coverage
    ---------
    :param feature_t:
    :return: feat_int_len_of_prdct_acc
    """
    feat_int_len_of_prdct_acc = {}
    with open(feature_t, newline='') as csv_file:
        spam_reader = csv.reader(csv_file, delimiter=';')
        for row in spam_reader:
            product_accession = row[10]
            feature_interval_length = row[-3]
            if len(product_accession) != 0:
                feat_int_len_of_prdct_acc[product_accession] = feature_interval_length

    feat_int_len_of_prdct_acc.pop('product_accession')
    return feat_int_len_of_prdct_acc


def identity_cov_filter_hit_table(len_dict, hit_tab):
    """
    Function which compute the coverage and store the hits only with identity and coverage > 70%
    ---------
    :param hit_tab:
    :param len_dict
    :return: best_cds_coverage
    """
    best_cds_coverage = {}

    print(hit_tab)
    with open(hit_tab, newline='') as csv_file:
        spam_reader = csv.reader(csv_file, delimiter=',')

        for row in spam_reader:
            query_accession = row[0].split('_')[2]
            target_accession = row[1].split('_')[2]
            alignment_length = row[3]
            identity = row[2]
            e_value = row[-2]
            for accession, cds_length in len_dict.items():
                if accession == query_accession:
                    coverage = (int(alignment_length) / int(cds_length)) * 100
                    if coverage > 70 and float(identity) > 70:
                        # Possible BBH (They are filtered)
                        best_cds_coverage[query_accession, target_accession] = e_value
    return best_cds_coverage


def bbh(best_cds_1, best_cds_2):
    """
    Function to get the BBH
    ---------
    :param best_cds_1:
    :param best_cds_2
    :return: bbh
    """
    bbh = {}
    for accessions_1, e_value_1 in best_cds_1.items():
        for accessions_2, e_value_2 in best_cds_2.items():
            if accessions_1[0] == accessions_2[1] and accessions_1[1] == accessions_2[0]:
                bbh[accessions_1] = e_value_1
    return bbh


def create_bbh_hit_table(list_bbh, hit_tab, new_file):
    """
    Function which create a file only with BBH
    ---------
    :param list_bbh:
    :param hit_tab
    :param new_file
    """
    new_opened_file = open(new_file, "a")
    with open(hit_tab, newline='') as csv_file:
        spam_reader = csv.reader(csv_file, delimiter=',')
        for row in spam_reader:
            accession_q = row[0].split('_')[2]
            accession_t = row[1].split('_')[2]
            for bbh, e_value in list_bbh.items():
                query_accession = bbh[0]
                target_accession = bbh[1]
                if query_accession in accession_q:
                    if target_accession in accession_t:
                        new_row = str(row)[1:-1]
                        new_opened_file.write(new_row + '\n')
    new_opened_file.close()


def get_bbh_gene_length(bbh_list, feat_tab):
    """
    Function to extract the BBH length to compute the ANI
    ---------
    :param bbh_list:
    :param feat_tab
    :return: bbh_length_dict
    """
    bbh_length_dict = {}

    with open(feat_tab, newline='') as csv_file:
        spam_reader = csv.reader(csv_file, delimiter=';')
        for row in spam_reader:
            for bbh, e_value in bbh_list.items():
                if bbh[0] == row[10] or bbh[1] == row[10]:
                    gene_length = int(row[8]) - int(row[7]) + 1
                    bbh_length_dict[bbh[0], bbh[1]] = gene_length

    return bbh_length_dict


def bbh_length_addition(bbh_len_dict):
    """
    Function to add all bbh length to compute the ANI
    ---------
    :param bbh_len_dict:
    :return: bbh_length_dict
    """
    bbh_len_addition = 0
    for val in bbh_len_dict.values():
        bbh_len_addition += val
    return bbh_len_addition


def best_bbh(bbh_hit_tab):
    """
    Function to get the best BH
    ---------
    :param bbh_hit_tab:
    :return: best_bbh
    """
    best_bbh = {}

    with open(bbh_hit_tab, newline='') as csv_file:
        spam_reader = csv.reader(csv_file, delimiter=',')
        for row in spam_reader:
            accession = row[0][2:-1], row[1][2:-1]
            e_value = float(row[-2][2:-1])
            percentage_identity = float(row[2][2:-1])
            alignment_length = float(row[3][2:-1])
            if accession in best_bbh:
                if e_value == 0.0 or e_value < best_bbh[accession][0]:
                    best_bbh[accession] = e_value, percentage_identity, alignment_length
                else:
                    pass
            else:
                best_bbh[accession] = e_value, percentage_identity, alignment_length

    return best_bbh


def compute_ani(best_bbh, bbh_len_addition):
    """
    Function to compute the ANI
    ---------
    :param best_bbh:
    :param bbh_len_addition:
    :return: ani
    """
    addition_bbh_multiplication = 0

    for accession, value in best_bbh.items():
        percentage_identity = value[1]
        alignment_length = value[2]
        identity_alignment_multiplication = percentage_identity * alignment_length
        addition_bbh_multiplication += identity_alignment_multiplication

    if bbh_len_addition != 0:
        ani = addition_bbh_multiplication / bbh_len_addition
        return ani
    else:
        return "0 BBH found"


if __name__ == "__main__":

    global_path = '/Users/marinegirardey/Documents/m2/GECO/ANIwBBH/'

    hit_and_feature_table_list = [['g37_vs_M2321', 'M2321_vs_g37', 'g37', 'M2321'],
                                  ['g37_vs_girerdii', 'girerdii_vs_g37', 'g37', 'girerdii'],
                                  ['g37_vs_M6320', 'M6320_vs_g37', 'g37', 'M6320'],
                                  ['g37_vs_NEAQ', 'NEAQ_vs_g37', 'g37', 'NEAQ'],
                                  ['girerdii_vs_M2321', 'M2321_vs_girerdii', 'girerdii', 'M2321'],
                                  ['girerdii_vs_M6320', 'M6320_vs_girerdii', 'girerdii', 'M6320'],
                                  ['girerdii_vs_NEAQ', 'NEAQ_vs_girerdii', 'girerdii', 'NEAQ'],
                                  ['M2321_vs_M6320', 'M6320_vs_M2321', 'M2321', 'M6320'],
                                  ['M2321_vs_NEAQ', 'NEAQ_vs_M2321', 'M2321', 'NEAQ'],
                                  ['M6320_vs_NEAQ', 'NEAQ_vs_M6320', 'M6320', 'NEAQ']]

    for genome_pairs in hit_and_feature_table_list:
        hit_table_1 = global_path + 'hit_table/' + genome_pairs[0] + '_alignment_hit_table.csv'
        hit_table_2 = global_path + 'hit_table/' + genome_pairs[1] + '_alignment_hit_table.csv'
        feature_table_1 = global_path + 'feature_table/' + genome_pairs[2] + '_feature_table.csv'
        feature_table_2 = global_path + 'feature_table/' + genome_pairs[3] + '_feature_table.csv'
        bbh_hit_table = global_path + 'bbh_file/' + genome_pairs[2] + '_' + genome_pairs[3] + '_bbh_hit_table.csv'

        dict_feat_int_len_of_prdct_acc_1 = get_feature_intervals_length(feature_table_1)
        dict_feat_int_len_of_prdct_acc_2 = get_feature_intervals_length(feature_table_2)
        best_cds_dict_1 = identity_cov_filter_hit_table(dict_feat_int_len_of_prdct_acc_1, hit_table_1)
        best_cds_dict_2 = identity_cov_filter_hit_table(dict_feat_int_len_of_prdct_acc_2, hit_table_2)

        list_of_bbh = bbh(best_cds_dict_1, best_cds_dict_2)

        create_bbh_hit_table(list_of_bbh, hit_table_1, bbh_hit_table)

        bbh_length_dict = get_bbh_gene_length(list_of_bbh, feature_table_1)
        bbh_length_add = bbh_length_addition(bbh_length_dict)
        b_bbh = best_bbh(bbh_hit_table)
        ani = compute_ani(b_bbh, bbh_length_add)
        print('-----START-----')
        print('For ' + genome_pairs[2] + ' and ' + genome_pairs[3] + ' genomes')
        print('ANI : ', ani)
        print('---------------')
