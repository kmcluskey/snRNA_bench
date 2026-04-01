import os
import numpy as np
import pandas as pd
import json
from loguru import logger
from IPython.display import display
from pathlib import Path
from scipy.stats import kendalltau
from sklearn.preprocessing import RobustScaler, StandardScaler, MinMaxScaler

"""
A set of helper files that can be used throughout the various TubAtlas apps.
"""


def get_dros_genes_df():

    # Get the DF that holds the FBgn codes, annotations and synonyms for use in this helper class
    DATA_FOLDER = Path('/Users/Karen/TubAtlas/data/Orthologues/DrosophilaMelanogaster')
    dros_genes_df = pd.read_csv(os.path.join(DATA_FOLDER, 'dros_annotations_2023_05.csv'),
                                escapechar='\\').set_index('FBgn')

    # Replace the NaN with the gene name 'nan' for FBgn0036414
    dros_genes_df.loc['FBgn0036414', 'Symbol'] = 'nan'
    dros_genes_df.loc['FBgn0036414', 'Sec_Annot'] = 'nan'

    return dros_genes_df


def get_gene_code_df(fbgn_codes):
    """
    For a list of FBgn codes return a DF returning the Gene ID, Gene Symbol and Gene Name for drosophila
    :param fbgn_codes: List of strings : list of fbgn_codes:
    :return:DF with FBgn as index and columns with Gene, Symbol and name
    """
    logger.info("Producing the gene code df")
    m_dict = get_gene_code_dict(fbgn_codes)
    name_df = pd.DataFrame(m_dict).T
    name_df.columns = ["Annotation", "Symbol", "Name"]
    name_df.index.name = 'FBgn'
    # name_df.reset_index()
    return name_df


def get_gene_code_dict(fbgn_codes):
    """
    For a list of FBgn codes return a dict with a tuple of identifiers ( Gene ID, Gene Symbol and Gene Name)
    :param fbgn_codes: List of strings : list of fbgn_codes:
    :return: Dict of format {'FBgn0263002': ('CR43310', 'CR43310', gene_name)}
    """
    # TODO: write out each dict as JSON and store in a DB for easy access for each set of orthologues.

    logger.info("Getting the gene code dictionary")
    dros_genes_df = get_dros_genes_df()
    fb_genes = {}
    for fbgn in fbgn_codes:
        fb_id, fb_sy, fb_name = None, None, None
        if str(fbgn) != 'nan':
            try:
                row = dros_genes_df[dros_genes_df.index == fbgn]
                fb_id = row.Annotation[0]
                fb_sy = row.Symbol[0]
                fb_name = row.Name[0]
            except Exception as e:
                logger.warning("Returning None for {} as {} ".format(fbgn, e))
            fb_genes[fbgn] = fb_id, fb_sy, fb_name

    return fb_genes


def add_dros_names(species_df, identifier):
    """
    :param species_df: A data frame with FBgn codes that we want to add other gene identifiers to
    :param identifier: String, taken from species:id dict can be FBgb or other identifier.
    :return: species DF with new naming columns
    """
    logger.info("Adding dros Annotations, Names and Symbols")

    # Add on the gene names and symbols
    name_df = get_gene_code_df(list(species_df.index.get_level_values('FBgn')))
    name_df.index = name_df.index.set_names(['FBgn'])
    name_df = name_df.reset_index()

    display(name_df)

    # # Merge the two DFs together - not sure we want this here - these need to be in a DB structure
    species_df = species_df.reset_index()

    columns_to_drop = ['Symbol', 'Name', 'Annotation']

    existing_cols_to_drop = [col for col in columns_to_drop if col in species_df.columns]
    if existing_cols_to_drop:
        species_df = species_df.drop(columns=existing_cols_to_drop)

    species_df = species_df.merge(name_df, on='FBgn', how='left')
    species_df.drop_duplicates(inplace=True)

    # Put the index back to a multi-index or just FBgn if no identifier passed

    if identifier == 'FBgn':
        species_df = species_df.set_index([identifier])
    else:
        species_df = species_df.set_index([identifier, 'FBgn'])

    return species_df


def split_coding_non_coding_df(ab_en_df):
    """
    This splits the abundance/enrichemnt DF into non_coding and everything else (if not CR - it is kept)
    :param ab_en_df: The abundance/enrichment df with a 'Gene' column
    :return: two dataframes with non-coding genes and everything else (coding)
    """
    ab_en_df['Annotation'] = ab_en_df['Annotation'].fillna('-')

    coding = ab_en_df[ab_en_df.Annotation.str.startswith(("CG", "-"))]
    non_coding = ab_en_df[ab_en_df.Annotation.str.startswith("CR")]

    # create a new DataFrame where values in column A are greater than 2 and values in column B are less than 8

    logger.info("Returning {} coding and {} non-coding genes".format(coding.shape[0], non_coding.shape[0]))

    return coding, non_coding


def remove_ids(start_list, to_remove):
    """
    A method to simply remove elements from a list and return the list
    :param start_list: List with offending elements
    :param to_remove: List of elements you want to remove
    :return: Cleaned list
    """

    if to_remove is not None:
        return_list = [x for x in start_list if x not in to_remove]
    else:
        return_list = start_list

    return return_list


def remove_dict_items(my_dict, to_remove):
    """
    A method to remove elements from any dictionary keys list and return cleaned dictionary
    :param my_dict: Dictionary that the elements should be removed from
    :param to_remove: A list of elements that should be removed from the dictionary keys
    :return:
    """
    # Iterate through dict1
    for key, values in my_dict.items():
        # Check if values exist in dict2 and remove them
        my_dict[key] = [value for value in values if value not in to_remove]

    return my_dict


def write_csv(df, filename, tissue, RESULTS_FOLDER):
    logger.info("We are attempting to write out: %s" % filename)
    display(df.head())

    write_dir = os.path.join(RESULTS_FOLDER, tissue)
    write_csv_to_dir(df, write_dir, filename)

    # Fixme: This can probably be deleted - code moved below for usability
    #     if os.path.exists(write_dir):
    #         df.to_csv(os.path.join(write_dir, filename))
    #         logger.info("Wrote out %s" % os.path.join(write_dir, filename))
    #     else:
    #         os.makedirs(write_dir)
    #         df.to_csv(os.path.join(write_dir, filename))
    #         logger.info("Wrote out %s" % os.path.join(write_dir, filename))
    #
    # except Exception as e:
    #     logger.warning("Writing failed because of %s " % e)

def write_csv_to_dir(df, write_dir, filename):
    """
    :param df: DF to save as CSV
    :param write_dir: The full write directory
    :param filename: String The fileName of the CSV_file
    :return:
    """
    logger.info(f'Attempting to Write out {filename} to {write_dir}')
    try:
        if os.path.exists(write_dir):
            df.to_csv(os.path.join(write_dir, filename))
            logger.info("Wrote out %s" % os.path.join(write_dir, filename))
        else:
            os.makedirs(write_dir)
            df.to_csv(os.path.join(write_dir, filename))
            logger.info("Wrote out %s" % os.path.join(write_dir, filename))

    except Exception as e:
        logger.warning("Writing failed because of %s " % e)
        raise e

def write_json_to_dir(dct, write_dir, filename):
    """
    :param dct: dict to save as json
    :param write_dir: The full write directory
    :param filename: String The fileName of the CSV_file
    :return:
    """

    write_dir = Path(write_dir)
    full_write = write_dir / filename
    logger.info(f'Attempting to Write out {filename} to {write_dir}')

    if not write_dir.exists():
        logger.info(f"Making dir {write_dir}")
        write_dir.mkdir(parents=True)

    try:
        with open(full_write, 'w') as json_file:
            json.dump(dct, json_file)
        logger.info(f"Wrote out {full_write}")
    except Exception as e:
        logger.error(f"Error writing {full_write}: {e}")


def create_directory(directory_path):
    """
    Check if a directory exists, if not create it
    :param directory_path: The path to use
    :return:
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")


# Given a list of gene names get the FBgn code
def get_fbgn_from_name_old(identifiers):

    print (identifiers)
    dros_annot_df = get_dros_genes_df()
    fbgn_list = [dros_annot_df[dros_annot_df.Symbol == i].index.values[0] for i in identifiers]
    print ("returning", fbgn_list)
    return fbgn_list

def get_fbgn_from_name(identifiers):
    "Takes a list of identifiers and returns a list of fbgn codes"

    fbgn_list = []
    dros_annot_df = get_dros_genes_df()

    for i in identifiers:
        try:
            fbgn = dros_annot_df[dros_annot_df.Symbol == i].index.values[0]
            fbgn_list.append(fbgn)
        except IndexError as e:
            logger.warning(f'error: {e}; no fbgn code found for {i}, skipping')
    if not fbgn_list:
        fbgn_list = [None]

    return fbgn_list



def get_kendall_score(list1, list2):
    """
    1 indicates a perfect agreement in rankings (direct correlation).
    -1 indicates a perfect disagreement in rankings (inverse correlation).
    0 indicates no association between rankings.
    :param list1: List, to be compared
    :param list2: List, to be compared
    """
    # # Create a mapping from strings to numerical values
    mapping = {value: rank for rank, value in enumerate(sorted(set(list1 + list2)))}
    # Convert the lists to numerical ranks
    ranked_list1 = [mapping[item] for item in list1]
    ranked_list2 = [mapping[item] for item in list2]

    # Calculate Kendall Tau
    tau, _ = kendalltau(ranked_list1, ranked_list2)

    return tau

def check_same_order(l1, l2):
    """
    Check if common elements in two lists appear in the same order and return True or False accordingly
    :param l1: List
    :param l2: List
    :return: same_order: Boolean
    """

    # Find common elements in both lists while maintaining order
    common_elements_a = [x for x in l1 if x in l2]
    common_elements_b = [x for x in l2 if x in l1]

    same_order = common_elements_a == common_elements_b

    return same_order # All elements are in the same order across lists

def normalise_vector(lfc_column, scaler_type):
    """
    Method to take a colum from a df and return a normalised vector
    :param scaler_type: String, Type of normalisation to perform, standard, robust, minmax
    :param lfc_column: DF column with values to be normalised
    :return: normalised_vector, np.array, normalised values
    """
    vector_array = np.array(lfc_column)

    # Reshape the vector to a make sure itsa single column
    vector = vector_array.reshape(-1, 1)

    # Create MinMaxScaler instance
    if scaler_type == "standard":
        scaler = StandardScaler()
    elif scaler_type == "robust":
        scaler = RobustScaler()
    elif scaler_type == "minmax":

        # Define the scaler
        scaler = MinMaxScaler(feature_range=(-10, 10))

    else:
        logger.warning("please pass standard, robust or minmax for scaler type")
        raise

    # Fit and transform the vector using MinMaxScaler
    normalized_vector = scaler.fit_transform(vector)

    return normalized_vector

def remove_duplicates_in_order(original_list):
    """
    :param original_list: List of elements with (or without duplicates)
    :return: List of elements without duplicates in order
    """

    # Create a set to track seen elements
    seen = set()

    # Use list comprehension to iterate through the original list
    # Only add elements to the resulting list if they have not been seen before
    unique_list = [element for element in original_list if element not in seen and not seen.add(element)]

    return unique_list

def write_data_to_json(dict_or_list, output_filename, data_dir):
    """
    Write a Python object to a JSON file in the project data directory.

    Parameters
    ----------
    dict_or_list : any JSON-serialisable object e.g. list, dict
    output_filename : str - Output filename (e.g. 'RPE_high_markers.json')
    data_dir : str or Path - Path to the project data directory
    """

    data_dir = Path(data_dir)
    out_path = data_dir / output_filename

    with open(out_path, "w") as json_file:
        json.dump(dict_or_list, json_file)

    return out_path
