import os
import sys
import requests
import numpy as np
import pandas as pd
import scipy.stats as spstats
from tqdm.notebook import tqdm
from joblib import Parallel, delayed


def check_dir(dir: str):
    """
    Creates a given path driectory "dir" if it does not exist.
    Args:
        dir (str): Path to the directory.
    """
    if os.path.exists(dir) and os.path.isdir(dir):
        pass
    else:
        os.makedirs(dir)


def string_protein_search(string_ids: list, hgnc_table: pd.DataFrame):
    """
    Retrieve gene symbols for a list of STRING protein IDs.

    This function takes a list of STRING protein IDs and a pandas DataFrame containing HGNC gene data.
    It then uses the Ensembl REST API to look up the corresponding transcript and gene IDs for each protein ID.
    Finally, it merges the gene information with the HGNC data to return a DataFrame containing the original
    STRING protein IDs and their corresponding gene symbols.

    Parameters:
    string_ids (list): A list of STRING protein IDs.
    hgnc_table (pd.DataFrame): A pandas DataFrame containing HGNC gene data, with columns 'ensembl_gene_id' and 'symbol'.

    Returns:
    pd.DataFrame: A pandas DataFrame with columns 'string_ids' and 'symbol', containing the original STRING protein IDs and their corresponding gene symbols.
    """
    proteins = {prot: prot.split(".")[1] for prot in string_ids}
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # The first query returned transcript symbols
    r = requests.post(
        server + ext, headers=headers, json={"ids": list(proteins.values())}
    )
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    transcripts = {}
    for k, v in proteins.items():
        if decoded[v] is not None:
            transcripts[k] = decoded[v]["Parent"]

    # The second query will return gene symbols
    r = requests.post(
        server + ext, headers=headers, json={"ids": list(transcripts.values())}
    )
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    genes = []
    for k, v in transcripts.items():
        if decoded[v] is not None:
            genes.append((k, decoded[v]["Parent"]))

    genes = pd.DataFrame(genes, columns=pd.Index(["string_ids", "ensembl_gene_id"]))
    genes = pd.merge(genes, hgnc_table)[["string_ids", "symbol"]]

    return genes


def spearman_corr(expression, mutation_burden):
    rho, pval = spstats.spearmanr(
        expression,
        mutation_burden,
    )
    return pd.Series(
        {
            "rho": rho,
            "pval": pval,
        }
    )


def driver_neighbour_corr(
    driver: int,
    neighbourlist: np.ndarray,
    mutationtab: np.ndarray,
    expressiontab: np.ndarray,
    ctfilter: np.ndarray | None = None,
) -> tuple:
    if ctfilter is not None:
        filt = ctfilter[driver]
    else:
        filt = np.array([True] * len(mutationtab))

    exp = expressiontab[np.ix_(filt, neighbourlist)]
    mask = np.flatnonzero(exp.sum(axis=0) > 0)

    if len(mask) == 0:
        # if driver has no expression in any neighbour
        return [np.nan], [np.nan], [np.nan]

    elif len(mask) == 1:
        # if driver has only one neighbour with expression
        rho, pvalue = spstats.spearmanr(exp[:, mask], mutationtab[filt, driver])
        return [rho], [pvalue], neighbourlist[mask]

    else:
        rho, pvalue = spstats.spearmanr(exp[:, mask], mutationtab[filt, driver])
        return rho[-1, :-1], pvalue[-1, :-1], neighbourlist[mask]


def between_cancer_corr(
    exparray: np.ndarray,
    mutationarray: np.ndarray,
    graph: pd.DataFrame,
    cancertype_filter: np.ndarray | None = None,
    n_jobs: int = -1,
    progressbar: bool = True,
):
    driverlist = graph.columns.tolist()
    if progressbar:
        neighbourlist = tqdm([np.flatnonzero(graph[driver]) for driver in driverlist])
    else:
        neighbourlist = [np.flatnonzero(graph[driver]) for driver in driverlist]

    corr_results = Parallel(n_jobs=n_jobs)(
        delayed(driver_neighbour_corr)(
            driver, neighbours, mutationarray, exparray, cancertype_filter
        )
        for driver, neighbours in enumerate(neighbourlist)
    )

    results = {
        "driver": [],
        "neighbour": [],
        "rho": [],
        "rho_pval": [],
    }
    for driver in range(len(driverlist)):
        results["driver"].extend([driverlist[driver]] * len(corr_results[driver][2]))
        results["neighbour"].extend(corr_results[driver][2])
        results["rho"].extend(corr_results[driver][0])
        results["rho_pval"].extend(corr_results[driver][1])

    results = (
        pd.DataFrame(results)
        .set_index("neighbour")
        .merge(
            pd.Series(graph.index).rename("neighbour"),
            left_index=True,
            right_index=True,
        )
        .dropna()
        .reset_index(drop=True)
    )
    return results
