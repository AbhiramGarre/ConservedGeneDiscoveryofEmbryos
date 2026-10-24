#!/usr/bin/env python3
"""Cross-species candidate screening from precomputed within-species gene variability.

This exploratory script does not prove developmental expression conservation:
variation in unpaired populations can be driven by cell composition or stage.
Requires a verified Ensembl one-to-one ortholog mapping table.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def screen_orthologs(mouse_csv, zebrafish_csv, ortholog_file, output_dir, fraction=0.10):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    mouse = pd.read_csv(mouse_csv)
    zebra = pd.read_csv(zebrafish_csv)
    required_variability = {'ensembl_gene_id', 'log1p_variance'}
    if not required_variability.issubset(mouse.columns) or not required_variability.issubset(zebra.columns):
        raise ValueError('Variability CSVs must contain ensembl_gene_id and log1p_variance')
    mouse = mouse[['ensembl_gene_id', 'log1p_variance']].rename(
        columns={'ensembl_gene_id':'mouse_gene_id','log1p_variance':'mouse_log1p_variance'})
    zebra = zebra[['ensembl_gene_id', 'log1p_variance']].rename(
        columns={'ensembl_gene_id':'zebrafish_gene_id','log1p_variance':'zebrafish_log1p_variance'})

    ortho = pd.read_csv(ortholog_file, sep=None, engine='python', dtype=str)
    required = {'mouse_gene_id', 'zebrafish_gene_id'}
    if not required.issubset(ortho.columns):
        raise ValueError('Ortholog file must have mouse_gene_id and zebrafish_gene_id headers')
    if 'orthology_type' in ortho.columns:
        ortho = ortho[ortho['orthology_type'].str.lower().eq('ortholog_one2one')]
    else:
        raise ValueError('Require orthology_type so one-to-one assignments can be verified')
    ortho = ortho.dropna(subset=['mouse_gene_id', 'zebrafish_gene_id'])
    ortho = ortho[ortho.mouse_gene_id.str.startswith('ENSMUSG') &
                  ortho.zebrafish_gene_id.str.startswith('ENSDARG')]
    ortho = ortho.drop_duplicates(subset=['mouse_gene_id', 'zebrafish_gene_id'])
    # One-to-one guarantee: drop any problematic IDs that map to >1 partner.
    ortho = ortho[~ortho.mouse_gene_id.duplicated(keep=False)]
    ortho = ortho[~ortho.zebrafish_gene_id.duplicated(keep=False)]

    merged = ortho.merge(mouse, on='mouse_gene_id', how='inner', validate='one_to_one').merge(
        zebra, on='zebrafish_gene_id', how='inner', validate='one_to_one')
    n = len(merged)
    if n < 100:
        raise ValueError(f'Only {n} one-to-one mapped genes were found in both matrices. '
                         'Check Ensembl genome-build compatibility and ortholog IDs.')
    if not 0 < fraction < 1:
        raise ValueError('--top-fraction must be between 0 and 1')

    # Percentiles computed only among genes that have a verified, observable ortholog.
    # Each species is ranked separately, avoiding direct comparison of expression units.
    for sp in ['mouse', 'zebrafish']:
        merged[f'{sp}_variance_percentile'] = merged[f'{sp}_log1p_variance'].rank(
            method='average', pct=True)
        merged[f'{sp}_high_variability'] = merged[f'{sp}_variance_percentile'] >= 1 - fraction
    merged['high_variability_in_both'] = (merged.mouse_high_variability &
                                          merged.zebrafish_high_variability)
    merged['min_variance_percentile'] = merged[['mouse_variance_percentile',
                                                'zebrafish_variance_percentile']].min(axis=1)
    merged['mean_variance_percentile'] = merged[['mouse_variance_percentile',
                                                 'zebrafish_variance_percentile']].mean(axis=1)
    merged = merged.sort_values(['high_variability_in_both', 'min_variance_percentile',
                                 'mean_variance_percentile'], ascending=False)
    merged.to_csv(out/'ranked_one_to_one_ortholog_candidates.csv', index=False)

    k = int(merged['high_variability_in_both'].sum())
    m = int(merged['mouse_high_variability'].sum())
    z = int(merged['zebrafish_high_variability'].sum())
    expected = m*z/n
    p = float(hypergeom.sf(k-1, n, m, z))
    results = {
      'jointly_observed_one_to_one_orthologs': n,
      'top_fraction_per_species': fraction,
      'mouse_high_variability_genes': m,
      'zebrafish_high_variability_genes': z,
      'overlap_high_variability_genes': k,
      'expected_overlap_under_random_independent_selection': expected,
      'hypergeometric_upper_tail_pvalue': p,
      'warning': ('Exploratory screen only: shared variability does not establish '
                  'conserved developmental functions. Cell composition, processing, '
                  'and developmental-stage differences require follow-up validation.')
    }
    (out/'overlap_summary.json').write_text(json.dumps(results, indent=2))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mouse', required=True, help='mouse_gene_variability.csv')
    parser.add_argument('--zebrafish', required=True, help='zebrafish_gene_variability.csv')
    parser.add_argument('--orthologs', required=True, help='Verified Ensembl ortholog TSV/CSV')
    parser.add_argument('--output', default='conserved_gene_results')
    parser.add_argument('--top-fraction', type=float, default=0.10)
    args = parser.parse_args()
    print(json.dumps(screen_orthologs(args.mouse, args.zebrafish,
                                     args.orthologs, args.output, args.top_fraction), indent=2))
