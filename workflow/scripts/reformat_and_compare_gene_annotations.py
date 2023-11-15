#!/usr/bin/env python

'''
Reformat annotated BED file and compare gene annotations to provided genes from FASTA descriptions
Lourdes Rosano, Nov 2023
'''

import sys
import os
import argparse

'''
function definitions
'''

# Parse input file listing one region per row (there can be duplicates) with annotations
# Assumes no header and 4 columns: chrom, start, end, annotation
def read_input(regions_new_mapping, regions_original_mapping):
    # Keep track of original sorting in input file
    sorted_regions = []  
    infile = open(args.infile, "r")
    for line in infile:
        line_split = line.strip().split("\t")
        # Sanity check for having exactly 9 columns
        if len(line_split) != 9:
            raise ValueError("Input BED file does not contain 9 columns. Please check that you provide (in order): chrom, start, end, original_annotation, intersected_chrom, intersected_start, intersected_end, intersected_strand, intersected_annotation.")

        chrom = line_split[0].strip()
        start = int(line_split[1].strip())
        end = int(line_split[2].strip())
        # Use uppercase for consistency of the gene annotations
        annot = line_split[3].strip().upper()
        new_annot = line_split[8].strip().upper()

        # Sanity check for special case of original gene annotations with format '[GENE]_PosCtrl'
        original_annot = annot.split("_")[0].strip()

        # Create unique region ID
#         region_name = str(chrom) + "__" + str(start) + "__" + str(end) + "__" + original_annot
        region_name = str(chrom) + "__" + str(start) + "__" + str(end) + "__" + annot

        # Sanity check that at least 1 gene annotation was found for the current region
        if not new_annot:
            raise ValueError("Region %s contains no annotations!" %(region_name))

        # Keep track of the order of regions in the BED file
        if region_name not in sorted_regions:
            sorted_regions.append(region_name)

        # Keep track of each unique region name and all the new (intersected) gene annotations associated to them
        if region_name not in regions_new_mapping.keys():
            regions_new_mapping[region_name] = []
        if new_annot not in regions_new_mapping[region_name]:
            regions_new_mapping[region_name].append(new_annot)

        # Keep track of the originally assigned genes for each unique region name
        if region_name not in regions_original_mapping.keys():
            regions_original_mapping[region_name] = []
        if original_annot not in regions_original_mapping[region_name]:
            regions_original_mapping[region_name].append(original_annot)

    infile.close()
    return (sorted_regions, regions_new_mapping, regions_original_mapping)


'''
main script
'''

parser = argparse.ArgumentParser(description='Reformat annotated BED file (from bedtools intersect): regions overlapping with more than one annotation will be collapsed into one single row. Also, compare intersected gene annotations to provided gene names from FASTA descriptions: output file of regions where no match could be found between gene annotations (see more details in argument description).')
parser.add_argument('--infile', dest='infile', required=True, help='Input BED file possibly containing duplicated regions due to having >1 overlapping gene annotations per region. Assumes no header and 9 columns: chrom, start, end, original_annotation, intersected_chrom, intersected_start, intersected_end, intersected_strand, intersected_annotation.')
parser.add_argument('--outfile', dest='outfile', required=True, help='Output BED file containing 1 row per unique region with grouped annotations (separated with ";").')

args = parser.parse_args()

# Parse input BED file
regions_new_mapping = {}
regions_original_mapping = {}
(sorted_regions, regions_new_mapping, regions_original_mapping) = read_input(regions_new_mapping, regions_original_mapping)

# Write reformatted file content to output file
# At the same time:
#  - Compare new (intersected) gene annotations assigned to each region to the original gene annotations of the BED file
#  - Keep track of all regions where original and new gene annotations could not be matched (if any)
outfile = open(args.outfile, "w")

# Keep track of the number of regions where new and original gene annotations could be matched
n_matched = 0
# Keep track of the number of regions where new and original gene annotations could NOT be matched
n_unmatched = 0
# Keep track of the number of regions where new gene annotations were not available
n_empty = 0

# Use same order of regions as in the input BED file
for region in sorted_regions:
    new_annotations = regions_new_mapping[region]
    original_annotations = regions_original_mapping[region]

    # Keep track of all the matched genes found for the current region
    region_found_genes = []
    # Keep track of all the unmatched genes found for the current region
    region_not_found_genes = []

    final_annotation = None
    match_type = None

    # Iterate individual gene annotations for the current region
    for gene in new_annotations:
        # Skip cases where no gene annotation could be intersected with the region
        if (gene == "*"):
            continue

        # Check whether current gene annotation was also provided as the original gene annotation
        if gene in original_annotations:
            # Keep track of unique gene annotations matched for the current regions
            if gene not in region_found_genes:
                region_found_genes.append(gene)
        else:
            # Keep track of unique gene annotations not matched for the current regions
            if gene not in region_not_found_genes:
                region_not_found_genes.append(gene)


    # If current region had at least 1 matched gene annotation, keep it as the only annotation
    if region_found_genes:
        if len(region_found_genes) > 1:
            raise ValueError("Several gene annotations could be matched for region %s!" %(region))
        final_annotation = region_found_genes[0]
        match_type = "matched"
        n_matched += 1

    # If current region did not have at least 1 matched gene annotation, keep track of the unmatched region
    # Assign all new available gene annotations found for the region
    elif (not region_found_genes) and region_not_found_genes:
        final_annotation = ";".join(new_annotations)
        match_type = "not_matched"
        n_unmatched += 1

    # If no gene annotations were parsed for the current region, then none were overlapped and region is empty
    elif (not region_found_genes) and (not region_not_found_genes):
        final_annotation = "*"
        match_type = "empty"
        n_empty += 1

    else:
        raise ValueError("Unexpected condition was met for region %s!" %(region))

    # Retrieve again chrom, start, end and original gene from region name
    region_split = region.split("__")
    chrom = region_split[0]
    start = region_split[1]
    end = region_split[2]
    original_gene = region_split[3]

    outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(chrom, start, end, original_gene, final_annotation, match_type))

outfile.close()

print("Found %s regions where new and original gene annotations could be matched." %(str(n_matched)))
print("Found %s regions where new and original gene annotations could not be matched." %(str(n_unmatched)))
print("Found %s regions where new gene annotations were not available." %(str(n_empty)))
