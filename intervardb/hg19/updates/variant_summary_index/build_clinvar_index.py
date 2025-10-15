#!/usr/bin/env python
"""
Build ClinVar Gene Pathogenic LOF Variant Index
Creates a pre-computed index for fast InterVar PVS1 lookups
"""

import gzip
import re
import pickle
import sys
import os

def build_clinvar_gene_index(variant_summary_file, output_index_file):
    """
    Build complete gene index from variant_summary.txt.gz
    Saves to pickle file for fast loading
    """
    
    if not os.path.isfile(variant_summary_file):
        print("ERROR: variant_summary.txt.gz not found at: %s" % variant_summary_file)
        sys.exit(1)
    
    try:
        file_size_mb = os.path.getsize(variant_summary_file) / (1024 * 1024)
        print("=" * 70)
        print("ClinVar Pathogenic LOF Variant Index Builder")
        print("=" * 70)
        print("Input file: %s (%.0f MB)" % (variant_summary_file, file_size_mb))
        print("Output index: %s" % output_index_file)
        print("=" * 70)
    except OSError as e:
        print("ERROR: Cannot access file: %s" % str(e))
        sys.exit(1)
    
    gene_data = {}
    line_count = 0
    genes_found = set()
    variants_indexed = 0
    
    try:
        print("\nPhase 1: Scanning ClinVar variant_summary.txt.gz...")
        print("-" * 70)
        
        with gzip.open(variant_summary_file, 'rt', encoding='utf-8', errors='ignore') as f:
            # Skip header
            header = f.readline()
            
            for line in f:
                line_count += 1
                
                # Progress indicator every 500k lines
                if line_count % 500000 == 0:
                    print("Progress: %d lines | %d genes | %d variants" % (line_count, len(genes_found), variants_indexed))
                
                try:
                    cols = line.strip().split('\t')
                    if len(cols) < 35:
                        continue
                    
                    # Column indices: 4=GeneSymbol, 6=ClinicalSignificance, 2=Name
                    gene_symbol = cols[4]
                    clinsig = cols[6].lower()
                    name = cols[2]
                    
                    # Filter: Must be pathogenic (not conflicting)
                    if 'pathogenic' not in clinsig or 'conflicting' in clinsig:
                        continue
                    
                    # Check if it's LOF type (nonsense, frameshift, splice, etc.)
                    is_lof = False
                    lof_keywords = ['nonsense', 'frameshift', 'splice', 'stop', 'ter', 'fs', 
                                    'deletion', 'truncat', 'del', 'dup']
                    name_lower = name.lower()
                    for keyword in lof_keywords:
                        if keyword in name_lower:
                            is_lof = True
                            break
                    
                    if not is_lof:
                        continue
                    
                    # Initialize gene data if first time
                    if gene_symbol not in gene_data:
                        gene_data[gene_symbol] = {
                            'exon_counts': {},
                            'position_list': [],
                            'variant_names': []
                        }
                        genes_found.add(gene_symbol)
                    
                    # Extract exon number from name
                    exon_match = re.search(r'exon[\s:]*(\d+)', name, re.IGNORECASE)
                    if exon_match:
                        exon_num = int(exon_match.group(1))
                        gene_data[gene_symbol]['exon_counts'][exon_num] = gene_data[gene_symbol]['exon_counts'].get(exon_num, 0) + 1
                    
                    # Extract protein position from p. notation
                    # Formats: p.Arg714Ter, p.Gln918fs, p.R714*, p.Q918Tfs*24
                    prot_match = re.search(r'p\.[A-Z][a-z]{2}(\d+)', name)
                    if not prot_match:
                        prot_match = re.search(r'p\.[A-Z](\d+)', name)
                    
                    if prot_match:
                        prot_pos = int(prot_match.group(1))
                        gene_data[gene_symbol]['position_list'].append(prot_pos)
                    
                    # Store variant name for debugging
                    gene_data[gene_symbol]['variant_names'].append(name)
                    variants_indexed += 1
                
                except (ValueError, IndexError) as e:
                    continue
        
        print("-" * 70)
        print("Phase 1 Complete: Scanned %d lines" % line_count)
        print("  - Found %d genes with pathogenic LOF variants" % len(genes_found))
        print("  - Indexed %d pathogenic LOF variants" % variants_indexed)
        
    except (IOError, OSError) as e:
        print("ERROR: Failed to read variant_summary file: %s" % str(e))
        sys.exit(1)
    
    # Phase 2: Calculate truncated region counts
    print("\nPhase 2: Calculating truncated region counts...")
    print("-" * 70)
    
    for gene_symbol in gene_data:
        truncated_counts = {}
        position_list = gene_data[gene_symbol]['position_list']
        
        for pos in sorted(set(position_list)):
            count = sum(1 for p in position_list if p >= pos)
            truncated_counts[pos] = count
        
        gene_data[gene_symbol]['truncated_counts'] = truncated_counts
    
    print("Phase 2 Complete: Calculated truncated regions for %d genes" % len(gene_data))
    
    # Phase 3: Save to pickle file
    print("\nPhase 3: Saving index to disk...")
    print("-" * 70)
    
    try:
        with open(output_index_file, 'wb') as f:
            pickle.dump(gene_data, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        index_size_mb = os.path.getsize(output_index_file) / (1024 * 1024)
        print("Index saved successfully: %s (%.2f MB)" % (output_index_file, index_size_mb))
    
    except (IOError, OSError) as e:
        print("ERROR: Failed to save index file: %s" % str(e))
        sys.exit(1)
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("INDEX BUILD COMPLETE")
    print("=" * 70)
    print("Total genes indexed: %d" % len(gene_data))
    print("Total pathogenic LOF variants: %d" % variants_indexed)
    
    # Show top 10 genes by variant count
    gene_counts = [(gene, len(data['position_list'])) for gene, data in gene_data.items()]
    gene_counts.sort(key=lambda x: x[1], reverse=True)
    
    print("\nTop 10 genes by pathogenic LOF variant count:")
    for i, (gene, count) in enumerate(gene_counts[:10]):
        exon_count = sum(gene_data[gene]['exon_counts'].values())
        print("  %2d. %-10s : %d variants (%d with exon annotation)" % (i+1, gene, count, exon_count))
    
    print("\n" + "=" * 70)
    print("Index ready for use with InterVar!")
    print("=" * 70)
    
    return gene_data

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python build_clinvar_index.py <variant_summary.txt.gz> [output_index.pkl]")
        print("\nExample:")
        print("  python build_clinvar_index.py variant_summary.txt.gz clinvar_lof_index.pkl")
        sys.exit(1)
    
    variant_summary_file = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_index_file = sys.argv[2]
    else:
        # Default output name
        output_index_file = os.path.join(os.path.dirname(variant_summary_file), 'clinvar_lof_index.pkl')
    
    build_clinvar_gene_index(variant_summary_file, output_index_file)
