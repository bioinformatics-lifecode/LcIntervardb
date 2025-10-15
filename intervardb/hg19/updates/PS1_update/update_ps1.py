import gzip
import re

def parse_protein_change(name):
    """Extract protein change from Name field"""
    match = re.search(r'\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)', name)
    if match:
        ref_aa_3letter = match.group(1)
        position = match.group(2)
        alt_aa_3letter = match.group(3)
        
        # Convert 3-letter to 1-letter amino acid codes
        aa_map = {
            'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C',
            'Gln':'Q', 'Glu':'E', 'Gly':'G', 'His':'H', 'Ile':'I',
            'Leu':'L', 'Lys':'K', 'Met':'M', 'Phe':'F', 'Pro':'P',
            'Ser':'S', 'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V'
        }
        
        ref_aa = aa_map.get(ref_aa_3letter)
        alt_aa = aa_map.get(alt_aa_3letter)
        
        # Return None if either amino acid is not in map (filters Ter, X, etc.)
        if not ref_aa or not alt_aa:
            return None, None, None
            
        protein_notation = f"p.{ref_aa}{position}{alt_aa}"
        return ref_aa, alt_aa, protein_notation
    return None, None, None

# Process ClinVar file
with gzip.open('variant_summary.txt.gz', 'rt') as f:
    header = f.readline().strip().split('\t')
    
    with open('PS1.AA.change.patho.hg19.updated', 'w') as out:
        for line in f:
            fields = line.strip().split('\t')
            
            # Filter criteria
            variant_type = fields[1]
            clinical_sig = fields[6]
            assembly = fields[16]
            chromosome = fields[18]
            start = fields[19]
            ref_vcf = fields[32]
            alt_vcf = fields[33]
            name = fields[2]
            rs_id = fields[9]
            var_id = fields[30]
            
            # Apply filters
            if assembly != 'GRCh37':
                continue
            if variant_type != 'single nucleotide variant':
                continue
            if 'Pathogenic' not in clinical_sig:
                continue
            if ref_vcf == 'na' or alt_vcf == 'na':
                continue
                
            # Parse protein change - will return None for stop codons
            ref_aa, alt_aa, protein_notation = parse_protein_change(name)
            if not ref_aa or not alt_aa:
                continue
            
            # Use RS# if available, otherwise VariationID
            variant_id = rs_id if rs_id != '-1' else var_id
            
            # Write output
            out.write(f"{chromosome}\t{start}\t{start}\t{ref_vcf}\t{alt_vcf}\t{ref_aa}\t{alt_aa}\t{protein_notation}\t{variant_id}\n")

print("PS1 file created: PS1.AA.change.patho.hg19.updated")
