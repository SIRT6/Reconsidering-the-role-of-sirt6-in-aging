#!/usr/bin/env python3
"""
Map academic papers to aging hallmarks based on PyAlex concept scores.
"""

import json
from collections import defaultdict

# Define aging hallmarks
HALLMARKS = [
    "genome instability",
    "telomere attrition",
    "epigenetic alterations",
    "loss of proteostasis",
    "disabled macroautophagy",
    "deregulated nutrient sensing",
    "mitochondrial dysfunction",
    "cellular senescence",
    "stem cell exhaustion",
    "intercellular communication",
    "chronic inflammation",
    "dysbiosis"
]

# Define concept mappings to hallmarks
# Maps concept display names (lowercase) to hallmark categories
CONCEPT_TO_HALLMARK = {
    # Genome instability
    "genome instability": ["genome instability"],
    "genomic instability": ["genome instability"],
    "dna damage": ["genome instability"],
    "dna repair": ["genome instability"],
    "double-strand break": ["genome instability"],
    "mutation": ["genome instability"],
    "chromosomal instability": ["genome instability"],
    "instability": ["genome instability"],
    "dna": ["genome instability"],
    
    # Telomere attrition
    "telomere": ["telomere attrition"],
    "telomerase": ["telomere attrition"],
    "telomere shortening": ["telomere attrition"],
    "telomere length": ["telomere attrition"],
    
    # Epigenetic alterations
    "epigenetics": ["epigenetic alterations"],
    "histone": ["epigenetic alterations"],
    "histone deacetylase": ["epigenetic alterations"],
    "acetylation": ["epigenetic alterations"],
    "methylation": ["epigenetic alterations"],
    "chromatin": ["epigenetic alterations"],
    "dna methylation": ["epigenetic alterations"],
    "histone modification": ["epigenetic alterations"],
    "sirtuin": ["epigenetic alterations"],
    "sirt1": ["epigenetic alterations"],
    "sirt2": ["epigenetic alterations"],
    "sirt3": ["epigenetic alterations"],
    "sirt6": ["epigenetic alterations"],
    "chromatin remodeling": ["epigenetic alterations"],
    
    # Loss of proteostasis
    "proteostasis": ["loss of proteostasis"],
    "proteasome": ["loss of proteostasis"],
    "ubiquitin": ["loss of proteostasis"],
    "protein folding": ["loss of proteostasis"],
    "unfolded protein response": ["loss of proteostasis"],
    "heat shock protein": ["loss of proteostasis"],
    "chaperone": ["loss of proteostasis"],
    "protein aggregation": ["loss of proteostasis"],
    "protein degradation": ["loss of proteostasis"],
    "endoplasmic reticulum stress": ["loss of proteostasis"],
    
    # Disabled macroautophagy
    "autophagy": ["disabled macroautophagy"],
    "macroautophagy": ["disabled macroautophagy"],
    "mitophagy": ["disabled macroautophagy"],
    "lysosome": ["disabled macroautophagy"],
    "beclin": ["disabled macroautophagy"],
    "atg": ["disabled macroautophagy"],
    "autophagosome": ["disabled macroautophagy"],
    
    # Deregulated nutrient sensing
    "mtor": ["deregulated nutrient sensing"],
    "ampk": ["deregulated nutrient sensing"],
    "insulin signaling": ["deregulated nutrient sensing"],
    "igf-1": ["deregulated nutrient sensing"],
    "nutrient sensing": ["deregulated nutrient sensing"],
    "caloric restriction": ["deregulated nutrient sensing"],
    "metabolism": ["deregulated nutrient sensing"],
    "glucose metabolism": ["deregulated nutrient sensing"],
    "insulin": ["deregulated nutrient sensing"],
    
    # Mitochondrial dysfunction
    "mitochondria": ["mitochondrial dysfunction"],
    "mitochondrial": ["mitochondrial dysfunction"],
    "oxidative phosphorylation": ["mitochondrial dysfunction"],
    "reactive oxygen species": ["mitochondrial dysfunction"],
    "oxidative stress": ["mitochondrial dysfunction"],
    "mitochondrial dna": ["mitochondrial dysfunction"],
    "electron transport chain": ["mitochondrial dysfunction"],
    "atp": ["mitochondrial dysfunction"],
    
    # Cellular senescence
    "senescence": ["cellular senescence"],
    "cellular senescence": ["cellular senescence"],
    "senescent": ["cellular senescence"],
    "sasp": ["cellular senescence"],
    "p16": ["cellular senescence"],
    "p21": ["cellular senescence"],
    "cell cycle arrest": ["cellular senescence"],
    
    # Stem cell exhaustion
    "stem cell": ["stem cell exhaustion"],
    "stem cells": ["stem cell exhaustion"],
    "progenitor cell": ["stem cell exhaustion"],
    "hematopoietic stem cell": ["stem cell exhaustion"],
    "regeneration": ["stem cell exhaustion"],
    "tissue regeneration": ["stem cell exhaustion"],
    
    # Intercellular communication
    "cell signaling": ["intercellular communication"],
    "signaling pathway": ["intercellular communication"],
    "cytokine": ["intercellular communication"],
    "extracellular matrix": ["intercellular communication"],
    "cell adhesion": ["intercellular communication"],
    
    # Chronic inflammation
    "inflammation": ["chronic inflammation"],
    "inflammatory": ["chronic inflammation"],
    "inflammaging": ["chronic inflammation"],
    "nf-κb": ["chronic inflammation"],
    "nfkb": ["chronic inflammation"],
    "interleukin": ["chronic inflammation"],
    "tnf": ["chronic inflammation"],
    "tumor necrosis factor": ["chronic inflammation"],
    "immune response": ["chronic inflammation"],
    
    # Dysbiosis
    "microbiome": ["dysbiosis"],
    "microbiota": ["dysbiosis"],
    "gut microbiome": ["dysbiosis"],
    "dysbiosis": ["dysbiosis"],
    "bacteria": ["dysbiosis"],
    "bacterial": ["dysbiosis"],
}


def load_papers(filepath):
    """Load papers with concepts from JSON file."""
    with open(filepath, 'r') as f:
        return json.load(f)


def map_paper_to_hallmarks(paper_title, concepts):
    """
    Map a paper to aging hallmarks based on its concepts.
    
    Returns:
        dict: Dictionary with 'hallmarks' (dict of hallmark->score) and 
              'not-related to aging' (bool)
    """
    hallmark_scores = defaultdict(float)
    
    for concept in concepts:
        concept_name = concept['display_name'].lower()
        concept_score = concept['score']
        
        # Check if this concept maps to any hallmarks
        if concept_name in CONCEPT_TO_HALLMARK:
            for hallmark in CONCEPT_TO_HALLMARK[concept_name]:
                hallmark_scores[hallmark] += concept_score
    
    # Prepare result
    if hallmark_scores:
        return {
            "hallmarks": dict(hallmark_scores),
            "not-related to aging": False
        }
    else:
        return {
            "hallmarks": {},
            "not-related to aging": True
        }


def process_all_papers(papers_data):
    """
    Process all papers and map them to hallmarks.
    
    Args:
        papers_data: Dictionary mapping paper titles to concept lists
        
    Returns:
        dict: Processed results for all papers
    """
    results = {}
    
    for paper_title, concepts in papers_data.items():
        results[paper_title] = map_paper_to_hallmarks(paper_title, concepts)
    
    return results


def main():
    """Main execution function."""
    # Load input data
    input_file = '/mnt/user-data/uploads/papers_with_concepts.json'
    papers_data = load_papers(input_file)
    
    print(f"Loaded {len(papers_data)} papers")
    
    # Process papers
    results = process_all_papers(papers_data)
    
    # Save results
    output_file = '/mnt/user-data/outputs/papers_mapped_to_hallmarks.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    
    # Print summary statistics
    total_papers = len(results)
    aging_related = sum(1 for r in results.values() if not r['not-related to aging'])
    not_aging_related = total_papers - aging_related
    
    print(f"\n=== Summary Statistics ===")
    print(f"Total papers: {total_papers}")
    print(f"Aging-related papers: {aging_related} ({aging_related/total_papers*100:.1f}%)")
    print(f"Not related to aging: {not_aging_related} ({not_aging_related/total_papers*100:.1f}%)")
    
    # Count papers per hallmark
    hallmark_counts = defaultdict(int)
    for result in results.values():
        for hallmark in result['hallmarks'].keys():
            hallmark_counts[hallmark] += 1
    
    print(f"\n=== Papers per Hallmark ===")
    for hallmark in sorted(hallmark_counts.keys(), key=lambda x: hallmark_counts[x], reverse=True):
        print(f"{hallmark}: {hallmark_counts[hallmark]} papers")
    
    # Show a few examples
    print(f"\n=== Sample Results ===")
    for i, (paper_title, result) in enumerate(list(results.items())[:3]):
        print(f"\n{i+1}. {paper_title}")
        if result['not-related to aging']:
            print("   → Not related to aging")
        else:
            print(f"   → Hallmarks:")
            for hallmark, score in sorted(result['hallmarks'].items(), key=lambda x: x[1], reverse=True):
                print(f"      • {hallmark}: {score:.3f}")


if __name__ == "__main__":
    main()
