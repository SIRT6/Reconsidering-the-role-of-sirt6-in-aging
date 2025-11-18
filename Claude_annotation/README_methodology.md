# Aging Hallmarks Paper Mapping - Methodology and Results

## Overview
This analysis maps 2,147 academic papers to aging hallmarks based on PyAlex concept scores. Each paper is assigned to one or more of 12 predefined aging hallmarks based on the presence and relevance of associated concepts.

## Aging Hallmarks
The analysis uses the following 12 hallmarks of aging:
1. **Genome instability**
2. **Telomere attrition**
3. **Epigenetic alterations**
4. **Loss of proteostasis**
5. **Disabled macroautophagy**
6. **Deregulated nutrient sensing**
7. **Mitochondrial dysfunction**
8. **Cellular senescence**
9. **Stem cell exhaustion**
10. **Intercellular communication**
11. **Chronic inflammation**
12. **Dysbiosis**

## Methodology

### Concept-to-Hallmark Mapping
Papers are mapped to hallmarks based on their PyAlex concepts. Each concept is associated with one or more hallmarks through a comprehensive mapping dictionary. Key mappings include:

#### Genome Instability
- DNA damage, DNA repair, double-strand breaks, mutations, chromosomal instability

#### Telomere Attrition
- Telomere, telomerase, telomere shortening, telomere length

#### Epigenetic Alterations
- Histone modifications, acetylation, methylation, chromatin remodeling
- Sirtuins (SIRT1-6), histone deacetylases

#### Loss of Proteostasis
- Proteasome, ubiquitin, protein folding, unfolded protein response
- Heat shock proteins, chaperones, protein aggregation

#### Disabled Macroautophagy
- Autophagy, macroautophagy, mitophagy, lysosomes, autophagosomes

#### Deregulated Nutrient Sensing
- mTOR, AMPK, insulin signaling, IGF-1, caloric restriction, metabolism

#### Mitochondrial Dysfunction
- Mitochondria, oxidative phosphorylation, reactive oxygen species
- Oxidative stress, electron transport chain

#### Cellular Senescence
- Senescence, SASP, p16, p21, cell cycle arrest

#### Stem Cell Exhaustion
- Stem cells, progenitor cells, tissue regeneration

#### Intercellular Communication
- Cell signaling, cytokines, extracellular matrix, cell adhesion

#### Chronic Inflammation
- Inflammation, inflammaging, NF-κB, interleukins, TNF, immune response

#### Dysbiosis
- Microbiome, microbiota, gut microbiome, bacteria

### Scoring System
For each paper:
1. Extract all PyAlex concepts with their relevance scores
2. Match concept names (case-insensitive) to the hallmark mapping dictionary
3. Sum scores for all concepts mapped to each hallmark
4. Papers with no matching concepts are classified as "not-related to aging"

## Results Summary

### Overall Statistics
- **Total papers analyzed**: 2,147
- **Aging-related papers**: 1,568 (73.0%)
- **Not related to aging**: 579 (27.0%)
- **Papers with multiple hallmarks**: 729 (46.5% of aging-related papers)

### Papers per Hallmark
| Hallmark | Number of Papers | Percentage |
|----------|-----------------|------------|
| Epigenetic alterations | 1,203 | 56.0% |
| Genome instability | 267 | 12.4% |
| Chronic inflammation | 267 | 12.4% |
| Mitochondrial dysfunction | 204 | 9.5% |
| Deregulated nutrient sensing | 173 | 8.1% |
| Cellular senescence | 135 | 6.3% |
| Disabled macroautophagy | 100 | 4.7% |
| Stem cell exhaustion | 81 | 3.8% |
| Telomere attrition | 55 | 2.6% |
| Loss of proteostasis | 39 | 1.8% |
| Intercellular communication | 31 | 1.4% |
| Dysbiosis | 7 | 0.3% |

### Key Findings

1. **Epigenetic alterations** is the most represented hallmark, appearing in over half of all papers. This reflects the strong focus on sirtuins and histone modifications in the dataset.

2. **Multi-hallmark papers** represent nearly half of aging-related papers, indicating that many studies explore interconnections between different aging mechanisms.

3. **Comprehensive aging studies**: Some papers map to 6-7 hallmarks simultaneously, particularly review papers covering multiple aspects of aging biology.

### Example High-Impact Papers

#### Paper with Most Hallmarks (7)
**"The Role and Molecular Pathways of SIRT6 in Senescence and Age‐related Diseases"**
- Epigenetic alterations: 2.219
- Genome instability: 1.251
- Telomere attrition: 1.142
- Disabled macroautophagy: 0.785
- Cellular senescence: 0.563
- Mitochondrial dysfunction: 0.544
- Deregulated nutrient sensing: 0.443

#### Strong Single-Hallmark Examples
**"SIRT6 Promotes DNA Repair Under Stress by Activating PARP1"**
- Genome instability: 2.328 (focused study on DNA repair)

**"SIRT6 is a histone H3 lysine 9 deacetylase that modulates telomeric chromatin"**
- Epigenetic alterations: 3.860
- Telomere attrition: 0.654
- Genome instability: 0.112

## Output Format

The output JSON file (`papers_mapped_to_hallmarks.json`) has the following structure:

```json
{
  "Paper Title": {
    "hallmarks": {
      "hallmark_name1": score1,
      "hallmark_name2": score2
    },
    "not-related to aging": false
  }
}
```

For papers not related to aging:
```json
{
  "Paper Title": {
    "hallmarks": {},
    "not-related to aging": true
  }
}
```

## Limitations and Considerations

1. **Concept coverage**: The mapping is based on predefined concept-to-hallmark associations. Papers discussing aging through concepts not in our mapping dictionary may be misclassified.

2. **Score interpretation**: Scores represent cumulative relevance of concepts to each hallmark, not necessarily the paper's primary focus.

3. **Multi-disciplinary papers**: Papers covering multiple topics may have lower scores per hallmark despite being relevant to multiple aging mechanisms.

4. **SIRT-heavy dataset**: The dataset appears to be enriched for sirtuin-related papers, leading to high representation of epigenetic alterations.

## Usage for Downstream Analysis

The output can be used for:
- **Identifying research gaps**: Hallmarks with fewer papers may represent underexplored areas
- **Finding cross-cutting studies**: Papers with multiple hallmarks reveal mechanistic connections
- **Focused literature review**: Filter papers by specific hallmarks of interest
- **Trend analysis**: Examine how research focus has shifted across hallmarks over time
- **Collaboration discovery**: Identify research groups working on specific hallmark combinations

## Files Included

1. `papers_mapped_to_hallmarks.json` - Main output file with all paper mappings
2. `map_papers_to_hallmarks.py` - Python script used for the analysis
3. `README_methodology.md` - This documentation file

## Contact and Citation

This analysis was performed using PyAlex concept data and custom hallmark mapping logic. For questions about the methodology or to report issues with specific paper classifications, please review the concept-to-hallmark mapping dictionary in the source code.
