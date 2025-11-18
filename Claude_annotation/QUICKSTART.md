# Quick Start Guide - Aging Hallmarks Mapping Results

## Files Included

### 1. Main Output Files
- **`papers_mapped_to_hallmarks.json`** (506 KB)
  - Primary output file with complete hallmark mapping
  - JSON format for programmatic access
  - Structure: `{"Paper Title": {"hallmarks": {...}, "not-related to aging": bool}}`

### 2. CSV Files for Analysis
- **`papers_to_hallmarks.csv`** (458 KB)
  - Long format: one row per paper-hallmark pair
  - Columns: Paper Title, Hallmark, Score, Is_Aging_Related
  - Best for: filtering by hallmark, visualization, pivot tables

- **`paper_summary_stats.csv`** (321 KB)
  - One row per paper with aggregate statistics
  - Columns: Paper Title, Number_of_Hallmarks, Total_Score, Primary_Hallmark, Primary_Hallmark_Score, Is_Aging_Related
  - Best for: sorting by score, finding multi-hallmark papers, identifying paper focus

- **`hallmark_statistics.csv`** (843 bytes)
  - One row per hallmark with aggregate statistics
  - Columns: Hallmark, Paper_Count, Average_Score, Max_Score, Min_Score, Total_Score
  - Best for: understanding hallmark representation, identifying research gaps

### 3. Documentation
- **`README_methodology.md`** (6.6 KB)
  - Complete methodology documentation
  - Detailed results and analysis
  - Limitations and use cases

## Quick Examples

### Example 1: Find all papers on "cellular senescence"
```python
import json

with open('papers_mapped_to_hallmarks.json', 'r') as f:
    data = json.load(f)

senescence_papers = [
    (title, info['hallmarks']['cellular senescence'])
    for title, info in data.items()
    if 'cellular senescence' in info.get('hallmarks', {})
]

# Sort by score
senescence_papers.sort(key=lambda x: x[1], reverse=True)
print(f"Found {len(senescence_papers)} papers on cellular senescence")
```

### Example 2: Find papers spanning multiple hallmarks
```python
import json

with open('papers_mapped_to_hallmarks.json', 'r') as f:
    data = json.load(f)

multi_hallmark = [
    (title, len(info['hallmarks']), info['hallmarks'])
    for title, info in data.items()
    if len(info.get('hallmarks', {})) >= 4
]

multi_hallmark.sort(key=lambda x: x[1], reverse=True)
print(f"Found {len(multi_hallmark)} papers with 4+ hallmarks")
```

### Example 3: Using CSV in Excel/R/Python pandas
```python
import pandas as pd

# Load summary statistics
df = pd.read_csv('paper_summary_stats.csv')

# Get top 10 papers by total score
top_papers = df.nlargest(10, 'Total_Score')

# Count papers by number of hallmarks
hallmark_dist = df['Number_of_Hallmarks'].value_counts().sort_index()

# Filter aging-related papers only
aging_papers = df[df['Is_Aging_Related'] == 'Yes']
```

### Example 4: Find research gaps
```python
import pandas as pd

# Load hallmark statistics
hallmarks = pd.read_csv('hallmark_statistics.csv')

# Sort by paper count to see underrepresented areas
gaps = hallmarks.sort_values('Paper_Count')
print("Least studied hallmarks:")
print(gaps[['Hallmark', 'Paper_Count']].head())
```

## Key Insights from the Data

### Coverage
- **73%** of papers are related to aging hallmarks
- **27%** have no identified hallmark associations
- **46.5%** of aging-related papers map to multiple hallmarks

### Top Hallmarks
1. **Epigenetic alterations** (1,203 papers) - Dominant focus
2. **Genome instability** (267 papers)
3. **Chronic inflammation** (267 papers)
4. **Mitochondrial dysfunction** (204 papers)

### Underrepresented Areas
1. **Dysbiosis** (7 papers) - Major research gap
2. **Intercellular communication** (31 papers)
3. **Loss of proteostasis** (39 papers)
4. **Telomere attrition** (55 papers)

### Common Combinations
Papers often combine:
- Epigenetic alterations + Genome instability (214 papers)
- Chronic inflammation + Epigenetic alterations (153 papers)
- Epigenetic alterations + Mitochondrial dysfunction (121 papers)

## Recommended Workflows

### Workflow 1: Focused Literature Review
1. Load `papers_mapped_to_hallmarks.json`
2. Filter by your hallmark(s) of interest
3. Sort by score to find most relevant papers
4. Review papers with multiple related hallmarks for comprehensive coverage

### Workflow 2: Identifying Cross-cutting Research
1. Load `paper_summary_stats.csv`
2. Filter for papers with 3+ hallmarks
3. Analyze which hallmark combinations are most common
4. Use for understanding mechanistic connections

### Workflow 3: Trend Analysis
1. If you have publication dates, join with `papers_to_hallmarks.csv`
2. Analyze temporal trends in hallmark focus
3. Identify emerging vs declining research areas

### Workflow 4: Collaboration Discovery
1. Extract author information from source data
2. Join with `papers_to_hallmarks.csv`
3. Identify researchers working on specific hallmark combinations
4. Map collaboration networks

## Quality Control

### Validation Checks Performed
✓ All 2,147 papers processed
✓ Score consistency verified (all scores > 0)
✓ No duplicate papers
✓ All hallmarks from predefined list
✓ Boolean flag consistency checked

### Known Limitations
- Some papers about aging may be classified as "not-related" if they use concepts not in the mapping dictionary
- Sirtuin-heavy dataset leads to high epigenetic alterations representation
- Scores represent cumulative concept relevance, not paper focus intensity

## Support and Feedback

### Common Issues
- **Q: Why is my paper marked "not-related to aging"?**
  - A: The paper's concepts don't match our hallmark mapping. This may indicate the paper uses different terminology or focuses on aspects not captured in the 12 hallmarks.

- **Q: Why does my paper have a low score despite being relevant?**
  - A: Scores are cumulative from concept scores. Papers with few concepts or concepts with lower PyAlex relevance scores will have lower totals.

- **Q: Can I add my own hallmark mappings?**
  - A: Yes! Edit the CONCEPT_TO_HALLMARK dictionary in `map_papers_to_hallmarks.py` and re-run the analysis.

## Next Steps

1. **Explore the data** using your preferred tools (Python, R, Excel)
2. **Validate classifications** for papers in your area of expertise
3. **Customize mappings** if needed by editing the source script
4. **Combine with metadata** (dates, authors, citations) for deeper analysis
5. **Share findings** with the research community

---

For detailed methodology and full analysis, see `README_methodology.md`
