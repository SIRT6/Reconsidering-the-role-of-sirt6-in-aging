# Aging Hallmarks Mapping - Complete Output Package

## 📋 Start Here

**New users**: Read [FINAL_SUMMARY.txt](FINAL_SUMMARY.txt) for a quick overview, then see [QUICKSTART.md](QUICKSTART.md) for usage examples.

**Detailed methodology**: See [README_methodology.md](README_methodology.md)

## 📁 File Directory

### Primary Data Files
| File | Size | Description | Best For |
|------|------|-------------|----------|
| [papers_mapped_to_hallmarks.json](papers_mapped_to_hallmarks.json) | 506 KB | Complete mapping results in JSON format | Programmatic access, complete data structure |
| [papers_to_hallmarks.csv](papers_to_hallmarks.csv) | 458 KB | Long format: one row per paper-hallmark pair | Filtering by hallmark, pivot tables, visualization |
| [paper_summary_stats.csv](paper_summary_stats.csv) | 321 KB | One row per paper with aggregate stats | Sorting, identifying multi-hallmark papers |
| [hallmark_statistics.csv](hallmark_statistics.csv) | 843 B | Aggregate statistics per hallmark | Overview, identifying research gaps |

### Documentation Files
| File | Size | Description |
|------|------|-------------|
| [FINAL_SUMMARY.txt](FINAL_SUMMARY.txt) | 5.0 KB | Executive summary with key findings |
| [QUICKSTART.md](QUICKSTART.md) | 6.3 KB | Quick reference guide with code examples |
| [README_methodology.md](README_methodology.md) | 6.6 KB | Complete methodology and detailed analysis |
| [INDEX.md](INDEX.md) | This file | Navigation guide |

### Source Code
| File | Size | Description |
|------|------|-------------|
| [map_papers_to_hallmarks.py](map_papers_to_hallmarks.py) | 8.9 KB | Python script used for mapping (reproducible, customizable) |

## 🎯 Quick Access by Use Case

### "I want to find papers on a specific hallmark"
→ Use `papers_to_hallmarks.csv` and filter by Hallmark column

### "I want to analyze multi-hallmark papers"
→ Use `paper_summary_stats.csv` and filter where Number_of_Hallmarks > 1

### "I want to build my own analysis pipeline"
→ Load `papers_mapped_to_hallmarks.json` in your preferred language

### "I want to understand the methodology"
→ Read `README_methodology.md`

### "I want to modify the concept mappings"
→ Edit `map_papers_to_hallmarks.py` and re-run

### "I need quick statistics"
→ Open `hallmark_statistics.csv` or read `FINAL_SUMMARY.txt`

## 📊 Key Statistics at a Glance

- **Total papers**: 2,147
- **Aging-related**: 1,568 (73%)
- **Multiple hallmarks**: 729 papers (47% of aging-related)
- **Most common**: Epigenetic alterations (1,203 papers)
- **Least common**: Dysbiosis (7 papers)

## 🔬 The 12 Aging Hallmarks

1. Genome instability
2. Telomere attrition
3. Epigenetic alterations
4. Loss of proteostasis
5. Disabled macroautophagy
6. Deregulated nutrient sensing
7. Mitochondrial dysfunction
8. Cellular senescence
9. Stem cell exhaustion
10. Intercellular communication
11. Chronic inflammation
12. Dysbiosis

## 🚀 Getting Started Examples

### Python
```python
import json
import pandas as pd

# Load JSON data
with open('papers_mapped_to_hallmarks.json', 'r') as f:
    data = json.load(f)

# Or load CSV
df = pd.read_csv('paper_summary_stats.csv')
```

### R
```r
library(jsonlite)
library(readr)

# Load JSON
data <- fromJSON('papers_mapped_to_hallmarks.json')

# Or load CSV
df <- read_csv('paper_summary_stats.csv')
```

### Excel
Simply open any of the `.csv` files directly in Excel or Google Sheets.

## ⚙️ Output Format Details

### JSON Structure
```json
{
  "Paper Title": {
    "hallmarks": {
      "hallmark_name": score
    },
    "not-related to aging": false
  }
}
```

### CSV Columns

**papers_to_hallmarks.csv**:
- Paper Title, Hallmark, Score, Is_Aging_Related

**paper_summary_stats.csv**:
- Paper Title, Number_of_Hallmarks, Total_Score, Primary_Hallmark, Primary_Hallmark_Score, Is_Aging_Related

**hallmark_statistics.csv**:
- Hallmark, Paper_Count, Average_Score, Max_Score, Min_Score, Total_Score

## 📈 Analysis Ideas

1. **Temporal trends**: Join with publication dates to track hallmark focus over time
2. **Citation analysis**: Combine with citation data to identify influential papers per hallmark
3. **Author networks**: Map collaboration patterns across hallmark combinations
4. **Gap analysis**: Identify underexplored hallmark combinations
5. **Literature review**: Generate reading lists for specific hallmark interests

## ⚠️ Important Notes

- Scores represent cumulative concept relevance, not paper quality
- 27% of papers have no hallmark matches (may use different terminology)
- Dataset appears enriched for sirtuin-related papers
- Concept-to-hallmark mappings can be customized in the source code

## 🔧 Customization

To modify concept mappings:
1. Open `map_papers_to_hallmarks.py`
2. Edit the `CONCEPT_TO_HALLMARK` dictionary
3. Re-run: `python3 map_papers_to_hallmarks.py`

## 📞 Questions?

- **Methodology questions**: See `README_methodology.md`
- **Usage examples**: See `QUICKSTART.md`
- **Quick reference**: See `FINAL_SUMMARY.txt`

---

**Generated**: 2025-11-17  
**Dataset**: 2,147 papers from PyAlex  
**Hallmarks**: 12 predefined aging hallmarks  
**Status**: Complete ✓
