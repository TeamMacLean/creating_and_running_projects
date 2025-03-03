Project-Oriented Bioinformatics Workflow

A well-organized bioinformatics project structure is essential for reproducibility, collaboration, and efficient development. A good project structure is based around a well organised project directory.

## Core Directory Structure
```
project_name/
├── README.md
├── analysis/
├── scripts/
├── lib/
├── data/
├── Snakemake
└── software.def
```

### Analysis Directory

The `analysis/` directory contains notebooks and documents detailing the analytical process:

	* Jupyter/R notebooks documenting analyses with embedded visualizations
	* Reports generated from analyses
	* Organized chronologically (e.g., `0001_quality_control`, `0002_alignment`)
	* Each analysis should be self-contained but reference shared components

Effective naming conventions in your analysis directory are crucial for project organization and clarity. Prefer leading `00..` chronological naming. E.G.

```
analysis/

   ├── 0001_fastqc_raw_reads.ipynb
   ├── 0002_trimming_assessment.ipynb
   ├── 0003_fastqc_trimmed_reads.ipynb
   ├── 0004_assembly_qc.Rmd
   └── 0005_annotation.Rmd
```


This structure provides a clear progression of analysis steps, is easy to follow for new team members and enforces logical workflow ordering.

#### Delineating HPC steps and Notebook steps

It is often not practical to run a step in a Notebook because of computational resource or time requirements. It is therefore necessary and legitimate to reference a Snakemake step or script run at the top of a notebook and work on the output from the prior HPC step in a notebook. E.G 

```
	0001_assembly_checking.Rmd
	
	In this document I'm going to assess the quality of assemblies from various programs on my reads. Reads were assembled
	using the snakemake step `assemble_all` and the output .fasta stored in `data/assembly_1.fa`. The snakemake step `run_busco` was automatically run and dropped in the `data/busco_assembly_1.fa`. 
```

This helps contextualise the current analysis.

#### Applying Scientific Thinking Explicitly

In the group we believe that `GoHREp` is an excellent tool for planning and expressing scientific thinking. If you didn't `GoHREp` you didn't science. Therefore in these documents `GoHREp` is the prefferred way to explain what you're about to get up to. Here's a deeper read on `GoHREp` by it's creator [https://kamounlab.medium.com/](https://kamounlab.medium.com/gohrep-and-plesi-guides-to-navigate-through-your-research-projects-ae354a21f1bc)

#### An example `GoHRep` in an `000X_my_analysis.Rmd`

Here's how it might go in practice:

```
#' # Differential Expression Analysis of Potato Genes During Late Blight Infection
#' 
#' ## GoHREP Framework
#' 
#' ### Goal
#' Identify key potato genes responding to Phytophthora infestans infection at 48 hours post-inoculation.
#' 
#' ### Hypothesis
#' Resistant potato cultivars will upregulate defense-related genes (PR proteins, WRKY transcription factors) 
#' earlier and more strongly than susceptible cultivars.
#' 
#' ### Resources
#' - Count matrices: ../data/processed/counts/potato_48hpi.csv
#' - Metadata: ../data/metadata/potato_samples.csv
#' - Potato genome annotation: resources/spud_v4.1_annotation.gff
#' 
#' ### Experimental Plan
#' 1. Filter low-count genes (<10 counts in >50% samples)
#' 2. Normalize counts using variance stabilizing transformation
#' 3. Perform DESeq2 analysis comparing resistant vs. susceptible cultivars
#' 4. Identify enriched biological pathways in DE genes
#' 5. Visualize top genes with pathogen-responsive elements

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(biomaRt)
library(readr)
library(dplyr)

# Source custom functions from lib directory
source("../scripts/deseq_functions.R")
source("../scripts/visualization_tools.R")

# Load data
counts <- read_csv("../data/processed/counts/potato_48hpi.csv", row.names = 1)
metadata <- read_csv("../data/metadata/potato_samples.csv")

# Ensure sample order matches between counts and metadata
metadata <- metadata[match(colnames(counts), metadata$sample_id),]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData = metadata,
  design = ~ cultivar_resistance
)

...
```

#### Temporal considerations

The `0001_my_analysis.ipynb` documents are not daily records, they may take less or more than a day to create/complete. They are tied to the scientific intent described in their `GoHREp` so when the aims in there are done, start a new one. Even if you have to go on holiday in the middle! 


### Snakemake 

The Snakemake file naturally complements the directory structure we discussed. Your `scripts/` directory provides the individual tools, your  `data/` and`lib/` directories contain reusable data, and Snakemake orchestrates how and when each component runs on the HPC. It formalizes the relationships between the different parts of your analysis pipeline. Most importantly, it provides reproducibility by capturing the complete workflow, from raw data to final results, in a way that can be shared, versioned, and rerun. This addresses one of the core challenges in bioinformatics: ensuring that complex multi-step analyses can be reproduced exactly. As such the Snakemake file should exist from the start and be updated incrementally to make most efficient use. You can read a bit about what it is and how to use it here: [Snakemake tutorial](https://danmaclean.github.io/snakemake/) 

#### Scope of Snakemake

Not all steps need go in the Snakemake, especially the analysis in `analysis` does not generally need to go in there. In most cases Snakemake will be good for tying together the HPC steps. 

### Scripts Directory

The `scripts/` (or `src/`) directory contains executable code for data processing:

Pipeline scripts that coordinate multiple analysis steps
Command-line tools and utilities
Batch processing scripts
Submission scripts for HPC environments
Shell scripts for environment setup and data preparation

Ideally these will be called from a master Snakemake document that starts off small and gets bigger as the project evolves, each new step being integrated into the last ones capturing the build process.

### Library Directory (`lib/`)

The `lib/` directory contains reusable data, likely from external sources, like reference genomes or annotations that you will use over and over, but basically just downloaded from external providers.
It can also contain defs for software installation or other largely static stuff. You wouldn't write in a library book, so these are files you don't write in. But perhaps you'd read them often. A useful `lib` file is a `samples_to_directory.txt` file that definitively matches the paths to large data files (see `data/`) and the sample name you would like to refer to them by.


### The `data/` directory

This directory is for data generated within the project, typically by the Snakemake file and anything from your scripts that are co-ordinated by the Snakemake and you `00X..` analysis files. Not all data can or should be stored here. Data that is transient and can easily be regenerated should be deleted when done with and data that are large or have an explicit place in data stores in the file system (sequence read data and similar) should be referred to. It's often good to use a `samples_to_directory.txt` file for this.

### Dependency Management:

Document software versions in README.md
Document daily changes in README.md
Consider containerization (Singularity) for complete reproducibility

### Version Control:

Track with Git, but not large data files
Use .gitignore to exclude data and results directories
Refer to large e.g sequencing data with their fixed path in the HPC
Use a `lib/samples_to_directory.txt` file

### Best Practices

Modularity: Design components that can be reused and combined
Documentation: Maintain comprehensive documentation in code and separate files
Testing: Implement unit tests for library components and integration tests for workflows
Configuration: Use configuration files rather than hardcoded parameters
Logging: Implement detailed logging for troubleshooting and tracking progress


