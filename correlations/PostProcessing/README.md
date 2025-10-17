# Post-processing macros for angular correlation part

## Table of contents
1. [Code description](#project_overview)
2. [Execution](#execution)

---

## Code description <a name="project_overview"></a>

### Project description

This project provides a framework to perform fits on angular correlation distributions. The code is organized into two main components:
  * Angular correlation extraction:
      * ExtractOutputCorrel.C
      * DhCorrelationExtraction.cxx
      * DhCorrelationExtraction.h
  * Angular correlation distribution fit:
      * FitCorrel.C
      * DhCorrelationFitter.cxx
      * DhCorrelationFitter.h
Each component handles a specific stage of the post-processing workflow, from data extraction to fitting of the final correlation distributions.

## Execution <a name="execution"></a>
The configuration parameters are defined in dedicated config.json files, organized by centrality bin.

To run the full analysis, execute:

<pre> ```bash ./AnalysisExecution.sh ``` </pre>

To merge the output files corresponding to positive and negative Δη, use:

<pre> ```bash ./run_hadd.sh ``` </pre>