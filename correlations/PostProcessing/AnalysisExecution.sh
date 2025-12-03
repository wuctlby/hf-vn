#!/bin/bash

# \brief Bash script to run the azimuthal correlation analysis
# \usage ./AnalysisExecution.sh

# Create the output directory
OUTDIR="Output_AnalysisExecution"
if [ ! -d "$OUTDIR" ]; then
    mkdir Output_AnalysisExecution
fi

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_010_posDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_010_posDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_010_posDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_010_posDeta.json")
# .q
# EOF

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_1020_posDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_1020_posDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_1020_posDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_1020_posDeta.json")
# .q
# EOF

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_2050_posDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_2050_posDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_2050_posDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_2050_posDeta.json")
# .q
# EOF

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_5090_posDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_5090_posDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_5090_posDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_5090_posDeta.json")
# .q
# EOFc

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_010_negDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_010_negDeta.json")
# .q
# EOF

# Correlation fit
log_file="Output_AnalysisExecution/stdoutFitCorrel_010_negDeta.log"
root -b -l <<'EOF' 2>&1 | tee "$log_file"
Printf("[INFO] FIT CORRELATIONS");
gSystem->Exec("rm -f DhCorrelationFitter_cxx*")
gSystem->Exec("rm -f FitCorrel_C*")
gSystem->SetBuildDir(".", kTRUE)
.L DhCorrelationFitter.cxx++
.L FitCorrel.C++
.x FitCorrel.C("config_CorrAnalysis_v2_010_negDeta.json")
.q
EOF

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_1020_negDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_1020_negDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_1020_negDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_1020_negDeta.json")
# .q
# EOF

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_2050_negDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_2050_negDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_2050_negDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_2050_negDeta.json")
# .q
# EOF

# # Correlation extraction
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutExtractCorrel_5090_negDeta.log
# Printf("[INFO] EXTRACT CORRELATIONS");
# .L DhCorrelationExtraction.cxx+
# .x ExtractOutputCorrel.C("config_CorrAnalysis_v2_5090_negDeta.json")
# .q
# EOF

# # Correlation fit
# root -b -l <<EOF &> Output_AnalysisExecution/stdoutFitCorrel_5090_negDeta.log
# Printf("[INFO] FIT CORRELATIONS");
# .L DhCorrelationFitter.cxx+
# .x FitCorrel.C("config_CorrAnalysis_v2_5090_negDeta.json")
# .q
# EOF

