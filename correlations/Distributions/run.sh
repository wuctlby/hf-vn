OPTION="-b --configuration json://cfg.json"
LOGFILE="log.txt"

o2-analysis-hf-correlator-flow-charm-hadrons-reduced $OPTION --aod-writer-json OutputDir.json --aod-file @input_derived_data.txt 1 > $LOGFILE 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi

touch ao2ds_to_merge.txt
echo AO2D.root > ao2ds_to_merge.txt
o2-aod-merger --input ao2ds_to_merge.txt --max-size 1 --output Tree.root
rm ao2ds_to_merge.txt
rm localhost*
