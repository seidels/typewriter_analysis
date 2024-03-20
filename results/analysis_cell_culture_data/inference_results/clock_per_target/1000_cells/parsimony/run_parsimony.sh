
for seed in $(seq 1 649)
do
   java -jar /Users/azwaans/all_beasts2.7/out/artifacts/all_beasts2_7_sciphy_jar/all_beasts2.7.jar -version_file /Users/azwaans/all_beasts2.7/sciphy/version.xml -seed ${seed} typewriter_13Sites_1000Cells_DataSet1_PARSIMONY.xml
done
