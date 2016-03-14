./EdGen -i input_test2.dat
root -l -q run_analysis.C
./EdGen -i input_test.dat
root -l -q run_newAnalysis.C
./EdGen -i input_test5.dat
root -l -q run_analysis_5.C
./EdGen -i input_test6.dat
root -l -q run_analysis_6.C
root -l  analysis_output.root analysis_new_output.root analysis_5_output.root analysis_6_output.root