./EdGen -i input_test2.dat
root -l -q run_analysis.C
./EdGen -i input_test.dat
root -l -q run_newAnalysis.C
./EdGen -i input_test5.dat
root -l -q run_analysis_5.C
./EdGen -i input_test6.dat
root -l -q run_analysis_6.C
./EdGen -i input_test7.dat
root -l -q run_analysis_7.C
./EdGen -i input_t_test1.dat
root -l -q run_analysis_t_test1.C
./EdGen -i input_t_test2.dat
root -l -q run_analysis_t_test2.C
./EdGen -i input_t_test3.dat
root -l -q run_analysis_t_test3.C
./EdGen -i input_t_test4.dat
root -l -q run_analysis_t_test4.C
root -l  analysis_output.root analysis_new_output.root analysis_5_output.root analysis_6_output.root analysis_7_output.root analysis_t_test1_output.root  analysis_t_test2_output.root  analysis_t_test3_output.root  analysis_t_test4_output.root 
