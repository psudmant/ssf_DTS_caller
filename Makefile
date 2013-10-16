ALL: 
	python setup.py build_ext --inplace --include-dirs /net/eichler/vol7/home/psudmant/EEE_Lab/projects/common_code//ssf_DTS_caller/heap/
	
	#:/net/eichler/vol7/home/psudmant/EEE_Lab/projects/common_code/ssf_DTS_caller/ssf_DTS_caller/c_edge/ 
	
clean:
	rm ./c_hierarchical_edge_merge.cpp ./traverse_contours.c ./c_hierarchical_edge_merge.so ./traverse_contours.so get_windowed_variance.c get_windowed_variance.so 
	
