all:
	python setup.py build_ext --inplace
	rm -rf build
	mkdir -p ../bin
	cp *.so ../bin/
	cp sr_sam ../bin/
	cp find_best_split ../bin/
	cp combine_alignments ../bin/
	cp align.sh ../bin/
	cp align_reads_with_lz ../bin/
	cp version.py ../bin/
	rm *.c
	rm *.so

.PHONY : clean

clean:
	-rm -rf ../bin
