all:
	make test n=complete o=+
	make test n=complete o=-
	make test n=incomplete o=+
	make test n=incomplete o=-
	make test n=split o=++
	make test n=split o=+-
	make test n=split o=-+
	make test n=split o=--

test:
	./../bin/sr_sam -t 2 -1 reference/reference.fa inputs/${n}.${o}.fq.gz \
	> tmp.sam
	diff tmp.sam outputs/${n}.${o}.sam
	rm tmp.sam
