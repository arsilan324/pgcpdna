name=pg29-plastid
ref=NC_021456

# Picea abies chloroplast complete genome
edirect_query='Picea abies[Organism] chloroplast[Title] complete genome[Title] RefSeq[Keyword]'

all: $(name).gff.gene $(name).gbk.png

clean:
	rm -f $(name).orig.gff $(name).gff $(name).orig.gbk $(name).gbk $(name).gbk.png

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Fetch data from NCBI

cds_aa.orig.fa cds_na.orig.fa: %.fa:
	esearch -db nuccore -query $(edirect_query) \
		|efetch -format fasta_$* >$@

cds_aa.fa cds_na.fa: %.fa: %.orig.fa
	sed -E 's/^>(.*gene=([^]]*).*)$$/>\2|\1/' $< |seqtk seq >$@

# asn,faa,ffn,fna,frn,gbk,gff,ptt,rnt,rpt,val
plastids/%:
	mkdir -p plastids
	curl -fsS http://ftp.cbi.pku.edu.cn/pub/database/Genome/Chloroplasts/$@ >$@

%.faa: plastids/%.gbk
	bin/gbk-to-faa <$< >$@

%.frn: plastids/%.frn
	sed 's/^>.*\[gene=/>/;s/\].*$$//' $< >$@

%.maker.output/stamp: maker_opts.ctl %.fa $(ref).frn cds_aa.fa
	maker -fix_nucleotides
	touch $@

%.orig.gff: %.maker.output/stamp
	gff3_merge -s -g -n -d $*.maker.output/$*_master_datastore_index.log >$@

%.gff: %.orig.gff
	gsed '/rrn/s/mRNA/rRNA/; \
		/trn/s/mRNA/tRNA/' $< \
	|gt gff3 -addintrons - >$@

%.orig.gbk: %.gff %.fa
	bin/gff_to_genbank.py $^
	mv $*.gb $@

%.gbk: %-header.gbk %.orig.gbk
	(cat $< && sed -En '/^FEATURES/,$${ \
		s/Name=/gene=/; \
		s/gene="([^|"]*)\|[^"]*"/gene="\1"/; \
		p;}' $*.orig.gbk) >$@

%.gbk.png: %.gbk %.ircoord
	drawgenemap --format png --infile $< --outfile $< \
		--gc --ircoord `<$*.ircoord`

%.gff.png: %.gff
	gt sketch $@ $<

# Report the annotated genes

%.gff.gene: %.gff
	sort -k1,1 -k4,4n -k5,5n $< \
	|gsed -nE '/\tgene\t/!d; \
		s/^.*Name=([^|;]*).*$$/\1/; \
		s/-gene//; \
		p' >$@

# Six-frame translation

%-orf.fa: %.fa
	sixpack -table 11 -sequence $< -outfile $*-orf.sixpack -outseq /dev/stdout |seqtk seq >$@

%.fa.fm: %.fa
	abyss-index -a '-ACDEFGHIKLMNPQRSTVWXY' $<

$(ref).%.abyss-map.sam: %.fa $(ref).fa.fm
	abyss-map --no-rc -k6 $< $(ref).fa >$@
