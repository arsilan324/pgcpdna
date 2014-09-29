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
	sed -E 's/^>(.*gene=([^]]*).*)$$/>\2|\1/' $< >$@

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
	gt gff3 -addintrons $< \
	|gsed '/rrn/s/mRNA/rRNA/; \
		/trn/s/mRNA/tRNA/' >$@

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
		--ircoord `<$*.ircoord`

%.gff.png: %.gff
	gt sketch $@ $<

%.all.maker.proteins.fasta %.all.maker.transcripts.fasta: %.maker.output/stamp
	fasta_merge -d $*.maker.output/$*_master_datastore_index.log

$(ref).%.p2g: %.fa $(ref).fa
	exonerate $< $(ref).fa \
		-m p2g -i -11 \
		--proteinwordlen 3 --proteinhspthreshold 13 \
		--geneticcode 11 \
		--splice5 splice0 --splice3 splice0 \
		--bestn 1 >$@

$(ref).%.p2g.gff2: %.fa $(ref).fa
	exonerate $< $(ref).fa \
		-m p2g -i -11 \
		--proteinwordlen 3 --proteinhspthreshold 13 \
		--geneticcode 11 \
		--splice5 splice0 --splice3 splice0 \
		--bestn 1 \
		--showtargetgff true --showalignment false --showvulgar false \
		>$@
		#|sed -n '/^##gff-version/,/END OF GFF DUMP/p' >$@

%.p2g.gff: %.p2g.gff2
	bin/process_exonerate_gff3 -t Protein $< >$@

# Report the annotated genes

%.gff.gene: %.gff
	sort -k1,1 -k4,4n -k5,5n $< \
	|gsed -nE '/\tgene\t/!d; \
		s/^.*Name=([^|;]*).*$$/\1/; \
		s/-gene//; \
		p' >$@
