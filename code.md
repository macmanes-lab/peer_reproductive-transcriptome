### Rcorrector (v1.0.1)
```
perl /share/Rcorrector/run_rcorrector.pl -k 31 -t 30 \
-1 testes.R1.fastq,epidiymus.R1.fastq,vasdefs.R1.fastq \
-2 testes.R2.fastq,epidiymus.R2.fastq,vasdefs.R2.fastq
```

### Trinity and Trimmomatic  (Trinity v2.1.1)
```
Trinity --SS_lib_type RF --seqType fq --max_memory 40G --trimmomatic --CPU 30 --full_cleanup --output reproductive_trinity \
--left testes.R1.cor.fq,epidiymus.R1.cor.fq,vasdefs.R1.cor.fq \
--right testes.R2.cor.fq,epidiymus.R2.cor.fq,vasdefs.R2.cor.fq \
--quality_trimming_params "ILLUMINACLIP:/share/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

### BUSCO (v1.1b1)
```
python3 /share/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 16 -l /share/BUSCO_v1.1b1/vertebrata \
-o reproductive_busco -in reproductive_trinity.Trinity.fasta
```

### Transrate (v1.0.1) on Original Trinity Assembly
```
transrate --output reproductive_transrate -t 16 \
--assembly reproductive_trinity.Trinity.fasta \
--left testes.R1.cor.fq,epidiymus.R1.cor.fq,vasdefs.R1.cor.fq \
--right testes.R2.cor.fq,epidiymus.R2.cor.fq,vasdefs.R2.cor.fq \
--reference /mnt/data3/lauren/NEWtranscriptome/mousePEP/Mus_musculus.GRCm38.pep.all.fa
```

### Kallisto (v0.42.4) INDEX

```
kallisto index -i reproductive.idx reproductive_trinity.Trinity.fasta
```

### Kallisto QUANT
```
kallisto quant -t 32 -i reproductive.idx -o reproductive.output testes.R1.cor.fq testes.R2.cor.fq epidiymus.R1.cor.fq epidiymus.R2.cor.fq vasdefs.R1.cor.fq vasdefs.R2.cor.fq
```

### Salmon (v0.5.1) INDEX
```
/share/salmon-0.5.1/build/src/salmon index -t reproductive_trinity.Trinity.fasta -i reproductive_salmon.idx --type quasi -k 31
```

### Salmon (v0.5.1) QUANT 
```
/share/salmon-0.5.1/build/src/salmon quant -p 32 -i reproductive_salmon.idx -l IU -1 testes.R1.cor.fq epidiymus.R1.cor.fq vasdefs.R1.cor.fq -2 testes.R2.cor.fq epidiymus.R2.cor.fq vasdefs.R2.cor.fq -o salmon_reproductive
```
### dammit v0.2.7.1

```
dammit databases --install --database-dir /mnt/data3/macmanes/dammit_databases/ --full --busco-group vertebrata \
dammit annotate /mnt/data3/lauren/NEWtranscriptome/RERUNtpmHALF_transrate/reproductive_RERUNhighexpHALF.trinity.Trinity/good.reproductive_RERUNhighexpHALF.trinity.Trinity.fasta --busco-group vertebrata --n_threads 36 --database-dir /mnt/data3/macmanes/dammit_databases/ --full
```

### non-coding RNA (ncRNA) from Ensembl (Mus musculus)

```
wget ftp://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz
gzip -d Mus_musculus.GRCm38.ncrna.fa.gz
makeblastdb -in /mnt/data3/lauren/NEWtranscriptome/ncMouse/Mus_musculus.GRCm38.ncrna.fa -out ncrna -dbtype nucl

blastn -query /mnt/data3/lauren/NEWtranscriptome/ncMouse/reproductive.annotated.header.fixed.fasta -db /mnt/data3/lauren/NEWtranscriptome/ncMouse/ncrna -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > ncrna.blastn.outfmt6
```

### TCDB blastx
```
wget http://www.tcdb.org/public/tcdb
makeblastdb -in /mnt/data3/lauren/NEWtranscriptome/ncMouse/tcdb -out tcdb -dbtype prot
blastx -query /mnt/data3/lauren/NEWtranscriptome/ncMouse/reproductive.annotated.header.fasta \
-db /mnt/data3/lauren/NEWtranscriptome/ncMouse/tcdb \
-max_target_seqs 1 \
-outfmt '6 qseqid pident evalue stitle' \
-evalue 1e-5 -num_threads 10 | tee NEWtcdb.txt
```
















