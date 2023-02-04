# ConSurf-debug
Stand Alone version of ConSurf with detailed tutorial - (installation + database + usage)

# About ConSurf
- https://consurf.tau.ac.il/consurf_index.php
- stand alone python wrapper: https://consurf.tau.ac.il/STANDALONE/stand_alone_consurf-1.00.rar
- ConSurf DB: https://consurfdb.tau.ac.il/ [not useful, file cannot be downloaded]
- ConSurf Github: https://github.com/Rostlab/ConSurf [not update anymore since 2015. Just take a glimpse, we won't use it]

# installation
1. download python script (wrapper) from https://consurf.tau.ac.il/STANDALONE/stand_alone_consurf-1.00.rar

2. download Rate4Site from https://www.tau.ac.il/~itaymay/cp/rate4site.html
- ftp://rostlab.org/rate4site/
```
tar zxvf rate4site-3.0.0.tar.gz
cd rate4site-3.0.0
./configure --prefix=/home/YOUR/PATH
make
make install
```
Note1: It could be successfully installed, but with some bugs like below. We can ignore it. 
Note2: The default algorithm `-im` is broken, we must use `-ib`, so we will add `--Maximum_Likelihood` when running `python stand_alone_consurf.py`. 
error:
rate4site: errorMsg.cpp:41: static void errorMsg::reportError(const string&, int): Assertion `0' failed.
Aborted (core dumped)

3. create conda env [see instructions.txt file]
```
# [official suggestion] Load modules: module load python/python-3.8 hmmr/hmmr-3.1b2 clustalw/2.1
conda create --name consurf python=3.8
conda activate consurf
conda install -c bioconda hmmer=3.1 clustalw=2.1 cd-hit mafft=7 prank muscle=3
# additional python packages
pip install biopython requests
```

4. download the databases you need
```
# I just need this, 67G after decompression
wget -b https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gzip -d uniref90.fasta.gz
```

4. modify configuration file of 
`GENERAL_CONSTANTS.py`: the path of all tools and databases are stored in this file.
```
CD_HIT_DIR = "/home/zz950/softwares/miniconda3/envs/consurf/bin/cd-hit"
UNIREF90_DB_FASTA = "/home/zz950/softwares/ConSurf/db/uniref90.fasta"
```
some tool path need to be modified in python script, `stand_alone_consurf.py`
```
PRT_JAR_FILE = "/home/zz950/softwares/ConSurf/prottest-3.4.2/prottest-3.4.2.jar"
```

5. Finally, run the main python script `stand_alone_consurf.py`
```
# these are the minimun parameters, don't remove any of them.
python stand_alone_consurf.py --algorithm HMMER --Maximum_Likelihood --seq /home/zz950/softwares/ConSurf/test.fasta --dir /home/zz950/softwares/ConSurf/test
```

6. Compare results from stand alone and web server versions.
Confirmed. The same.

# others
- all log info is in `test/log.txt` file. The full command and error info. It's useful for debug.
- a similar tool using deep learning. `vespa_emb Input_protein_seq -o data/embeddings.h5  --prott5_weights_cache data/cache` [don't like it!]
- One module (-im) of rate4site is broken. `Likelihood after optimization is -0x1.c39f46ba773c5p+14`. The log(likelihoodsV[pos]) is negative and the script will stop. So we will use (-ib) module [rate inference method].
```
Assertion failed: (log(likelihoodsV[pos])>0.0), function computeML_siteSpecificRate, file siteSpecificRate.cpp, line 27.
Abort trap: 6
```

# How I debug?
- see `log.txt` file, check where the program stoped;
- use sublime, search the folder to locate the stop code;
- insert print or log code to output the logic of the program;
- soon, you will find the bug.

