# Purpose
Run ConSurf for all human and mouse proteins (Download from NCBI database).

# Steps
- download protein fasta `wget file.link`
- merge all fasta `cat * > all.fasta`
- seperate to individual fasta `seperate_seq_in_fasta.py`
- qsub `get_conserv_score.sh` to server nodes, run in parallel

# Design
- for each individual fasta, check if its result file (same name in `results` folder) exist
- if exist, go to next individual fasta
- if not exist, run Consurf to generate the result file
- [in this way, you can run the jobs in multiple nodes, no need to worry about program abort]
