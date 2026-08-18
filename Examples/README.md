# RNA--RNA interactions examples

## KHP

### HP1:HP2

rikin_pipeline.py --fastaA HP1.fasta --fastaB HP2.fasta -c khp.cfg -o KHP --name HP1-HP2

### HP3:HP2

rikin_pipeline.py --fastaA HP3.fasta --fastaB HP2.fasta -c khp.cfg -o HP3-HP2


## MicA-MicA homo-dimer
rikin_pipeline.py --fastaA MicA.fasta --seqB MicA.fasta -c micA_micA.cfg -o MicA_MicA

## MicA-OmpA hetero-dimer
rikin_pipeline.py --fastaA MicA.fasta --seqB OmpA.fasta -c ompA_micA.cfg -o MicA_OmpA
