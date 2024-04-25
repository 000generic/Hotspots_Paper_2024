#! python

input_fasta = open( 'output/17-AGS_blastp-rgs-TRPs_X_species16.aa', 'r' )

species_seqid_seqs = {}
for data in input_fasta:
    if data[ 0 ] == '>':
        header = data[1:-1]
        info = header.split( '__')
        if len( info ) >1:
            gspp = info[ 0 ]
            if gspp not in species_seqid_seqs.keys():
                species_seqid_seqs[gspp] = {}
            species_seqid_seqs[ gspp ][ header ] = ''
            header_test = True
        else:
            header_test = False
    else:
        if header_test:
            species_seqid_seqs[ gspp ][ header ] += data.strip()

for gspp in species_seqid_seqs.keys():
    new_fasta = 'output/18-' + gspp + '-all_TRPs.aa'
    output_fasta = open( new_fasta, 'w' )

    for header in species_seqid_seqs[ gspp ]:
        sequence = species_seqid_seqs[ gspp ][ header ]
        output_fasta.write( '>' + header  + '\n' + sequence + '\n')

    output_fasta.close()

input_fasta.close()
