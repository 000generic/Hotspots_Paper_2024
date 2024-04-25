#! python

input_reports = open('output/1-blast-reports-self_X_self', 'r')
output_keep = open('output/2-species16_ALL_X_blast-hits-better-than-1e-60-KEEP', 'w')
output_drop = open('output/2-species16_ALL_X_blast-hits-worse-than-1e-60-DROP', 'w')

# ../../../blast/output/blast-report-1e-3-Aplysia_californica_X_Aplysia_californica
# Aplysia_californica__LOC118478065__XP_035826908.1       Aplysia_californica__LOC118478065__XP_035826908.1       100.000 148     0       0       1       148     1       148     2.13e-106       299

for report in input_reports:
    report = report.strip()
    input_report = open(report, 'r')
    info_report = report.split('/')[ -1 ].split( '-' )
    genus = info_report[ 4 ].split( '_' )[ 0 ]
    species = info_report[ 4 ].split( '_' )[ 1 ]
    gspp = genus + '_' + species
    print( gspp )  # monitoring progress

    for blast in input_report:
        blast = blast.strip()
        info_blast = blast.split('\t')
        query = info_blast[0]
        subject = info_blast[1]
        evalue = float(info_blast[-2])

        if evalue == 0 or evalue <= float( 1e-60) :
            output_keep.write( gspp + '\t' + query + '\t' + subject + '\t' + str(evalue) + '\n' )
        else:
            output_drop.write( gspp + '\t' + query + '\t' + subject + '\t' + str(evalue) + '\n' )
    input_report.close()

input_reports.close()
output_keep.close()
output_drop.close()

