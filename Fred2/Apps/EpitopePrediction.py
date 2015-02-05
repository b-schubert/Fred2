__author__ = 'schubert'

import argparse, sys, os
import pandas
from tempfile import NamedTemporaryFile
sys.path.append('/nfs/wsi/abi/nicotin/etk/Fred2')

from Fred2.IO import FileReader
from Fred2.Core.Generator import generate_peptides_from_protein
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.MartsAdapter import MartsAdapter

import sys, getopt, urllib, urllib2
uniprot_url = 'http://www.uniprot.org/uniprot/'
ncbi_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='

def uniprot_createFasta(input, output):
    """
    REST adapert for UNIPROT
    :param: str input: defines the input uniprot ids
    :param: str output: defines the output file to which the querried IDs are written (in fasta)
    """
    geneList = formatToList(input)
    for gene in geneList:
        website = urllib2.urlopen(uniprot_url + gene + '.fasta')
        result = website.read()
        output.write(result + '\n')

def ncbi_createFasta(input, output):
        """
         REST adapert for RefSeq
        :param: str input: defines the input uniprot ids
        :param: str output: defines the output file to which the querried IDs are written (in fasta)
        """
        geneList = formatToList(input)

        for gene in geneList:
                website = urllib2.urlopen(ncbi_url + gene + '&rettype=fasta&retmode=text')
                result = website.read()
                output.write(result + '\n')

def formatToList(arg):
    return map(lambda x: x.replace(";","").replace(",","").strip(), arg)

#def table_to_html(table):
#    return table.to_html(classes=['display','etk-resultsT'],float_format='{0:.3f}'.format))

def table_to_html(df):
    header=["<th>Sequence</th><th>Method</th>"]
    body=[]
    table="""
    <table class="display etk-resultsT">
        <thead>
            %s
        </thead>
        <tbody>
            %s
        </tbody>
    </table>    
    """
    header.extend( "<th>%s</th>"%a.name for a in df.columns)
    for s in list(df.index.levels[0]):
        methods = list(df.index.levels[1])
	row=[]
        #row=["<tr>\n\t<td rowspan="+str(len(methods))+">"+str(s)+"</td>\n"+"<td>"+methods[0]+"</td>\n"+"\n".join("<td>%.02f</td>"%df.loc[(s,methods[0]),a] for a in df.columns)+"</tr>"]
        for m in methods:
            #print s,m
            #ro ="<tr><td>"+m+"</td>\n%s</tr>"%"\n".join("<td>%.02f</td>"%df.loc[(s,m),a] for a in df.columns)
            ro="<tr>\n\t<td>"+str(s)+"</td>\n"+"<td>"+m+"</td>\n"+"\n".join("<td>%.02f</td>"%df.loc[(s,m),a] for a in df.columns)+"</tr>"
	    row.append(ro)
        body.append("\n".join(row))
    return table%("\n".join(header),"\n".join(body))

def to_json(df, out, pep_to_prot=None, web_adress=None):
    """
        returns the result injson format for DataTable's ajax approach
    """            
    with open(out, "w") as f:
        f.write('{\n"data":[\n')
        rows = []
        for i in df.index:
            j= map(str,i)
            if pep_to_prot is None:
                row = '[\n\t"%s",\n\t"%s",\n%s\n]'%(j[0],j[1],",\n".join('\t"%0.3f"'%df.loc[i,a] for a in df.columns))
            else:
                #print "web_adress",web_adress
                if web_adress is None:
                    row = '[\n\t"%s",\n\t"%s",\n%s,\n\t"%s"\n]'%(j[0],j[1],",\n".join('\t"%0.3f"'%df.loc[i,a] for a in df.columns)," ".join(p.transcript_id for p in pep_to_prot[j[0]]))
                else:
                    row = '[\n\t"%s",\n\t"%s",\n%s,\n\t"%s"\n]'%(j[0],j[1],",\n".join('\t"%0.3f"'%df.loc[i,a] for a in df.columns)," ".join(web_adress%(p.transcript_id,p.transcript_id) for p in pep_to_prot[j[0]]))
            #f.write("[\n")
            #for a in df.columns:
            #    row+='\t"'+str(df.loc[i,a])+'",\n'
            #row +="\n]"
            rows.append(row)
        f.write(",\n".join(rows))
        f.write('\n]\n}')
        
def to_output_table(df,out, pep_to_prot=None):
    """
        generates an output table of the predictions which can be used as input for optitope
    """
    with open(out, "w") as f:
        f.write("Alleles:\t"+"\t".join(a.name for a in df.columns)+"\n")
        method = df.index.levels[1][0]
        for pep in df.index.levels[0]:
            if pep_to_prot is None:
                f.write(str(pep)+"\t"+"\t".join("%.3f"%df.loc[(pep,method),a] for a in df.columns)+"\n")
            else:
                f.write(str(pep)+"\t"+"\t".join("%.3f"%df.loc[(pep,method),a] for a in df.columns)+"\t"+"\t".join(p.transcript_id for p in pep_to_prot[pep])+"\n")

def main():
    parser = argparse.ArgumentParser(description="Reads protein or peptide sequences and predicts peptides "+
                                                 "for a specified prediction method and HLA alleles.")
    parser.add_argument("-i", "--input",
                        nargs="+",
                        required=True,
                        help="Input data can be RefSeq ID, UniProt ID, fasta file, peptide file (one peptide per line),"
                             +" or peptide sequences as sequences (max 50)"
                        )
    input_types = parser.add_mutually_exclusive_group(required=True)
    input_types.add_argument("-r","--refseq",
                             action="store_true",
                             help= "Specifies the input as RefSeq IDs")
    input_types.add_argument("-u","--uniprot",
                             action="store_true",
                             help= "Specifies the input as UniProt IDs")
    input_types.add_argument("-f","--fasta",
                             action="store_true",
                             help= "Specifies the input as protein (multi-)Fasta file")
    input_types.add_argument("-pf","--pepfile",
                             action="store_true",
                             help= "Specifies the input as peptide file")
    input_types.add_argument("-p","--peptide",
                             action="store_true",
                             help= "Specifies the input as peptide sequences")
    parser.add_argument("-a", "--alleles",
                        nargs="+",
                        required=True,
                        help="Specifies for which alleles prediction should be made. " +
                             "Input either can be alleles as string (new nomenclature), or a file with one allele per line.")
    allele_types = parser.add_mutually_exclusive_group(required=True)
    allele_types.add_argument("-af", "--allelefile",
                               action="store_true",
                               help="Specifies the allele input as allele file.")
    allele_types.add_argument("-as", "--allelestring",
                               action="store_true",
                               help="Specifies the allele input as allele string.")
    parser.add_argument("-m", "--method",
                       required=True,
                       nargs="+",
                       help="Specifies the method used for prediction.")
    parser.add_argument("-l", "--length",
                        required=False,
                        type=int,
                        nargs="+",
                        default=[9],
                        help="Specifies the length of the peptides (default=9).")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Specifies the output path. Results will be written to CSV")
    parser.add_argument("-og", "--output_galaxy",
                        required=True,
                        help="galaxy specific output_galaxy")
    parser.add_argument("-am", "--available",
                        required=False,
                        action="store_true",
                        help="Returns all available methods and their allele models.")

    #COMMENT: These options are hidden and only used for ETK2
    parser.add_argument("-html", "--html",
                        required=False,
                        action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("-od", "--outdir",
                        required=False,
                        default="",
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.available:
        for pred, obj in AEpitopePrediction.registry.iteritems():
            if pred not in ["AEpitopePrediction", "APSSMEpitopePredictor", "ANetMHC", "ASVMEpitopePrediction"]:
                print "Method: ",pred
                print "Supported Alleles: ", " ".join(getattr(obj, "_"+pred+"__alleles" ))
                print "Supported Length: ", " ".join(map(str, getattr(obj,  "_"+pred+"__supported_length")))
                print
        sys.exit(0)
    #try:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    #print args
    '''
    Parser Input
    '''
    #RefSeq
    if args.refseq:
        tmp_fasta = NamedTemporaryFile(delete=False)
        ncbi_createFasta(args.input, tmp_fasta)
        tmp_fasta.close()
        proteins = FileReader.read_fasta(tmp_fasta.name, type="Protein",id_position=3)
        peptides = []
        for length in args.length:
            peptides.extend(generate_peptides_from_protein(proteins, length))	
        os.remove(tmp_fasta.name)

    #UniProt
    elif args.uniprot:
        tmp_fasta = NamedTemporaryFile(delete=False)
        #print "ARGS input", args.input
        uniprot_createFasta(args.input, tmp_fasta)
        tmp_fasta.close()
        proteins = FileReader.read_fasta(tmp_fasta.name, type="Protein",id_position=1)
        peptides = []
        for length in args.length:
            peptides.extend(generate_peptides_from_protein(proteins, length))
        os.remove(tmp_fasta.name)

    #fasta protein
    elif args.fasta:
        proteins = FileReader.read_fasta(args.input[0], type="Protein")
        peptides = []
        for length in args.length:
            peptides.extend(generate_peptides_from_protein(proteins, length))
    elif args.pepfile:
        peptides = FileReader.read_lines(str(args.input[0]), type="Peptide")
        
    elif args.peptide:
        if len(args.input) == 1:
            peps = args.input[0].replace("__cr____cn__", " ").replace(",", " ").replace(";"," ").split()
        else:
            peps = [p.replace(",", "").replace(";","") for p in args.input ]
        peptides = [Peptide(s.upper()) for s in peps if all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s.upper())]

    #read in alleles
    if args.allelefile:
        alleles = FileReader.read_lines(args.alleles[0], type="Allele")
    else:
        al = map(lambda x: x.upper(), args.alleles[0].split(","))
        alleles = [Allele( a[0]+"*"+a[1:3]+":"+a[3:] if "DR" not in a else a[:4]+"*"+a[4:6]+":"+a[6:]) for a in al]

    #print peptides 
    #print alleles
    if not peptides:
        raise ValueError("All specified peptides contain unaturall amino acids or the input was empty.")
        
    methods = args.method[0].split(",")
    supported_alleles = EpitopePredictorFactory(methods[0]).supportedAlleles
    for m in methods[1:]:
        supported_alleles = supported_alleles | EpitopePredictorFactory(m).supportedAlleles
        
    if supported_alleles.intersection(set(a.name for a in alleles)) == 0:
        print "Selected methods do not support a single allele that has been selected."
        sys.exit(6)
    
    result = []
    for m in methods:
        try:
            result.append(EpitopePredictorFactory(m).predict(peptides, alleles))
        except TypeError:
            continue
    
    if not result:
        Exception("No prediction results generated. Please check your configurations.")
        sys.exit(-1)
        
    r_df = result.pop()
    for r in result:
        r_df_a, r_a = r_df.align(r, fill_value=0)
        r_df = r_df_a + r_a

    output = args.outdir + "/epitope_prediction_result.csv"
    json_out = args.outdir+"/result.json"
    with open(output, "w") as out:
        r_df.to_csv(out)
    
    
    #print "Building result tables"
    pep_to_prot = {}
    #print list(r_df.index.levels[0])
    for p in list(r_df.index.levels[0]):
        #print p, p.proteins
        for pr in p.proteins.values():
            pep_to_prot.setdefault(str(p), []).append(pr)
    #print pep_to_prot
    if pep_to_prot:
        if args.refseq:
            to_json(r_df,json_out,pep_to_prot=pep_to_prot,web_adress="<a href='http://www.ncbi.nlm.nih.gov/protein/%s'>%s</a>")
        elif args.uniprot: 
            to_json(r_df,json_out,pep_to_prot=pep_to_prot,web_adress="<a href='http://www.uniprot.org/uniprot/%s'>%s</a>")
        else:
            to_json(r_df,json_out,pep_to_prot=pep_to_prot)
        to_output_table(r_df,args.output_galaxy,pep_to_prot=pep_to_prot)
    else:
        to_json(r_df,json_out)
        to_output_table(r_df,args.output_galaxy)



    #print r_df
    #generate Galaxy HTML output
    if args.html:
        begin_html = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <link rel="stylesheet" href="/static/style/etk.css" type="text/css" />
    <link rel="stylesheet" href="/static/style/datatable/css/jquery.dataTables.css" type="text/css" />
    <link rel="stylesheet" href="/static/style/datatable/dataTables.fixedColumns.css" type="text/css" />
    <link rel="stylesheet" href="/static/style/datatable/css/tabletools/dataTables.tableTools.css" type="text/css" />
    <script type="text/javascript" src="/static/scripts/libs/jquery/jquery.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/jquery/jquery-ui.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/jquery/datatable/media/js/jquery.dataTables.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/jquery/datatable/media/js/fnFakeRowspan.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/Scroller/js/dataTables.scroller.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/FixedColumns/js/dataTables.fixedColumns.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/TableTools/js/dataTables.tableTools.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/etk.js"></script>
</head>
<body>
    """

        setting = """  <div id="document"><h2 class="etk-heading">Input Summary</h2>
        <table class="etk-parameterT">
            <tr>
                <th>Predicton Methods:</th>
                <td>%s</td>
            </tr>
            <tr>
                <th>Peptide Length:</th>
                <td>%s</td>
            </tr>
        </table>
        <h2 class="etk-heading">Epitope Prediction Results</h2>
            <p><b>Notice:</b> The returned values represent prediction scores. They are not comparable between methods!<br/>
            Typically, prediction with higher scores mean stronger binding.</p><br/>
		  """%(" ".join(args.method),  " ".join(str(l) for l in args.length ))

        result_tables="""
            <div  id="tabs">
                <ul>
                    %s
                </ul>
                %s
            </div>  
        """

        table="""

        <!-- <input id="etk-search" placeholder="  filter"> -->
        <table id="result_table" class="display">

            <thead>
                <tr>
                    <th>Peptide</th><th>Method</th>"""+"".join("<th>%s</th>"%str(a) for a in r_df.columns)
        if pep_to_prot:
            table +="<th>Origin</th>" 
            
        table+="""
                </tr>
            </thead>
        </table>"""

        refs = """
        <h2 class="etk-heading">References</h2>
        <ol>
        <li>Pruitt K.D. , Tatusova T. , Maglott D.R. (2007) NCBI reference sequences (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins. Nucleic Acids Res 35:D61-D65</li>
        <li>The UniProt Consortium The Universal Protein Resource (UniProt). Nucleic Acids Res 35:D193-D197 (2007)</li>
        <li>Rammensee H. , Bachmann J. , Emmerich N.P. , Bachor O.A. , Stevanovic S. (1999) SYFPEITHI: database for MHC ligands and peptide motifs. Immunogenetics 50:213-219.</li>
        <li>Parker K.C. , Bednarek M.A. , Coligan J.E. (1994) Scheme for ranking potential HLA-A2 binding peptides based on independent binding of individual peptide side-chains. J Immunol 152:163-175.</li>
        <li>Doennes P., Kohlbacher O. (2006) SVMHC: a server for prediction of MHC-binding peptides. Nucleic Acids Res 34:W194-W197.</li>
        Lundegaard C, Lamberth K, Harndahl M, Buus S, Lund O, Nielsen M. (2008) NetMHC-3.0: accurate web accessible predictions of human, mouse and monkey MHC class I affinities for peptides of length 8-11. Nucleic Acids Res. 1;36(Web Server issue):W509-12.</li>
        <li>Nielsen M, Lundegaard C, Blicher T, Lamberth K, Harndahl M, et al. (2007) NetMHCpan, a Method for Quantitative Predictions of Peptide Binding to Any HLA-A and -B Locus Protein of Known Sequence. PLoS ONE 2(8): e796. doi: 10.1371/journal.pone.0000796</li>
        <li>Nielsen, M. and Lund, O. (2009) NN-align. An artificial neural network-based alignment algorithm for MHC class II peptide binding prediction. BMC bioinformatics, 10, 296.</li>
        <li>Karosiene E, Rasmussen M, Blicher T, Lund O, Buus S, and Nielsen M. (2013) NetMHCIIpan-3.0, a common pan-specific MHC class II prediction method including all three human MHC class II isotypes, HLA-DR, HLA-DP and HLA-DQ. Immunogenetics.</li>
        <li>Toussaint N. C, Feldhahn M, Ziehm M, Stevanovic M, and Kohlbacher O. (2011) T-cell epitope prediction based on self-tolerance. Proc. ICIW.</li>
        <li>Zhang L, Chen Y, Wong H-S, Zhou S, Mamitsuka H, et al. (2012) TEPITOPEpan: Extending TEPITOPE for Peptide Binding Prediction Covering over 700 HLA-DR Molecules. PLoS ONE 7(2): e30483. doi: 10.1371/journal.pone.0030483</li>
        </ol>
        """

        end_html = """
        </div>
        </body></html>"""

   
        with open(args.output, "w") as html_o:
            #r_df.index = r_df.index.map(lambda x: (str(x[0]),str(x[1])))
            #r_df.index = pandas.MultiIndex.from_tuples([tuple((str(p),m)) for p,m in r_df.index], names=['Seq','Method']) 
            #html_o.write(begin_html+setting+r_df.to_html(classes=['etk-resultsT'], float_format='{0:.3f}'.format)+end_html)
            html_o.write(begin_html+setting+table+refs+end_html)
    #except Exception as e:
    #    sys.stderr.write( "An error occured "+str(e))
    #    sys.exit(6)
if __name__ == "__main__":
    main()
