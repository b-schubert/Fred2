__author__ = 'schubert'


import argparse, sys, os

sys.path.append('/nfs/wsi/abi/nicotin/etk/Fred2')
#this is important for netChop to work
os.environ['NETCHOP'] = '/nfs/wsi/abi/nicotin/etk/software/netchop-3-1.1/'
os.environ['TMPDIR'] = '/nfs/wsi/abi/nicotin/etk/software/netchop-3-1.1/tmp/'

import itertools
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import ACleavageSitePrediction
from Fred2.IO.FileReader import read_lines
from Fred2.EpitopeAssembly.EpitopeAssembly import EpitopeAssembly
from Fred2.CleavagePrediction import CleavageSitePredictorFactory


def to_json(assembler,out):
    with open(out, "w") as f:
        f.write('{\n"data":[\n')
        rows = []
        for s,e in assembler.instance.w_ab:
             if s != "Dummy"  and e != "Dummy":
                 rows.append('[\n\t"%s",\n\t"%s",\n\t"%.3f"\n]'%(s,e,assembler.instance.w_ab[(s,e)]))
        f.write(",\n".join(rows)+"\n]\n}")

def main():
    parser = argparse.ArgumentParser(description="Reads peptide sequences and predicts best arrangement "+
                                                 "for a string-of-beats peptide vaccine based on proteasomal cleavage prediction.")
    parser.add_argument("-i", "--input",
                        nargs="+",
                        required=True,
                        help="peptide file (one peptide per line),"
                             +" or peptide sequences as sequences (max 50)"
                        )
    input_types = parser.add_mutually_exclusive_group(required=True)
    input_types.add_argument("-pf","--pepfile",
                             action="store_true",
                             help= "Specifies the input as peptide file")
    input_types.add_argument("-p","--peptide",
                             action="store_true",
                             help= "Specifies the input as peptide sequences")
    parser.add_argument("-m", "--method",
                        default="PCM",
                        help="Specifies the Cleavage Site prediction tool to use - default PCM."
                      )
    parser.add_argument("-s", "--solver",
                        default="glpk",
                        help="Specifies the ILP solver to use (must be installed) - default glpk"
                      )
    parser.add_argument("-adv", "--advanced",
                        type=float,
                        default=0.0,
                        help="Specifies if cleavage sites within epitopes should be considdered and with what strength [0,1]"
                      )
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Specifies the output path. Results will be written to CSV")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Specifies verbose run."
                      )
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
    parser.add_argument("-oi", "--out_internal",
                        required=False,
                        default="",
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    try:
        if args.available:
            for pred, obj in ACleavageSitePrediction.registry.iteritems():
                if pred not in ["ACleavageSitePrediction", "APSSMCleavageSitePredictor"]:
                    print "Method: ",pred
                    print "Supported Length: ", " ".join(map(str, getattr(obj,  "_"+pred+"__supported_length")))
                    print
            sys.exit(0)

        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

        if args.pepfile:
            peptides = read_lines(args.input[0], type="Peptide")

        else:
            peps = list(itertools.chain(*map(lambda x: x.replace("__cr____cn__"," ").replace(","," ").replace(";", " ").split(), args.input)))
        
            peptides = [Peptide(s.upper()) for s in peps if all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s.upper())]

        if not peptides:
            raise ValueError("All specified peptides contain unaturall amino acids or the input was empty.")
            
        cleav_pred = CleavageSitePredictorFactory(args.method)
        assembler = EpitopeAssembly(peptides, cleav_pred, solver=args.solver, weight=args.advanced, threads=4, verbosity=1 if args.verbose else 0)
        result = assembler.solve()
        with open(args.out_internal, "w") as io:
            io.write("\n".join( str(p) for p in result ))
        output = args.outdir+"/epitope_assembly.csv"
        json_out = args.outdir+"/result.json"
        with open(output, "w") as out:
            out.write("Start,Stop,Score\n")
            out.write("\n".join("%s,%s,%.2f"%(s,e,assembler.instance.w_ab[(s,e)]) for s,e in assembler.instance.w_ab if s != "Dummy"  and e != "Dummy"))

        #to json
        to_json(assembler,json_out)

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
        <div class="document">
        """
    
            setting = """  <h2 class="etk-heading">Input Summary</h2>

            <table class="etk-parameterT">
                <tr>
                    <th>Cleave Site Prediction Method:</th>
                    <td>%s</td>
                </tr>
                <tr>
                    <th> Influence of cleavage sites within epitopes:</th>
                    <td>%.2f</td>
                </tr>
            </table>"""%(args.method, args.advanced)


            table="""
            <h2 class="etk-heading">Cleavage Propability between Epitopes</h2>
            <table class="etk-resultsT">

                <thead>
                    <tr>
                        <th>Start</th><th>Stop</th><th>Score</th>
                    </tr>
                </thead>"""+"".join("<tr><td>%s</td><td>%s</td><td>%f</td></tr>"%(s,e,assembler.instance.w_ab[(s,e)]) for s,e in assembler.instance.w_ab if s != "Dummy"  and e != "Dummy")+"</table><br /><a href='%s' download='epitope_assembly'>Download as CSV</a>"%os.path.basename(output)
            
            
            table2="""
            <h2 class="etk-heading">Cleavage Propability between Epitopes</h2>
            <table id="result_table" class="display">

                <thead>
                    <tr>
                        <th>Start</th><th>Stop</th><th>Score</th>
                    </tr>
                </thead>
            </table>"""
            
            
            s = ["<b>N</b><font color='blue'>-</font> 0.0 <font color='blue'>-</font><b>"+str(result[0])+"</b>"]
            for i in xrange(1,len(result)):
                start=str(result[i-1])
                end=str(result[i])
                s.append("<font color='blue'>-</font> "+"<font color='#DF363D'>{:10.2f}</font>".format(assembler.instance.w_ab[(start,end)])+" <font color='blue'>-</font><b>"+end+"</b>")
            s.append("</b><font color='blue'>-</font> 0.0 <font color='blue'>-</font><b>C</b>")
            beats="""
                    <h2 class="etk-heading">String-of-Beads Sequence</h2>
                    <br />
                    <p>%s</p>
            """%"".join(s)
            refs = """
                <h2 class="etk-heading">References</h2>
                <ol>
                   <li>Toussaint N. C, Maman Y, Kohlbacher O, Louzoun Y. (2011) Universal peptide vaccines - Optimal peptide vaccine design based on viral sequence conservation. Vaccine 47(29), 8745-8753.</li>
                   <li>Doennes, P, and Kohlbacher, O. (2005). Integrated modeling of the major events in the MHC class I antigen processing pathway. Protein Sci, 14(8), 2132-2140. doi: 10.1110/ps.051352405</li>
                   <li>Nielsen, M, Lundegaard, C, Lund, O, and Kesmir, C. (2005). The role of the proteasome in generating cytotoxic T-cell epitopes: insights obtained from improved predictions of proteasomal cleavage. Immunogenetics, 57(1-2), 33-41. doi: 10.1007/s00251-005-0781-7</li>
                </ol>
                """
            end_html = "</div></body></html>"

            html_out = ".".join(output.split(".")[:-1])+".html"
            with open(args.output, "w") as html_o:
                html_o.write(begin_html+setting+table2+beats+refs+end_html)
    except Exception as e:
        sys.stderr.write("An Error occured: "+str(e))
        sys.exit(6)
        
if __name__ == "__main__":
    main()
