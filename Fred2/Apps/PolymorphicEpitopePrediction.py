"""
SNPv3 implementation with Fred2 as backend and BioMart as Database.
Does also support INDELS and Frame Shifts if so deisired

Input can be:
1) RefSeq, Ensamble, Uniprot IDs
2) vcf File
3) Polymorphic Fasta file


Annotation of variants is performed with ANNOVAR based on ensamble coordinates
"""
import sys
import os
import tempfile
import subprocess
import argparse
import itertools as itr

sys.path.append("/abi-projects/etk/Fred2/")

from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.IO import FileReader

from Fred2.Core.Base import COMPLEMENT
from Fred2.Core import Generator
from Fred2.Core.Variant import VariationType
from Fred2.Core.Allele import Allele
from Fred2.EpitopePrediction import EpitopePredictorFactory

ANOVA_PREP = "/abi-projects/etk/software/annovar/convert2annovar.pl -format vcf4 -withzyg %s > %s "
ANOVA = "/abi-projects/etk/software/annovar/annotate_variation.pl -out %s -build hg19 %s /abi-projects/etk/software/annovar/humandb/ -dbtype ensGene" #ensGene, refGene

def generate_proteins_from_polymorphic_fasta(poly_fasta):
	curr_seq = []
	curr_id = ""
	curr_vars = {}
	prots = []
	with open(poly_fasta, "r") as f:
		for l in f:
			if l.startswith(">"):
				if "|[" not in l.replace(" ",""):
					print "Specified file is not in polymorphic fasta format. See HELP for more information."
					sys.exit(-1)
				
				if curr_id == "":
					sp = l.split("|")
					curr_id = sp[0].strip()
					curr_vars = { int(var.split(":")[0]): var.split(":")[1].split(",")  for var in sp[1].replace("[","").replace("]","").split(";")}
				else: #generate Protein with Vars
					tmp_seqs = {"".join(curr_seq):[]}
					pos = curr_vars.keys()
					vars = curr_vars.values()
					for i,comb in enumerate(itr.product(*vars)):
						tmp = tmp_seqs[:]
						va = {}
						t_id=curr_id+"_"+str(i)
						for ps,v in itr.izip(pos,comb):
							tmp[ps-1]=v[0]
							va.append(Variant(random.randint(),VariationType.SNP,0,0,tmp_seq[ps-1],v[0],{t_id:MutationSyntax(t_id,0,ps*3+1,pos,"","")}))
						prots.append(Protein(tmp,t_id,t_id,_vars=va))
			else:
				curr_seq.append(l.strip().upper())
		tmp_seqs = {"".join(curr_seq):[]}
		pos = curr_vars.keys()
		vars = curr_vars.values()
		for i,comb in enumerate(itr.product(*vars)):
			tmp = tmp_seqs[:]
			vars = {}
			t_id=curr_id+"_"+str(i)
			for ps,v in itr.izip(pos,comb):
				tmp[pos-1]=v[0]
				vars.append(Variant(random.randint(),VariationType.SNP,0,0,tmp_seq[pos-1],v[0],{t_id:MutationSyntax(t_id,0,pos*3+1,pos,"","")}))
			prots.append(Protein(tmp,t_id,t_id,_vars=vars))
	return prots

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
				trans = set([ p.transcript_id.split(":")[0] for p in  pep_to_prot[j[0]]])
				vars = {}
				for tID, vs in i[0].vars.iteritems():
						for x in vs:
							if x.id.startswith("rs"):
								vars.setdefault(tID.split(":")[0], set()).add("<a href='http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s'>%s</a>"%(x.id,x.id))
							else:
								vars.setdefault(tID.split(":")[0], set()).add("%s%i%s"%(x.ref,x.genomePos,x.obs))
						
				if web_adress is None:
					row = '[\n\t"%s",\n\t"%s",\n%s,\n\t"%s",\n\t"%s"\n]'%(j[0],j[1],",\n".join('\t"%0.3f"'%df.loc[i,a] for a in df.columns),"<br/>".join(t for t in trans),
					 "<br/>".join( " ".join(vars.get(t, ""))  for t in trans))
				else:
					row = '[\n\t"%s",\n\t"%s",\n%s,\n\t"%s",\n\t"%s"\n]'%(j[0],j[1],",\n".join('\t"%0.3f"'%df.loc[i,a] for a in df.columns),"<br>".join(web_adress%(t,t) for t in trans),
					"<br/>".join( " ".join(vars.get(t,""))  for t in trans))
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
				trans = set([ p.transcript_id.split(":")[0] for p in  pep_to_prot[str(pep)]])
				f.write(str(pep)+"\t"+"\t".join("%.3f"%df.loc[(pep,method),a] for a in df.columns)+"\t"+"\t".join(trans)+"\n")
	
	
def to_html(out_file, csv_file_name, result, prot_to_pep,args):
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
	#######################################################################################
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
        <h2 class="etk-heading">Polymorphic Epitope Prediction Results</h2>
        <p><b>Notice:</b> The returned values represent prediction scores. They are not comparable between methods!<br/>
        Typically, prediction with higher scores mean stronger binding.</p><br/>
		  """%(" ".join(args.method),  str(args.length))
		  
	result_table="""<table id="result_table" class="display">

					<thead>
					<tr>
						<th>Peptide</th><th>Method</th>"""+"".join("<th>%s</th>"%str(a) for a in result.columns)
	if prot_to_pep:
		result_table +="<th>Origin</th><th>Influencing Variants</th>" 

	result_table+="""
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
	end_html = "</div></body></html>"
	with open(out_file, "w") as html_o:
		html_o.write(begin_html+setting+result_table+refs+end_html)
	
		

def main():
	parser = argparse.ArgumentParser(description="Polymorphic Epitope Prediction")
	parser.add_argument("-i", "--input",
						required=True,
						nargs='+',
						help="input"
						)
	parser.add_argument("-t", "--type",
						required=True,
						help="Specifies the input type")
	allele_types = parser.add_mutually_exclusive_group(required=True)
	allele_types.add_argument("-af", "--allelefile",
							 action="store_true",
							 help="Specifies the allele input as allele file.")
	allele_types.add_argument("-as", "--allelestring",
							action="store_true",
							help="Specifies the allele input as allele string.")
	parser.add_argument("-a", "--alleles",
						nargs="+",
						required=True,
						help="Allele file with frequencies (one allele and frequency per line)")
	parser.add_argument("-m", "--method",
						required=True,
						nargs="+",
						default="BIMAS",	
						help="Specifies the method used for prediction.")
	parser.add_argument("-l", "--length",
						required=False,
						type=int,
						default=9,
						help="Specifies the length of the peptides (default=9).")
	parser.add_argument("-s", "--snps",
						required=False,
						action="store_true",
						help="Filter for variations (excluding snps)")
	parser.add_argument("-in", "--indels",
						required=False,
						action="store_true",
						help="Filter for variations (excluding indels)")
	parser.add_argument("-fs", "--frame_shift",
						required=False,
						action="store_true",
						help="Filter for variations (excluding fram_shifts)")
	parser.add_argument("-o", "--output",
						required=True,
						help="Specifies the output path. Results will be written to CSV")


	#COMMENT: These options are hidden and only used for ETK2
	parser.add_argument("-og", "--output_galaxy",
						required=False,
						help=argparse.SUPPRESS)
	parser.add_argument("-od", "--outdir",
						required=False,
						default="",
						help=argparse.SUPPRESS)
	args = parser.parse_args()
	peps = []
	#print "Input ",args
	if not os.path.exists(args.outdir) and args.outdir != "":
		os.makedirs(args.outdir)

	if args.allelefile:
		alleles = FileReader.read_lines(args.alleles[0], type="Allele")
	else:
		al = map(lambda x: x.upper(), args.alleles[0].split(","))
		#print al
		alleles = [Allele( a[0]+"*"+a[1:3]+":"+a[3:] if "DR" not in a else a[:4]+"*"+a[4:6]+":"+a[6:]) for a in al]	
	
	
	#generate peptides
	genes = []
	rs_id_pos = None
	if args.type in ["refseq","swiss_accid","ensemble","swiss_gene","vcf"]:
		mart = MartsAdapter()
		if args.type in ["refseq","swiss_accid","ensemble","swiss_gene"]:
			rs_ids=[]
			rs_id_pos = 9
			tmp_annovar_input = tempfile.NamedTemporaryFile(delete=False)
			
			for pep_id in args.input:
				pep_id = pep_id.replace(",","").replace(";", "").strip()
				rows = mart.get_variant_id_from_protein_id(**{args.type:pep_id})
				#r = mart.get_transcript_information_from_protein_id(**{args.type:pep_id})
				#print pep_id
				if not rows:
					continue
					
				for row in rows:
					if row['Ensembl Gene ID'] == "error=true":
						raise Exception("BioMart query cased an error. BioMart.org might be malfunctioning.")
						
					if row['Protein location (aa)'] != "" or True:
						#print row
						if len(row['Variant Alleles'].split("/")) == 2:
							ref,obs=row['Variant Alleles'].split("/")
							if ref == "-":
								offset = 0
							elif obs == "-":
								offset = len(ref)-1
							else:
								offset = 0
						genes.append(row['Ensembl Transcript ID'])
						rs_ids.append("%s\t%s\t%i\t%s\t%s\t%s\t%s\n"%(row["Chromosome Name"],
																row["Chromosome Location (bp)"],
																int(row["Chromosome Location (bp)"])+offset,
																ref,
																obs,
																row['Ensembl Transcript ID'],
																row['Variation Name']))
			#print genes
			rs_ids = list(set(rs_ids))	
			if not rs_ids:
				raise Exception("No variants found for input or BioMart.org might not work.")
				
			tmp_annovar_input.write("".join(rs_ids))
			tmp_annovar_input.close()
		else: #if it is a vcf file convert it to annovar input
			tmp_annovar_input = tempfile.NamedTemporaryFile(delete=False)
			tmp_annovar_input.close()
			subprocess.call(ANOVA_PREP%(args.input[0],tmp_annovar_input.name), shell=True)
		
		#call anova main annotation function with prepared input file
		#print ANOVA%(tmp_annovar_input.name,tmp_annovar_input.name)
		subprocess.call(ANOVA%(tmp_annovar_input.name,tmp_annovar_input.name), shell=True)
		anovar_exonic = tmp_annovar_input.name+".exonic_variant_function"
		os.remove(tmp_annovar_input.name+".variant_function")
		os.remove(tmp_annovar_input.name)
	
		#read anovarfile
		#print "generating variants from annovar file"
		genes = set(genes)
		#print "Transcript filter ", genes

		vars = filter(lambda x: x.type != VariationType.UNKNOWN, FileReader.read_annovar_exonic(anovar_exonic,rs_id_pos=rs_id_pos))
		os.remove(anovar_exonic)
		#here we could filter the vars
		#print "Variants", vars

		if not args.snps:
			vars = filter(lambda x: x.type != VariationType.SNP, vars)
		
		if not args.indels:
			vars = filter(lambda x: x.type not in [VariationType.INS,VariationType.DEL,VariationType.FSINS,VariationType.FSDEL], vars)
			
		if not args.frame_shift:
			vars = filter(lambda x: x.type not in [VariationType.FSINS,VariationType.FSDEL], vars)
			
			
		if not vars:
			#to_html_error(args.output, "No variants could be found")
			raise Exception("No variants could be found. Please check whether your input is of human origin.")
		#for v in vars:
		#	print v.id, v, v.type
		#generate transcript from variants, proteins and peptides
		peps = Generator.generate_peptides_from_variants(vars, args.length, mart, trans_filter=genes)
		#prots = []
		#for t in trans:
	#		try:
	#			p = t.translate()
	#			prots.append(p)
	#		except ValueError as e:
	#			print e, "It will be ignored"
	#			continue
				
		#filter for peptides with varaints
		

	else:#polymorphic fasta input as file (always?)
		prot = generate_proteins_from_polymorphic_fasta(args.input[0])
	#print peps
	#print "Number of Variants", len(vars), " Heterozygious ", sum( 1 for v in vars if not v.isHomozygous), "Nof of Peps ", len(peps)
	#print "generate peptides"	
	#prots_to_pep = {}
	#for p in Generator.generate_peptides_from_protein(prots, args.length):
#		if p.vars:
#			peps.append(p)
#			for pr in p.proteins.values():
#				prots_to_pep.setdefault(str(p), []).append(pr)
	#make prediction
	#generate alleles 
	#print peps
	peps = filter(lambda x: len(x.get_all_variants())>0, peps)
	#print "prediction", len(peps)
	results = [ EpitopePredictorFactory(m).predict(peps,alleles=alleles) for m in args.method[0].split(",")]
	if not results:
		raise Exception("No prediction result was generated. Please verify that your input was correct.")
	
	r_df = results[0]
	for r in results[1:]:
		tmp_a,tmp_b = r_df.align(r, fill_value=0)
		r_df = tmp_a+tmp_b
	#print r_df
	#sys.exit()
	out_csv = args.outdir+"/snp_result.csv"
	json_out = args.outdir+"/result.json"
	pep_to_prot = {}
	for p in list(r_df.index.levels[0]):
		#print p, p.proteins
		for pr in p.proteins.values():
			pep_to_prot.setdefault(str(p), []).append(pr)
	
	if pep_to_prot:
		to_json(r_df,json_out,pep_to_prot=pep_to_prot,web_adress="<a href='http://www.ensembl.org/Homo_sapiens/Transcript/Idhistory?t=%s'>%s</a>")
		to_output_table(r_df,args.output_galaxy,pep_to_prot=pep_to_prot)
	else:
		to_json(r_df,json_out)
		to_output_table(r_df,args.output_galaxy)
	to_html(args.output, out_csv, r_df, pep_to_prot,args)

if __name__ == "__main__":
	main()
