# -*- coding: utf-8 -*-
"""
OptiType wrapper based on Fred2 and 
ETK1
"""
from __future__ import division
import numpy
import sys
import os
import argparse
import pandas
import itertools as itr
from collections import defaultdict
sys.path.append("/abi-projects/etk/Fred2/")

OPTITOPE_INPUTTYPE_UNKNOWN = 'unknown'
OPTITOPE_INPUTTYPE_EPITOPES = 'epitopes'
OPTITOPE_INPUTTYPE_EPITOPES_WITH_IMMUNOGENICITIES = 'epitopes_with_immunogenicities'
OPTITOPE_INPUTTYPE_CONSENSUS = 'consensus'
OPTITOPE_INPUTTYPE_MSA = 'msa'

from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Allele import Allele
from Fred2.Core.Result import EpitopePredictionResult

from Fred2.EpitopePrediction import EpitopePredictorFactory

from Fred2.EpitopeSelection.OptiTope import OptiTope

def extractEpitopeInformation(input, epitope_length):

	#determine input type
	epitope_info = []
	input_type = OPTITOPE_INPUTTYPE_UNKNOWN
	inp = open(input, "r")
	line2 = []
	inp2 = []

	for line in inp:
		line = str(line)
		line = line.lstrip()
		line = line.rstrip('\n')
		line = line.replace("__tc__", " ").replace("__gt__", ">")
		#print line
		line2.append(line.split("__cr____cn__"))
	#print "Input FIle:"
	for item in line2:
		#print item
		if item[0] != '':
			inp2.extend(item)

	

	for line in inp2:
		#print line
		if line.strip() == '' or line[0] == '#':
			continue

		if len(line) > 2 and line[0:2] == '>>':
			#print "optitope for MSA"
			input_type = OPTITOPE_INPUTTYPE_MSA
			info = extractEpitopeInformationFromMSA(input, epitope_length)
			#print "extracted infos ",info
			break
		
		if line[0] == '>':
			errors = "MSA has to start with '>>'."
			break
		#	input_type = OPTITOPE_INPUTTYPE_CONSENSUS
		#	info = extractEpitopeInformationFromConsensus(input, epitope_length)
		#	break

		fields = line.split()
		if len(fields) > 1 or (len(fields) == 1 and fields[0] == "Alleles:") : # epitopes and  immunogenicity
			#print "optitope from immuno"
			input_type = OPTITOPE_INPUTTYPE_EPITOPES_WITH_IMMUNOGENICITIES
			info = extractEpitopeInformationFromEpitopeAndImmunogenicities(input)
			break

		if len(fields) == 1: # epitopes
			input_type = OPTITOPE_INPUTTYPE_EPITOPES
			info = extractEpitopeInformationFromEpitopeInput(input)
			break

	#extract data
	if input_type == OPTITOPE_INPUTTYPE_UNKNOWN:
		errors = "Invalid input format."
	else:
		errors = info[0]
		epitope_info = info[1:]

	inp.close()
	return (errors, input_type, epitope_info)


def extractEpitopeInformationFromMSA(input, epitope_length):

	error = ''

	# >> antigen1
	# > seq1
	# ...
	# > seq2
	# ...
	# ...
	# >> antigen2
	# > ...
	NONE = 0
	ANTIGEN_HEADER = 1
	SEQUENCE_HEADER = 2
	SEQUENCE = 3


	last_read = NONE
	antigen = ''
	current_sequence = ''
	sequences = []
	
	consensus_info = {}
	epitope_info = []
	antigens = []
	conservation={}
	inp = open(input ,"r")
	
	line2 = []
	inp2 = []

	for line in inp:
		line = str(line)
		line = line.lstrip()
		line = line.rstrip('\n')
		line = line.replace("__tc__", " ").replace("__gt__", ">")
		line2.append(line.split("__cr____cn__"))

	for item in line2:
		if item[0] != '':
			inp2.extend(item)
	for line in inp2:
		line = line.strip()
		#if line.strip() == '' or line[0] == '#':
		if line == '' or line[0] == '#':
			continue
		#print 'in extract msa info',line
		if line[0:2] == '>>':
			if last_read == NONE or last_read == SEQUENCE: 
				if last_read == SEQUENCE: # done with previous antigen

					sequences.append(current_sequence)
					consensus = determineConsensusFromMSA(sequences)
					
					#print consensus
					if consensus[0] != '':
						consensus_info.setdefault(antigen, []).append((Protein(consensus[0].upper(), antigen, antigen), consensus[1]) )
					else: # different sequence lengths
						error = "MSA sequences of antigen " + antigen + " are of different lengths."
						break
			
				antigen = line[2:].strip()
				last_read = ANTIGEN_HEADER

				current_sequence = ''
				sequences = []
				
			else: # wrong format
				error = "Input format does not comply with required MSA format."
				break

		elif line[0] == '>':
			if last_read == ANTIGEN_HEADER or last_read == SEQUENCE:
				if last_read == SEQUENCE:
					sequences.append(current_sequence)
					current_sequence = ''
					
				last_read = SEQUENCE_HEADER

		else: # sequence

			if not isValidMSASequence(line):
				error = "Invalid amino acid sequence given for antigen " + antigen + ". "
				break
			
			if last_read == SEQUENCE_HEADER or last_read == SEQUENCE:
				#current_sequence += line[:-1]
				current_sequence += line
				last_read = SEQUENCE
			else:
				error = "Input format does not comply with required MSA format."
				break

	if error == '':
		if last_read == SEQUENCE:
			sequences.append(current_sequence)
			consensus = determineConsensusFromMSA(sequences)
			if consensus[0] != '':
				consensus_info.setdefault(antigen, []).append((Protein(consensus[0].upper(), antigen, antigen), consensus[1]) )
			else: # different sequence lengths
				error = "MSA sequences of antigen " + antigen + " are of different lengths."
		else:
			error = "Input format does not comply with required MSA format."

		if error == '':
			(error,  conservation) = extractEpitopesAndConservationFromConsensus(consensus_info, epitope_length)
			#print "infos form extractEpitopesAndConservationFromConsensus", conservation
			
	inp.close()			
	return (error, conservation)


# this function is not used at the moment - it might need some more work
def extractEpitopeInformationFromConsensus(input, epitope_length):

	error = ''

	# > antigen1
	# ...
	# ...
	# > antigen2
	# ...
	NONE = 0
	ANTIGEN_HEADER = 1
	SEQUENCE = 2

	last_read = NONE
	antigen = ''
	current_sequence = ''
	
	consensus_info = []

	antigens = []
	epitope_info = []

	inp = open(input, "r")
	line2=[]
	inp2=[]
	for line in inp:
		line = str(line)
		line = line.lstrip()
		line = line.rstrip('\n')
		line = line.replace("__tc__", " ")
		line2.append(line.split("__cr____cn__"))

	for item in line2:
		if item[0] != '':
			inp2.extend(item)

	for line in inp2:
		line = line.strip()
		#if line.strip() == '' or line[0] == '#':
		if line == '' or line[0] == '#':
			continue
		
		if line[0] == '>':
			if last_read == NONE or last_read == SEQUENCE: 
				if last_read == SEQUENCE: # done with previous antigen
					consensus = separateFrequenciesFromSequence(current_sequence)
					if consensus[0] == '': # wrong format
						error = 'Input format does not comply with required consensus format.'
						break
					else:
						consensus_info.append((antigen, consensus))
			
				antigen = line[2:].strip()
				last_read = ANTIGEN_HEADER
			else: # wrong format
					error = "Input format does not comply with required consensus format."
					break
		else: # sequence
			if last_read == ANTIGEN_HEADER or last_read == SEQUENCE:
				#current_sequence += line[:-1]
				current_sequence += line
				last_read = SEQUENCE
			else:
				error =  'Input format does not comply with required consensus format.'
				break

	if error == '':
		if last_read == SEQUENCE:
			consensus = separateFrequenciesFromSequence(current_sequence)
			if consensus[0] == '': # wrong format
				error =  'Input format does not comply with required consensus format.'
			else:
				consensus_info.append((antigen, consensus))
		else:
			error =  'Input format does not comply with required consensus format.'

		if error == '':
			#print "consensus_info ", consensus_info
			(error, conservation) = extractEpitopesAndConservationFromConsensus(consensus_info, epitope_length)
		
	
	inp.close()
	return (error, conservation)

#TODO: This has be modified to return an EpitopePredictionResult
def extractEpitopeInformationFromEpitopeAndImmunogenicities(input):
	# alleles a1 a2 ... 
  # epitope immunogenicities

	error = ''
	epitope_info = {}
	immunogenicities = {}
	alleles = []

	number_of_fields = -1
	number_of_alleles = 0

	read_alleles = 0
	inp = open(input, "r")
	line2=[]
	inp2=[]
	for line in inp:
		line = str(line)
		line = line.lstrip()
		line = line.rstrip('\n')
		line = line.replace("__tc__", " ")
		line2.append(line.split("__cr____cn__"))

	for item in line2:
		if item[0] != '':
			inp2.extend(item)
	
	for line in inp2:
		line = line.strip()
		#print "In Epitope from table",line
		if line == '' or line[0] == '#':
			continue
			
		fields = line.split()
		#print "Fields ",fields
		number_of_fields = len(fields) 
		if not read_alleles:
			if fields[0].strip() == "Alleles:":
				read_alleles = 1
				number_of_alleles = number_of_fields - 1
				if number_of_alleles > 0:
					alleles = map(lambda x: Allele(x) , fields[1:])
				else:
					error = "No alleles given."
					print error
					break
			else:
				error = "When entering a table of epitopes with immunogenicities the first line has to start with 'Alleles:' ."
				print error
				break
		else:
			
			if number_of_fields > number_of_alleles + 1:
				#print "In Origin part",fields,number_of_fields,number_of_alleles
				origin = {fields[j]:Protein(fields[0].upper(),fields[j],fields[j]) for j in xrange(number_of_alleles+1,number_of_fields)}
				epitope = Peptide(fields[0].upper(), proteins=origin)
				#print epitope
				if isValidAASequence(str(epitope)):
					conservation = 1.
					epitope_info[epitope] = conservation
				else:
					error = "Invalid amino acid sequence " + epitope + "." 
					break
				try:
					for i in xrange(number_of_alleles):
						#print alleles[i],float(fields[i+1])
						if alleles[i] not in immunogenicities:
							immunogenicities[alleles[i]] = {}
						immunogenicities[alleles[i]][epitope] = float(fields[i+1])
				except ValueError:
					error = "Invalid immunogenicity value for epitope " + epitope + " and allele " + alleles[i] + ": " + fields[i+1] 
					break
			elif number_of_fields < number_of_alleles:
				error = "Not enough immunogenicity entries."
				break
			else:
				epitope = Peptide(fields[0].upper())
				if isValidAASequence(str(epitope)):
					conservation = 1.
					epitope_info[epitope] = conservation
				else:
					error = "Invalid amino acid sequence " + epitope + "." 
					break
				
#				if not immunogenicities.has_key(epitope):
#					immunogenicities[epitope] = {}
#				else: # epitope appears twice
#					error = "Epitope " + epitope + " appears more than once." 
#					break
				
				try:
					for i in xrange(number_of_fields-1):
						if alleles[i] not in immunogenicities:
							immunogenicities[alleles[i]] = {}
						immunogenicities[alleles[i]][epitope] = float(fields[i+1])
				except ValueError:
					error = "Invalid immunogenicity value for epitope " + epitope + " and allele " + alleles[i] + ": " + fields[i+1] 
					break
	#print immunogenicities
	#print "Immunogenicity",immunogenicities
	immu = EpitopePredictionResult.from_dict(immunogenicities)
	immu.index = pandas.MultiIndex.from_tuples([tuple((i,"custom")) for i in immu.index], names=['Seq','Method'])
	inp.close()
	return (error, epitope_info, immu)


def extractEpitopeInformationFromEpitopeInput(input):
	# epitope

	error = ""
	epitope_info = {}
	epitope_length = -1
	inp=open(input, "r")

	line2 = []
	inp2 = []
	for line in inp:
		line = str(line)
		line = line.lstrip()
		line = line.lstrip('\n')
		line2.append(line.split("__cr____cn__"))
	for item in line2:
		if item[0] != '':
			inp2.extend(item)
	
	for line in inp2:
		line = line.strip()
		if line == '' or line[0] == '#':
				continue
		
		fields = line.split()
		if len(fields) != 1:
			error = "Input format does not comply with required format of epitope list. Please enter one epitope per line."
			break

		epitope = fields[0].upper()
		if isValidAASequence(epitope):
				if epitope_length == -1:
					epitope_length = len(epitope)
				else:
					if len(epitope) != epitope_length:
						error = "Epitopes have to be of equal length."
						break

				conservation = 1.
				epitope_info[Peptide(epitope)] = conservation
			 
		else:
				error = "Invalid amino acid sequence " + epitope + "." 
				break

	inp.close()
	return (error, epitope_info)

def separateFrequenciesFromSequence(sequence):
	consensus_sequence = ''
	frequencies = []
	error = ''

	open_bracket = 0
	current_aa = ''
	last_s = ''
	freq_string = ''

	for s in sequence:
		if s == '[':
			if open_bracket == 1 or last_s == '' or last_s == ']':
				error = 'Invalid consensus format.'
				break
			else:
				open_bracket = 1
				freq_string = ''
				current_aa = last_s
		elif s == ']':
			if open_bracket == 0:
				error = 'Invalid consensus format.'
				break
			else:
				open_bracket = 0
				consensus_sequence += current_aa
				freq = 1.
				try:
					freq = float(freq_string)
				except ValueError:
					error = "Invalid frequency string."
					break

				if freq <= 0 or freq > 1.0:
					error = "Invalid frequency."
					break
				else:
					frequencies.append(freq)
		else:
			if open_bracket == 0 and s != '-': #ignore gaps
				consensus_sequence += s
				frequencies.append[1.]
			else:
				freq_string += s

		last_s = s

	if error != '':
		consensus_sequence = ''
		frequencies = []

	return (consensus_sequence, numpy.array(frequencies))


def determineConsensusFromMSA(sequences):
	
	consensus_sequence = ''
	frequencies = numpy.array([])
	error = ''
	aa_tracker = []

	#print sequences
	number_of_sequences = len(sequences)
	if number_of_sequences > 0:
		seq_len = len(sequences[0])
		frequencies = numpy.zeros(seq_len)

		for i in xrange(seq_len):
			aa_tracker.append({})

		for seq in sequences:
			if len(seq) == seq_len:
					
				for p in xrange(seq_len):
					c = seq[p]
					if aa_tracker[p].has_key(c):
						aa_tracker[p][c] += 1
					else:
						aa_tracker[p][c] = 1

			else: # sequences of different lengths
				error = "MSA sequences are of different lengths."
				break

		if error == '':
				for i in range(len(aa_tracker)):
					max = 0
					max_aa = ''
					for aa, count in aa_tracker[i].items():
							if count > max:
								max = count
								max_aa = aa
					
					if max_aa != '-': #don't include gaps
						consensus_sequence += max_aa
						frequencies[i] = 1. * max/number_of_sequences

	return consensus_sequence, frequencies

def extractEpitopesAndConservationFromConsensus(consensus_info, epitope_length):

	#print "consensus_info", consensus_info
	error = ''
	epitope_data = []
	epitope_antigen_data = {}
	antigens = consensus_info.keys()
	conservation = {}
	epitopes = {}
	
	#print "type of consensus_info ",type(consensus_info)
	for ci in consensus_info.values():
		consensus = str(ci[0][0])
		frequencies = ci[0][1]

		#consensus = consensus.upper()    
		#if isValidAASequence(consensus):

		for i in xrange(len(consensus) - epitope_length + 1):
			epitope = Peptide(consensus[i:i+epitope_length], proteins={ci[0][0].transcript_id: ci[0][0]})
			#print "peptide", epitope
			co = numpy.product(frequencies[i:i+epitope_length])
			
			if not conservation.has_key(epitope):
				#print "test type prot in epitopes", epitope.proteins
				conservation[epitope] = co
				epitopes[epitope] = epitope
			else:
				epitope = epitopes[epitope]
				#print "test type prot in epitopes", epitope.proteins
				if conservation[epitope] < co:
					conservation[epitope] = co
				
				epitopes[epitope].proteins.update({ci[0][0].transcript_id: ci[0][0]})
					
				
		
			#epitope_data.append((epitope, conservation))
 
	#if error == '':
	#	for k, v in epitope_antigen_data.items():
	#		epitope_data.append((k, v[0], v[1]))

	return (error, conservation)

def isValidAASequence(epitope):

	epitope = epitope.strip()
	if epitope == '':
		return 1

	if not epitope.isalpha():
		return 0

	# check for invalid letters
	for aa in ['B','J','O','U','X','Z']:
		if aa in epitope:
			return 0

	return 1

def isValidMSASequence(seq):
	seq = seq.upper()

	return isValidAASequence(seq.replace('-',''))


def generate_alleles(allele_file, generated=None):
	"""
		generate allele objects from input
	"""
	result=[]
	with open(allele_file, "r") as f:

		for l in f:
			if l != '' and l != '\n' and l != ' ':
				li = l.lstrip()
				li = li.rstrip('\n')
				allele_pairs = li.split(',')
				if allele_pairs[0] == '':
					continue
				#print allele_pairs
				for a_p in allele_pairs:
					al_name,freq = a_p.split()
					if generated is None:
						result.append(Allele(al_name, prob=float(freq)))
					else:
						if al_name in generated:
							a = generated[al_name]
							a.prob = float(freq)
							result.append(a)
	#print "Results after file reading, ",result
	return result


def to_html(out_file, result, instance, pred_method):
	#print "In to html funktion"
	"""
		generates galaxy html output
		
		:param: string out_file: the output file
		:param: list(Peptide) result: the list of selected Peptides
		:param: Coopr.Solver instance: the coopre instance with the loaded results
	"""

	begin_html = """<?xml version="1.0" encoding="utf-8" ?>
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
	<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
	<head>
		<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
		<link rel="stylesheet" href="/static/style/etk.css" type="text/css" />
		<link rel="stylesheet" href="/static/style/datatable/css/jquery.dataTables.css" type="text/css" />
		<link rel="stylesheet" href="/static/style/datatable/dataTables.fixedColumns.css" type="text/css" />
		<link rel="stylesheet" href="/static/style/datatable/css/tabletools/dataTables.tableTools.css" type="text/css" />
		<link rel="stylesheet" href="/static/style/datatable/css/dataTables.colVis.css" type="text/css" />
		<script type="text/javascript" src="/static/scripts/libs/jquery/jquery.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/jquery-ui.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/datatable/media/js/jquery.dataTables.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/datatable/media/js/fnFakeRowspan.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/Scroller/js/dataTables.scroller.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/FixedColumns/js/dataTables.fixedColumns.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/TableTools/js/dataTables.tableTools.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/jquery/datatable/extensions/ColVis/js/dataTables.colVis.js"></script>
		<script type="text/javascript" src="/static/scripts/libs/optitope.js"></script>
	</head>

	<body>
		<div class="document">
	"""
	#######################################################################################
	setting = """<h2 class="etk-heading">Input Summary</h2>
				 	<table class="etk-parameterT">	
						<tr><th>Number of candidate epitopes:</th><td>%i</td></tr>
						<tr><th>Number of target alleles:</th><td>%i</td></tr>
						<tr><th>Prediction Method:</th><td>%s</td></tr>
					</table>	
			"""%(len(instance.E),len(instance.A),pred_method)
	#print setting
	
	#######################################################################################
	cons = ["<tr><th>Maximum number of epitopes to select </th><td>= "+ str(int(instance.k))+"</td></tr>"]
	if float(instance.t_c) > 0:
		cons.append("<tr><th>Epitope conservation </th><td>&ge; "+ str(float(instance.t_c)*100)+"%</td></tr>")
		
	if float(instance.t_allele) > 0:
		cons.append("<tr><th>Covered alleles </th><td>&ge; "+ str(int(instance.t_allele))+"</td></tr>")
		
	if float(instance.t_var) > 0:
		cons.append("<tr><th>Covered antigens</th><td>&ge; "+ str(int(instance.t_var))+"</td></tr>")
		
	constraints = """<h2 class="etk-heading">Constraints</h2>
				 	<table class="etk-parameterT">	
						%s
					</table>	
			"""%"\n".join(cons)
	#print constraints
	#######################################################################################
	res = ["<tr><th>Selected epitopes:</th><td>"+ str(len(result))+"</td></tr>"]
	
	#print instance.A
	if int(instance.t_var) > 0:
		cov_anti=[]
		for an in instance.Q:
			for e in result:
				if e in instance.E_var[an]:
					cov_anti.append(an)
		cov_anti = set(cov_anti)
		res.append("<tr><th>Covered antigens:</th><td>"+ str(len(cov_anti))+" of "+str(len(instance.Q))+"</td></tr>")
	cov_als = []
	res_set = set(result)
	locus = {}
	for a in instance.A:
		eps_of_all_i = list(instance.A_I[a])
		if res_set.intersection(set(eps_of_all_i)): 
			cov_als.append(a)
			locus.setdefault(str(a).split("*")[0],set()).add(a)
	cov_als = list(set(cov_als))
	res.append("<tr><th>Covered alleles:</th><td>"+ str(len(cov_als))+" of "+str(len(instance.A))+" </td></tr>")
	res.append("<tr><th>Locus coverage:</th></tr>")
	#print "Covered alleles", cov_als
	pop_cov= 1
	for k,g in sorted(locus.iteritems()):
		locus = list(g)
		#print "Locus, ",k,type(k)
		#print "pop parts,",[float(instance.p[a]) for a in locus ], (1 - sum( float(instance.p[a]) for a in locus))**2
		pop_cov *=(1.0 - sum( float(instance.p[a]) for a in locus))**2
		#print "Population Coverage locus",k,pop_cov 
		covered = 1.0-(1.0 - sum( float(instance.p[a]) for a in locus))**2
		res.append("<tr><th></th><td>%s</td><td>%.2f</td></tr>"%(k,covered*100))
		
	res.append("<tr><th>Population converage:</th><td></td><td>%.2f</td></tr>"%((1.0-pop_cov)*100.0))
	results_scalar = """<h2 class="etk-heading">Results</h2>
			<table class="etk-parameterT">
			%s
			</table>
	 	   """%"\n".join(res)
	#print results_scalar
	#######################################################################################
	
	table ="""<br/><div>
			<table id="result_table_optitope" class="display">
				<thead><tr>%s</tr></thead>
				<tbody>%s</tbody>
			</table><p style='float:left;'><a href='optitope_result.csv'>Download as CSV</a></p></div><br/>
	"""
	is_antigen_cons = int(instance.t_var) > 0
	header = "<tr><th>Epitope</th><th class='detail'>Conservation</th><th>Fraction of overall immunogenicity</th><th>Covered alleles</th><th class='detail'>Allele-wise immunogenicity</th>%s</tr>"%("<th>Covered antigens</th>"if is_antigen_cons else "" )
	
	rows = []
	overall_imm = sum( float(instance.i[e,a])*float(instance.p[a]) for e in result for a in instance.A)
	for e in result:
		row = "<tr><td>"+str(e)+"</td>"
		if float(instance.t_c) > 0:
			row += "<td class='detail'>"+str(float(instance.c[e])*100)+"%</td>"
		else:
			row += "<td class='detail'>100%</td>"
		row +="<td>%0.2f</td>"%(sum(float(instance.i[e,a])*float(instance.p[a]) for a in instance.A)/overall_imm)
		row +="<td>%s</td>"%" ".join( "%s"%str(a)  for a in instance.A if e in instance.A_I[a])
		row += "<td class='detail'>%s</td>"%"<br/>".join( "%s %.3f"%(str(a), float(instance.i[e,a])*float(instance.p[a]))  for a in instance.A if e in instance.A_I[a] ) 
		if is_antigen_cons:
			row += "<td>%s</td>"%" ".join( str(q) for q in instance.Q if e in instance.E_var[q])
		row+="</tr>"
		rows.append(row)
	table = table%(header,"\n".join(rows))
		
	
	refs = """
    <h2 class="etk-heading">References</h2>
    <ol>
<li>Toussaint N.C., Doennes P., Kohlbacher O. (2008) A Mathematical Framework for the Selection of an Optimal Set of Peptides for Epitope-Based Vaccines. PLoS Comp Biol 4:e1000246</li>
<li>Toussaint N. C. and Kohlbacher O. (2009) OptiTope - a web server for the selection of an optimal set of peptides for epitope-based vaccines. Nucleic Acids Research, 37:W617-22. doi:10.1093/nar/gkp293</li>
<li>Rammensee H. , Bachmann J. , Emmerich N.P. , Bachor O.A. , Stevanovic S. (1999) SYFPEITHI: database for MHC ligands and peptide motifs. Immunogenetics 50:213–219.</li>
<li>Parker K.C. , Bednarek M.A. , Coligan J.E. (1994) Scheme for ranking potential HLA-A2 binding peptides based on independent binding of individual peptide side-chains. J Immunol 152:163-175.</li>
<li>Doennes P., Kohlbacher O. (2006) SVMHC: a server for prediction of MHC-binding peptides. Nucleic Acids Res 34:W194-W197.</li>
<li>Lundegaard C, Lamberth K, Harndahl M, Buus S, Lund O, Nielsen M. (2008) NetMHC-3.0: accurate web accessible predictions of human, mouse and monkey MHC class I affinities for peptides of length 8-11. Nucleic Acids Res. 1;36(Web Server issue):W509-12.</li>
<li>Nielsen M, Lundegaard C, Blicher T, Lamberth K, Harndahl M, et al. (2007) NetMHCpan, a Method for Quantitative Predictions of Peptide Binding to Any HLA-A and -B Locus Protein of Known Sequence. PLoS ONE 2(8): e796. doi: 10.1371/journal.pone.0000796</li>
<li>Nielsen, M. and Lund, O. (2009) NN-align. An artificial neural network-based alignment algorithm for MHC class II peptide binding prediction. BMC Bioinformatics, 10, 296.</li>
<li>Karosiene E, Rasmussen M, Blicher T, Lund O, Buus S, and Nielsen M. (2013) NetMHCIIpan-3.0, a common pan-specific MHC class II prediction method including all three human MHC class II isotypes, HLA-DR, HLA-DP and HLA-DQ. Immunogenetics.</li>
<li>Zhang L, Chen Y, Wong H-S, Zhou S, Mamitsuka H, et al. (2012) TEPITOPEpan: Extending TEPITOPE for Peptide Binding Prediction Covering over 700 HLA-DR Molecules. PLoS ONE 7(2): e30483. doi: 10.1371/journal.pone.0030483</li>
    </ol>
	"""
	
	#######################################################################################
	end_html = """</div>
	</body></html>"""
	with open(out_file, "w") as html_o:
		html_o.write(begin_html+setting+constraints+results_scalar+table+refs+end_html)


def to_csv(out_file, result, instance, pred_method):
	"""
		Writes model to CSV
	"""
	with open(out_file, "w") as f:
		f.write("Prediction method: "+pred_method+"\n\n")
		
		cons = ["Maximum number of epitopes to select = "+ str(int(instance.k))+"\n"]
		if float(instance.t_c) > 0:
			cons.append("Epitope conservation =>"+ str(float(instance.t_c)*100)+"%\n")
		
		if float(instance.t_allele) > 0:
			cons.append("Covered alleles >= "+ str(int(instance.t_allele))+"\n")
		
		if float(instance.t_var) > 0:
			cons.append("Covered antigens >= "+ str(int(instance.t_var))+"\n")
		f.write("CONSTRAINTS\n"+"".join(cons)+"\n")
		
		res = ["Selected epitopes,"+ str(len(result))+""]
	
		if int(instance.t_var) > 0:
			cov_anti=[]
			for an in instance.Q:
				for e in result:
					if e in instance.E_var[an]:
						cov_anti.append(an)
			cov_anti = set(cov_anti)
			res.append("Covered antigens,"+ str(len(cov_anti))+" of "+str(len(instance.Q))+"")
		cov_als = []
		res_set = set(result)
		locus = {}
		for a in instance.A:
			eps_of_all_i = list(instance.A_I[a])
			if res_set.intersection(set(eps_of_all_i)): 
				cov_als.append(a)
				locus.setdefault(str(a).split("*")[0],set()).add(a)
		cov_als = set(cov_als)
		res.append("Covered alleles,"+ str(len(cov_als))+" of "+str(len(instance.A))+"")
		res.append("Locus coverage:\n")
		
		pop_cov= 1
		for k,g in locus.iteritems():
			locus = list(g)
			pop_cov *=(1.0 - sum( float(instance.p[a]) for a in locus))**2 
			covered = len(locus)/sum(1 for a in instance.A if a.split("*")[0] == k)
			res.append(",%s,%.2f"%(k,covered*100))
		res.append("Population converage:,,%.2f"%((1.0-pop_cov)*100))
		f.write("RESULTS\n"+"\n".join(res)+"\n\n")
		
		is_antigen_cons = int(instance.t_var) > 0
		header = "Epitope,Conservation,Fraction of overall immunogenicity,Covered alleles%s\n"%(",Covered antigens"if is_antigen_cons else "" )
	
		rows = []
		overall_imm = sum( float(instance.i[e,a])*float(instance.p[a]) for e in result for a in instance.A)
		for e in result:
			row = str(e)+","
			if float(instance.t_c) > 0:
				row += str(float(instance.c[e])*100)+","
			else:
				row += "100%,"
			row +="%0.2f,"%(sum(float(instance.i[e,a])*float(instance.p[a]) for a in instance.A)/overall_imm)
			row +="%s"%" ".join( str(a)  for a in instance.A if e in instance.A_I[a])
			#row += "%s"%"\n".join( str(float(instance.i[e,a])*float(instance.p[a]))  for a in instance.A if e in instance.A_I[a] )
			if is_antigen_cons:
				row += ",%s"%" ".join( str(q) for q in instance.Q if e in instance.E_var[q])
			#row+="\n"
			rows.append(row)
		f.write(header+"\n".join(rows)+"\n\n")
		f.write("TARGET POPULATION/INDIVIDUAL\n")
		f.write(",Allele,Probability\n")
		for a in instance.A:
			f.write(",%s,%.5f\n"%(str(a),float(instance.p[a])))
			
		
def main():
	'''
		some input stuff
	'''
	parser = argparse.ArgumentParser(description="Epitope Selection for vaccine design.")
	parser.add_argument("-i", "--input",
						required=True,
						help="Special MSA file, peptide list, or peptide with immunogenicity file"
						)
	parser.add_argument("-a", "--alleles",
						required=True,
						help="Allele file with frequencies (one allele and frequency per line)")
	parser.add_argument("-m", "--method",
						required=True,
						default="BIMAS",	
						help="Specifies the method used for prediction.")
	parser.add_argument("-l", "--length",
						required=False,
						type=int,
						default=9,
						help="Specifies the length of the peptides (default=9).")
	parser.add_argument("-k", "--k",
						required=False,
						type=int,
						default=10,
						help="Specifies the number of epitopes to select")
	parser.add_argument("-t", "--threshold",
						required=True,
						type=float,
						default=0,
						help="Specifies the binding threshold for all alleles")
	parser.add_argument("-o", "--output",
						required=True,
						help="Specifies the output path. Results will be written to CSV")
	parser.add_argument("-io", "--internal_output",
						required=True,
						help="Specifies the output path to internal output. Results will be written to list of peptides")
	parser.add_argument("-c_al", "--cons_allele",
						required=False,
						type=float,
						help="Acitvates allele converage constraint with specified threhsold")
	parser.add_argument("-c_ant", "--cons_antigen",
						required=False,
						type=float,
						help="Acitvates antigen converage constraint with specified threhsold")
	parser.add_argument("-c_con", "--cons_conservation",
						required=False,
						type=float,
						help="Acitvates conservation constraint with specified threhsold")
	parser.add_argument("-am", "--available",
						required=False,
						action="store_true",
						help="Returns all available methods and their allele models.")
	parser.add_argument("-cons", "--conservation",
						default="",
						required=False,
						help="Specifies a Conservation file. First colum is the peitope seq second column the conservation.")
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
	
	#print "DEBUG", args
	
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)	
	
	epitopePrediciton  = None
	alleles = None
	conservation = None

	if args.input == "test":
		args.input = "/abi-projects/etk/galaxy-etk/tool-data/hcv_msa_test_data.mmsa.txt"

	error, input_type, epi_infos = extractEpitopeInformation(args.input, args.length)
	#print "input_type", input_type
	if error != '':
		raise Exception(error)
		sys.exit(-1)
	
	if input_type == OPTITOPE_INPUTTYPE_UNKNOWN:
		Exception(error)
		sys.exit(-1)
		
	elif input_type in [OPTITOPE_INPUTTYPE_EPITOPES, OPTITOPE_INPUTTYPE_MSA]:
		# get alleles 
		alleles = generate_alleles(args.alleles)
		
		method = args.method
		conservation = epi_infos
		epitopePrediciton = EpitopePredictorFactory(method).predict(epi_infos[0].keys(), alleles)
		#print "predictions",
		#print epitopePrediciton
		
	elif input_type == OPTITOPE_INPUTTYPE_EPITOPES_WITH_IMMUNOGENICITIES:
		conservation = epi_infos[0]
		epitopePrediciton = epi_infos[1]
		alleles = generate_alleles(args.alleles, generated={a.name:a for a in epitopePrediciton.columns.values.tolist()})
	else:
		Exception(error)
		sys.exit(-1)
	
	#print "Alleles ", alleles
	#print epitopePrediciton
	thresh = {a.name:args.threshold for a in epitopePrediciton.columns}
	#print "outside ", thresh
	opti = OptiTope(epitopePrediciton, threshold=thresh, k=args.k, solver="cplex", threads=1, verbosity=0)
	
	#set constraints
	if args.cons_allele > 0.0:
		#print "allele constraiont enforced"
		opti.activate_allele_coverage_const(args.cons_allele/100)
		
	if args.cons_antigen > 0.0:
		opti.activate_antigen_coverage_const(args.cons_antigen/100)
		
	if args.cons_conservation > 0.0:
		if input_type == OPTITOPE_INPUTTYPE_MSA:
			opti.activate_epitope_conservation_const(args.cons_conservation/100, conservation=conservation)
		elif args.conservation != "":
			conservatio = {}
			with open(args.conservation, "r") as f:
				for l in f:
					seq,cons = l.replace(",", " ").replace(";", " ").split()
					conservation[seq.strip().upper()] = float(cons.strip())
			opti.activate_epitope_conservation_const(args.cons_conservation/100, conservation=conservation)
	try:
		result = opti.solve()
		#print result
		to_csv(args.outdir+"/optitope_result.csv",result, opti.instance, args.method)
		to_html(args.output,result, opti.instance, args.method)
		with open(args.internal_output, "w") as o:
			o.write("\n".join(str(e) for e in result))
	except ValueError as e:
		ValueError("Could not optimally solve the problem. Please modify your constraints.")
		sys.exit(-1)
	except Exception as e:
		raise Exception(e)
		sys.exit(-1) 
	
	
	
if __name__ == "__main__":
	main()
