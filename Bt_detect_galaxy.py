#!/global/conda/envs/artwork/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys
import argparse

def get_parser():
	"""
	Parse arguments
	@return: arguments list
	@rtype: parser object
	"""
	parser = argparse.ArgumentParser(description='Bacillus Thuringiensis serovars aizaiwai and kurstaki detection')
						
	parser.add_argument('-b', action="store", dest='blastx_output',
						type=str, required=True, help='blastx output tabular file')		

	parser.add_argument('-o', action="store", dest='output',
						type=str, default="output.tsv", help='output file (default:output.tsv')	

	parser.add_argument('-t', action="store", dest='decision_table',
						type=str, default="table_bt.txt", help='Decision table (default:table_bt.txt)')
						
	parser.add_argument('-min_id', action="store", dest='min_identity',
						type=str, default='90', help='minimum percent of blast identity (default:90)')					

	parser.add_argument('-min_cov', action="store", dest='min_coverage',
						type=str, default="90", help='minimum percent of blast coverage (default:90)')	

	return parser


		
		
def decision_table_parser(file) :	

	f = open(file, 'r')
	lines = f.readlines()
	f.close()
	
	flag = False
	dico_marker = {}
	header = True
	header_list = []
	dico_association = {}
	result_key = {}
	i = 0
	
	for line in lines :
		if '###' in line:
			flag = True
			continue
		if not flag :
			line = line.rstrip().split('\t')
			dico_marker[line[0]] = line[1:]
		if flag and header :
			header = False
			header_list = line.rstrip().split('\t')[1:]
		elif flag :
			line = line.rstrip().split('\t')
			dico_association[i] = line[1:]
			if line[0] not in result_key :
				result_key[line[0]] = []
			result_key[line[0]].append(i)
			i+=1
			
	return 	dico_marker, header_list, dico_association,result_key	
		
		
def seqLen(fasta):

	file = open(fasta,'r')
	lines = file.readlines()
	file.close()
	
	dico_sequence = {}
	for line in lines :
		line = line.rstrip()
		id = line.split('\t')[1]
		len = int(line.split('\t')[23])
		if id not in dico_sequence :
			dico_sequence[id] = len
	
	return dico_sequence
	
	
def blast_parser(blast_file, mincov, minid, dico_len):

	file = open(blast_file,'r')
	lines = file.readlines()
	file.close()
	
	marker_list = []
	
	for line in lines :
		line = line.rstrip().split('\t')
		if float(line[2]) < float(minid):
			continue
		cov = float(line[3])/float(dico_len[line[1]])*100
		if float(cov) < float(mincov):
			continue
		if line[1] not in marker_list :
			marker_list.append(line[1])
	
	return marker_list
	
	
def assign(marker_found,dico_marker,header_list,dico_association,genome_id,result_key):
	
	group_list = []
	result = genome_id
	for group in dico_marker :
		for marker in marker_found :
			if marker in dico_marker[group]:
				group_list.append(group)
				result = f'{result}\t{group}(+)'
				break
		if (group + '(') not in result :
			result = f'{result}\t{group}(-)'
					

	for Bacillus in dico_association :
		i = 0
		flag = True
		for element in dico_association[Bacillus]:
			if element == '1' and header_list[i] not in group_list :
				flag = False
				break
			elif element == '-1' and header_list[i] in group_list :
				flag = False
				break
			i+=1
		if flag :
			for element in result_key :
				if Bacillus in result_key[element]:
					result = f'{result}\t{element}'
					break
	
	return result
	
	
#main function	
def main():	

	
	# Get arguments 
	parser=get_parser()
	
	# Print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	Arguments=parser.parse_args()	
	
	# Process decision table
	dico_marker, header_list, dico_association,result_key = decision_table_parser(Arguments.decision_table)
	dico_len = seqLen(Arguments.blastx_output)
	
	output_file = open(Arguments.output,'w')
	genome_id = "sample"
	blast_out = Arguments.blastx_output
	marker_list = blast_parser(blast_out,Arguments.min_coverage,Arguments.min_identity, dico_len)
	result = assign(marker_list,dico_marker,header_list,dico_association,genome_id,result_key)
	output_file.write(result+'\n')
	output_file.close()	
		
		
if __name__ == "__main__":
	main()	   	
	
