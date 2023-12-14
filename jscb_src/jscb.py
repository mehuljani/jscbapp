
import sys 
import os
import random
from Bio import SeqIO

if len(sys.argv) == 1:
	print "\n\nError: Please provide input file in genbank format!!\n\n"
	sys.exit()

def createPic(title):
	os.system("perl cgview_xml_builder.pl -sequence JSCBinput --genes JSCB_formatted -global_label T -gene_labels T -labels_to_show label -custom featureThickness=45 labelFontSize=45 -output 1.xml")
	os.system("perl xmleditor.pl 1.xml 1")
	os.system("java -Djava.awt.headless=true -jar cgview/cgview.jar -i 1_2.xml -o %s.jpg -f jpg"%(title))
	return

def moveFiles(title):
	os.system("mkdir '%s'"%(title))
	os.system("mv '%s' '%s'"%(title+"_JSCB_output.tsv",title))
	os.system("mv '%s' '%s'"%("JSCB_output.clus",title))
	os.system("mv '%s' '%s'"%("JSCB_output.gi",title))
	os.system("mv '%s' '%s'"%("JSCB_coord",title))
	os.system("mv '%s' '%s'"%("JSCB_formatted",title))
	#os.system("mv 1.xml %s"%(title))
	#os.system("mv 1_2.xml %s"%(title))
	os.system("rm 1.xml")
	os.system("rm 1_2.xml")
	os.system("mv '%s' '%s'"%(title+".jpg",title))
	os.system("rm JSCBinput")
	#os.system("mv %s %s"%(title,"/media/Data/"))
	return

for gb_record in SeqIO.parse(sys.argv[1], "genbank"):
	f = open ("JSCBinput", "w")
	start_col, end_col, rna_start, rna_end, strand_col=[],[],[],[],[]
	actstart = []
	actend = []
	num_gene = 0
	nucs = ['A','T','C','G']
	
	accession = gb_record.id
	title = sys.argv[1].replace(".","_")+"_JSCB"
	#title = accession+"_"+gb_record.description.lower().replace(" ","_")
	print title
	if "plasmid" in gb_record.description.lower():
		continue
	seq=gb_record.seq
	seq=str(seq.upper())
	finseq=""
	for i in seq:
		if i in nucs:
			finseq+=i
		else:
			finseq+=random.choice(nucs)

	f.write(finseq)
	f.close()
	
	for feature in gb_record.features:
		start=feature.location.nofuzzy_start+1
		end=feature.location.nofuzzy_end
		if feature.strand == 1:
			strand ="+"
		else:
			strand="-"
		if feature.type == "CDS":
			if "join" not in str(feature.location):
				start_col.append(start)
				end_col.append(end)
				strand_col.append(strand)
		if feature.type == "rRNA":
			rna_start.append(start)
			rna_end.append(end)
	
	f1 = open ("JSCB_coord", "w")		
	for x, y, z in zip(strand_col,start_col,end_col):
		if (z-y+1)%3==0:
			f1.write(str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
			actstart.append(y)
			actend.append(z)
			num_gene+=1
	f1.close()			

	if num_gene == 0:
		#os.system("rm JSCBinput")
		os.system("rm JSCB_coord")
		outfile = open(title+"_JSCB_output.tsv", "w")
		outfile.write("GenBank file does not have gene co-ordinates\n")
		outfile.close()
		headers = ["seqname","source","feature","start","end","score","strand","frame"]
		headers = "\t".join(headers)+"\n"	
		formattedout = open("JSCB_formatted",'w')
		formattedout.write(headers)
		formattedout.close()
		createPic(title)
		moveFiles(title)		
		continue

	os.system("./jscb")

	gifile = open("JSCB_output.gi", "r")
	hashclussize = {}
	for line in gifile:
		try:
			hashclussize[line.split()[1]]+=1
		except:
			hashclussize[line.split()[1]]=1

	natclus = max(hashclussize, key=hashclussize.get)

	prv = -1
	gmode = 0
	gicounter=0
	gifile = open("JSCB_output.gi", "r")
	outfile = open(title+"_JSCB_output.tsv", "w")
	
	headers = ["seqname","source","feature","start","end","score","strand","frame"]
	headers = "\t".join(headers)+"\n"	
	formattedout = open("JSCB_formatted",'w')
	formattedout.write(headers)


	for line in gifile:
		if int(line.split()[1])!=int(natclus):
			if int(line.split()[0]) != prv+1:
				startgeneno = int(line.split()[0])
				startcoord = actstart[int(line.split()[0])-1]
				gmode = 1
			prv = int(line.split()[0])
		elif gmode == 1:
			endgeneno = int(line.split()[0])
			gilength = endgeneno-startgeneno
			endcoord = actend[int(line.split()[0])-2]
			if gilength >= 8:
				rnais = 0
				for x, y in zip(rna_start, rna_end):
					if startcoord < x and startcoord < y and endcoord > x and endcoord > y:
						rnais = 1
					
				if rnais == 0:
					gicounter+=1
					if gicounter == 1:
						outfile.write("Id\tStart\tEnd\n")
					outfile.write("GI-"+str(gicounter)+"\t"+str(startcoord)+"\t"+str(actend[int(line.split()[0])-2])+"\n")
					formattedline = ["gene"+str(gicounter),".","gene",str(startcoord),str(actend[int(line.split()[0])-2]),"1",".","."]
					formattedline = "\t".join(formattedline)
					formattedout.write(formattedline+"\n")
					
			gmode = 0

	if gicounter == 0:
		outfile.write("No Genomic Islands Identified\n")
	outfile.close()
	formattedout.close()
	
	createPic(title)
	moveFiles(title)


