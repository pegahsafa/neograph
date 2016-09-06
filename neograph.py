import pysam
import sys
from optparse import OptionParser

class Gene:
	'''corresponds to a gene ID from the reference annotation
	All the core functionality is called within the 
	getPaths method'''

	def __init__(self, geneID, chr, chrStart, chrEnd, strand, 
		reads, splicedReads, coverageDict, exons, 
		terminalExons, stopCodonStarts, TSSs, startCodonStarts, startCodonTSSdict):
		self.ID = geneID
		self.chr = chr
		self.chrStart = int(chrStart)
		self.chrEnd = int(chrEnd)
		self.strand = strand
		self.reads = reads
		self.splicedReads = splicedReads
		self.coverageDict = coverageDict
		self.canonicalExons = exons
		self.canonicalTerminalExons = terminalExons
		self.TSSs = TSSs
		if len(startCodonStarts) > 0:
			self.stopCodonStart = self.getMostDownstream(list(stopCodonStarts))
			self.startCodonStart = self.getMostUpstream(list(startCodonStarts))
			self.startCodonStarts = startCodonStarts
			self.startCodonTSSdict = startCodonTSSdict
			self.isCoding = True
		else:
			self.isCoding = False

	def removeUncoveredCanonical(self,percentCovMin=72):
		'''Remove exons from ref annotation that are not at least 72% (default) covered
		with bases with the minimum coverage set by the coverage dictionary (10 reads by default)'''
		canonicalExons = list(self.canonicalExons)
		for exon in canonicalExons:
			percentCov,medianCov = self.getPercentCovered(exon[0],exon[1])
			if percentCov < percentCovMin:
				self.canonicalExons.remove(exon)

	def defineIntrons(self):
		'''identify intragenic regions that are do not overlap with
		annotated exons'''
		exonLst = list(self.canonicalExons)
		exonLst.sort(key=lambda x: x[0])
		self.introns = set([])
		i = 0
		while i < (len(exonLst)-1):
			if (exonLst[i+1][0] - exonLst[i][1]) > 0:
				intron = (exonLst[i][1]+1,exonLst[i+1][0]-1)
				self.introns.add(intron)
			i += 1

	def getMostUpstream(self,lst):
		lst.sort()
		if self.strand == '+':
			return lst[0]
		elif self.strand == '-':
			return lst[-1]

	def getMostDownstream(self, lst):
		lst.sort()
		if self.strand == '+':
			return lst[-1]
		elif self.strand == '-':
			return lst[0]

	def filterCoverageDict(self, minCov=10):
		'''remove bases from coverage dictionary
		that fall below 10 (default) reads'''
		positions = self.coverageDict.keys()
		for pos in positions:
			if self.coverageDict[pos] < minCov:
				del self.coverageDict[pos]

	def findIslands(self):
		'''find fully covered islands (min coverage for each base = 10, default)
		coordinates are inclusive on both ends
		islands = [[start1,end1],[start2,end2]...]'''

		islands = []
		positions = self.coverageDict.keys()
		i = 0
		while i < len(positions):
			island_start = positions[i]
			try:
				while positions[i+1] == positions[i] + 1:
					i += 1
				island_end = positions[i]
			except IndexError:
				island_end = positions[i]
			i += 1
			islands.append([island_start,island_end])
		self.islands = islands

	def filterIslands(self, minMedianCoverage=10, minLength=100):
		'''filter out islands based on median coverage and length'''
		keptIslands = []
		for island in self.islands:
			coverage = [self.coverageDict[position] for position in range(island[0],island[1]+1)]
			median = self.getMedian(coverage)
			if (median >= minMedianCoverage) and ((island[1] - island[0]) > minLength):
				keptIslands.append(island)
		self.filteredIslands = keptIslands

	def findJunctions(self, minCoverage=5):
		''' find junctions based on split-mapped reads (Gene.splicedReads)
		format for junction coordinates:
		start and end correspond to the the first
		or last coordinate of the adjacent exons
		Dictionary scheme:
		{(junction start, junction end): coverage (int)}'''

		# build junction dictionary / compute coverage
		junctions = dict()
		for read in self.splicedReads:
			if len(read.blocks) == 2:
				junctionStart = read.blocks[0][1]
				junctionEnd = read.blocks[1][0] + 1
				try:
					junctions[(junctionStart,junctionEnd)] += 1
				except KeyError:
					junctions[(junctionStart,junctionEnd)] = 1

			elif len(read.blocks) == 3:
				for i in range(2):
					junctionStart = read.blocks[i][1]
					junctionEnd = read.blocks[i+1][0] + 1
					try:
						junctions[(junctionStart,junctionEnd)] += 1
					except KeyError:
						junctions[(junctionStart,junctionEnd)] = 1

		# filter according to minCoverage criteria
		for junction in junctions.keys():
			if junctions[junction] < minCoverage:
				del junctions[junction]

		self.junctionDict = junctions

	def makeJunctionBed(self):
		'''Print a junctions.bed file to stdout for
		genome browser visualization'''

		for junction in self.junctionDict.keys():
			toWrite = [self.chr, str(junction[0]), str(junction[1]), self.ID + '.junction',
			str(self.junctionDict[junction]), self.strand]
			print '\t'.join(toWrite)

	def getMaxCoverage(self,start,end):
		'''compute maximum coverage of a region'''
		highest = 0
		for position in range(start,end+1):
			try:
				if self.coverageDict[position] > highest:
					highest = self.coverageDict[position]
			except KeyError:
				continue
		return highest

	def findNearestCanonical(self,start,end):
		'''returns the nearest canonical exon within 10MB,
		returns None if nothing is found'''
		distance = 99999999
		for exon in self.canonicalExons:
			if not self.isOverlapping(exon,(start,end)):
				diff = abs(end - exon[1])
				if diff < distance:
					distance = diff
					closest = exon
		if distance == 99999999:
			return None
		return closest

	def isOverlapping(self,exon1,exon2):
		'''checks if two exons are overlapping
		takes in exon tuples, returns bool'''
		
		if exon1[0] < exon2[0] and exon1[1] < exon2[0]:
			return False
		elif exon2[0] < exon1[0] and exon2[1] < exon1[0]:
			return False
		else:
			return True

	def getMaxCoverageNearestExon(self,start,end):
		'''finds the nearest non-overlapping canonical exon,
		returns max coverage of that exon'''

		closestExon = self.findNearestCanonical(start,end)
		try:
			return self.getMaxCoverage(closestExon[0],closestExon[1])
		except TypeError:
			return None

	def getPercentCovered(self,start,end):
		'''returns 1. percentage of bases that are in the coverage dictionary,
		i.e., have a minimum coverage of 10 reads (default)
		2. median coverage of the bases that are covered'''

		allPos = set([i for i in range(start,end+1)])
		covered = set([i for i in self.coverageDict.keys()]).intersection(allPos)
		percentCov = (len(covered)*100)/float(len(allPos))
		coverage = [self.coverageDict[i] for i in covered]
		try:
			median = self.getMedian(coverage)
		except IndexError:
			# coverage is an empty list
			median = 1
		return percentCov,median

	def intronicRegionContaining(self,start,end):
		'''checks if a region contains a canonical intron,
		returns introns that it contains, and empty set if none'''

		containedIntrons = set([])
		for intron in self.introns:
			if self.isOverlapping(intron,(start,end)):
				containedIntrons.add(intron)
		return containedIntrons

	def getMedian(self,numbers):
	    return (sorted(numbers)[int(round((len(numbers) - 1) / 2.0))] + 
	    	sorted(numbers)[int(round((len(numbers) - 1) // 2.0))]) / 2.0

	def getMedianCoverage(self,start,end):
		'''returns median coverage of a region'''
		coverage = []
		for i in range(start,end+1):
			try:
				coverage.append(self.coverageDict[i])
			except KeyError:
				pass
		if len(coverage) > 0:
			return self.getMedian(coverage)
		else:
			return 1

	def coverageCriteriaMet(self,start,end, medianFrac=.15, percentCovMin=98):
		'''check if the region:
		a. is 98% (default) covered with reads
		b. contains any introns, if so:
			1. median coverage of intronic regions are at least 15% (default) of the max coverage of the neighboring exon
		c. if not, then require that the median coverage of the region is at least 15% (default) of the max coverage of
			the nearest non-overlapping exon
		'''

		percentCov,medianCov = self.getPercentCovered(start,end)
		if percentCov >= percentCovMin:
			try:
				maxCovNearestExon = self.getMaxCoverageNearestExon(start,end)
				if medianCov/float(maxCovNearestExon) >= medianFrac:
					if len(self.intronicRegionContaining(start,end)) > 0:
						introns = self.intronicRegionContaining(start,end)
						retainedIntrons = set([])
						for intron in introns:
							if (intron[1]-intron[0]) < (end-start):
								medianCov = self.getMedianCoverage(intron[0],intron[1])
								try:
									maxCovNearestExon = self.getMaxCoverageNearestExon(intron[0],intron[1])
									if medianCov/float(maxCovNearestExon) >= medianFrac:
										retainedIntrons.add(intron)
								except (TypeError,ZeroDivisionError):
									maxCov = self.getMaxCoverage(self.chrStart,self.chrEnd)
									if medianCov/float(maxCov) >= medianFrac:
										retainedIntrons.add(intron)
							else:
								medianCov = self.getMedianCoverage(start,end)
								try:
									maxCovNearestExon = self.getMaxCoverageNearestExon(start,end)
									if medianCov/float(maxCovNearestExon) >= medianFrac:
										retainedIntrons.add(intron)
								except (TypeError,ZeroDivisionError):
									maxCov = self.getMaxCoverage(self.chrStart,self.chrEnd)
									if medianCov/float(maxCov) >= medianFrac:
										retainedIntrons.add(intron)

						if retainedIntrons == introns:
							return True
						else:
							#exon = (start,end)
							#print exon

							# not all introns contained are truly retained by our criteria
							return False
					else:
						# doesn't contain intronic regions, but meets median / percent coverage criteria
						return True
				else:
					# overall median cov not 15% of nearest exon
					return False

			# if an exon within 10 megabases is not found, 
			# or max coverage comes out to 0 (no covered positions found
			# inside of the nearest canonical exon), use max cov of first canonical exon
			except (TypeError, ZeroDivisionError):
				exonlst = list(self.canonicalExons)
				exonlst.sort()
				if self.strand=='+':
					exon = exonlst[0]
				elif self.strand=='-':
					exon = exonlst[-1]
				maxCov = self.getMaxCoverage(exon[0],exon[1])
				try:
					if medianCov/float(maxCov) >= medianFrac:
						return True
					else:
						return False
				except ZeroDivisionError:
					maxCov = self.getMaxCoverage(self.chrStart,self.chrEnd)
					if medianCov/float(maxCov) >= medianFrac:
						return True
					else:
						return False
		else:
			return False

	def findInternalExons(self,i_medianFrac=.15,i_percentCovMin=98):
		'''identify boundaries of candidate non-terminal exons
		NOTE: "JUNCTION START" OR "JUNCTION END" REFERS TO
		START AND END WITH REFERENCE TO GENOMIC COORDINATES,
		NOT STRAND ORIENTATION OF THE GENE, i.e. any coordinate
		which is smaller in value than its partner will be considered
		a junction start'''

		junctionStarts = set([])
		junctionEnds = set([])

		# adding the TSS as the first junction end
		# may want to include (known) alternative TSSs, will come back to this...
		if self.strand == '+':
			junctionEnds = set([self.chrStart])
		elif self.strand == '-':
			junctionStarts = set([self.chrEnd])

		# build sets of junction ends and starts
		for junction in self.junctionDict.keys():
			junctionStarts.add(junction[0])
			junctionEnds.add(junction[1])

		# for each junction end, we want to find the NEXT junction
		# such that the region between the junction end, and the adjacent
		# junction start meets our coverage critera

		if self.strand == '+':
			exons = set([])
			for end in junctionEnds:
				for start in junctionStarts:
					if start > end:
						if self.coverageCriteriaMet(end,start,medianFrac=i_medianFrac,percentCovMin=i_percentCovMin):
							exons.add((end,start))
						else:
							continue
					else:
						continue

		elif self.strand == '-':
			exons = set([])
			for end in junctionEnds:
				for start in junctionStarts:
					if start > end:
						if self.coverageCriteriaMet(end,start,medianFrac=i_medianFrac,percentCovMin=i_percentCovMin):
							exons.add((end,start))
						else:
							continue
					else:
						continue

		self.internalExons = exons
		self.junctionEnds = junctionEnds
		self.junctionStarts = junctionStarts

	def findTerminalExons(self, wiggle=45,t_medianFrac=.15,t_percentCovMin=98):
		'''Identify terminal exons by:
		Find ends of coverage islands that are unique when compared to internal exons,
		as in outside the range of internal exon end +/- wiggle, and find regions between
		junction ends and unique coverage island ends that meet coverage criteria for novel
		exons'''

		terminalExons = set([])
		internalExonStarts = set([i[0] for i in self.internalExons])
		internalExonEnds = set([i[1] for i in self.internalExons])
		internalExonStartsWiggle = set([i for j in internalExonStarts for i in range(j-wiggle, j+wiggle+1)])
		internalExonEndsWiggle = set([i for j in internalExonEnds for i in range(j-wiggle, j+wiggle+1)])

		try:
			islandEnds = set([i[1] for i in self.filteredIslands])
			islandStarts = set([i[0] for i in self.filteredIslands])
		except NameError:
			islandEnds = set([i[1] for i in self.islands])
			islandStarts = set([i[0] for i in self.islands])
		
		if self.strand == '+':
			candidateTerminalEnds = islandEnds - internalExonEndsWiggle
			for end in self.junctionEnds:
				for terminalEnd in candidateTerminalEnds:
					if terminalEnd > end:
						if self.coverageCriteriaMet(end,terminalEnd,medianFrac=t_medianFrac,percentCovMin=t_percentCovMin):
							terminalExons.add((end,terminalEnd))

		elif self.strand == '-':
			candidateTerminalStarts = islandStarts - internalExonStartsWiggle
			for start in self.junctionStarts:
				for terminalStart in candidateTerminalStarts:
					if terminalStart < start:
						if self.coverageCriteriaMet(terminalStart,start,medianFrac=t_medianFrac,percentCovMin=t_percentCovMin):
							terminalExons.add((terminalStart,start))

		if self.isCoding:
			terminalList = list(terminalExons)
			for exon in terminalList:
				if self.terminalUniqueness(exon[0],exon[1]):
					pass
				else:
					terminalExons.remove(exon)
		else:
			pass

		self.terminalExons = terminalExons

	def terminalUniqueness(self,start,end):
		'''For coding genes, terminal exons that create variation
		in ways that do not affect protein sequence are collapsed
		into one single terminal exon'''
		true = 0
		false = 0
		# check if this exon will affect protein sequence
		# if not, return that it is non-unique
		if self.strand == '+':
			if start > (self.stopCodonStart + 2):
				return False
		elif self.strand == '-':
			if end < self.stopCodonStart:
				return False

		for exon in self.canonicalTerminalExons:
			if self.strand == '+':
				if start == exon[0]:
					if end > (self.stopCodonStart + 2):
						false += 1
					else:
						true += 1
				elif end == exon[1]:
					if start > (self.stopCodonStart + 2):
						false += 1
					else:
						true += 1
				else:
					true += 1
			elif self.strand == '-':
				if end == exon[1]:
					if start < self.stopCodonStart:
						false += 1
					else:
						true += 1
				elif start == exon[0]:
					if end < self.stopCodonStart:
						false += 1
					else:
						true += 1
				else:
					true += 1

		if true > 0:
			return True
		else:
			return False

	def buildGraph(self):
		'''Construct the splice graph, in the form a dicitonary
		{exon: set(neighboring exons)}'''

		# compute set union:
		allExons = self.internalExons | self.terminalExons | self.canonicalExons
		self.allExons = allExons

		graph = dict()
		for junction in self.junctionDict.keys():

			upstreamExons = set([])
			downstreamExons = set([])
			for exon in allExons:
				if self.strand == '+':
					# find exons that end in same position as junction begins
					if exon[1] == junction[0]:
						upstreamExons.add(exon)
						
					# find exons that begin in the same position as junction ends
					if junction[1] == exon[0]:
						downstreamExons.add(exon)
						
				elif self.strand == '-':
					if exon[0] == junction[1]:
						upstreamExons.add(exon)
						
					# find exons that begin in the same position as junction ends
					if junction[0] == exon[1]:
						downstreamExons.add(exon)
						
			# require that junction is paired with an exon on both ends
			if len(upstreamExons) < 1 or len(downstreamExons) < 1:
				pass
			else:
				for exon1 in upstreamExons:
					for exon2 in downstreamExons:
						try:
							graph[exon1].add(exon2)
						except KeyError:
							graph[exon1] = set([exon2])


		for exon in self.terminalExons:
			graph[exon] = set([exon])
		for exon in self.canonicalTerminalExons:
			graph[exon] = set([exon])

		self.spliceGraph = graph

	def getTSSExons(self):
		'''find exons which will be starting nodes
		by checking if they contain TSSs'''

		self.TSSExons = set([])
		for exon in self.allExons:
			if self.strand == '+':
				for TSS in self.TSSs:
					if exon[0] == TSS:
						self.TSSExons.add(exon)
			elif self.strand == '-':
				for TSS in self.TSSs:
					if exon[1] == TSS:
						self.TSSExons.add(exon)

		# pick the largest TSS-containing exon per start codon
		# prevent redundancy in peptide sequences

		if self.isCoding:
			startCodonDict = dict()
			for exon in self.TSSExons:
				for startCodon in self.startCodonTSSdict.keys():
					TSSs = self.startCodonTSSdict[startCodon]
					if self.strand == '+':
						for TSS in TSSs:
							if exon[0] == TSS:
								try:
									startCodonDict[startCodon].add(exon)
								except KeyError:
									startCodonDict[startCodon] = set([exon])
					elif self.strand == '-':
						for TSS in TSSs:
							if exon[1] == TSS:
								try:
									startCodonDict[startCodon].add(exon)
								except KeyError:
									startCodonDict[startCodon] = set([exon])

			uniquelyCodingTSSexons = set([])
			for startCodon in startCodonDict.keys():
				exonLst = list(startCodonDict[startCodon])
				largest = exonLst[0]
				# use most upstream TSS exon for each start codon
				for exon in exonLst:
					if (exon[1]-exon[0]) > (largest[1]-largest[0]):
						largest = exon
				uniquelyCodingTSSexons.add(largest)

			self.TSSExons = self.TSSExons & uniquelyCodingTSSexons
	
	def getPaths(self,minGeneCov=40,p_minCov=10,minIslandLength=100,minJunctionCov=5,
		exonCovFrac=.1,exonPercentCovMin=97,wiggle=45,codingOnly=True):
		'''for all combinations of 
		TSS-containing exons and terminal nodes,
		returns all possible paths using Gene.DFS method'''
		
		self.filterCoverageDict(minCov=p_minCov)
		self.defineIntrons()
		self.removeUncoveredCanonical()
		if len(self.canonicalExons) == 0:
			# if there are no canonical exons that are at least 85% covered with reads
			self.paths = None
			return None
		self.findJunctions(minCoverage=minJunctionCov)
		self.findIslands()
		self.filterIslands(minLength=minIslandLength)
		self.findInternalExons(i_medianFrac=exonCovFrac,i_percentCovMin=exonPercentCovMin)
		self.findTerminalExons(wiggle=wiggle,t_medianFrac=exonCovFrac,t_percentCovMin=exonPercentCovMin)
		self.buildGraph()
		self.getTSSExons()

		self.terminalExons = self.terminalExons | self.canonicalTerminalExons

		allPaths = []
		for startingNode in self.TSSExons:
			if startingNode in self.spliceGraph.keys():
				for terminalNode in self.terminalExons:
					paths = self.DFS(startingNode, terminalNode)
					if len(paths) > 0:
						allPaths.append(paths)
		allPaths = [lst for pathlst in allPaths for lst in pathlst]		
		self.paths = allPaths

	def checkPath(self,path):
		'''checks that at least one position in at least one exon in the path contains
		the three start codon coordinates'''

		posLst = []
		for startCodon in self.startCodonStarts:
			posLst.extend([i for i in range(startCodon,startCodon+3)])
		found = 0
		for exon in path:
			for pos in posLst:
				if exon[0] <= pos <= exon[1]:
					found += 1

		if found >= 3:
			return True
		else:
			return False

	def writeGTF(self):
		''' Use Gene.paths to write a GTF-formatted file to stdout'''

		self.pathCounts = []
		t = 1
		if self.paths == None:
			return None
		for path in self.paths:
			if len(path) > 0:
				if self.isCoding:
					if not self.checkPath(path):
						continue

				if self.strand == '-':
					TSS = path[0][1]
					path = path[::-1]
				elif self.strand == '+':
					TSS = path[0][0]
				
				transcript_str = self.ID + '.neo.' + str(t)
				field9 = 'gene_id "' + self.ID + '"; transcript_id "' + transcript_str + '";'
				
				if self.isCoding:
					# get start codon associated with TSS of the path
					for startCodon in self.startCodonTSSdict.keys():
						if TSS in self.startCodonTSSdict[startCodon]:
							start_codon = startCodon
					# print start codon line
					fields = [str(self.chr), 'NeoGraph.1', 'start_codon', str(start_codon), str(start_codon+2), '.',
					self.strand, '.', field9]
					print '\t'.join(fields)

				for exon in path:

					fields = [str(self.chr), 'NeoGraph.1', 'exon', str(exon[0]), str(exon[1]),
					'.', self.strand, '.', field9]
					print '\t'.join(fields)
				t += 1

		self.pathCounts.append([self.ID,len(self.canonicalExons),t])

	def DFS(self, start_vertex, end_vertex, path=[]):
		'''recursive depth-first search to find all paths
		between two nodes in the splice graph'''

		path = path + [start_vertex]
		if start_vertex == end_vertex:
			return [path]
		paths = []
		try:
			for vertex in self.spliceGraph[start_vertex]:
				if vertex not in path:
					extended_paths = self.DFS(vertex, end_vertex, path)
					for p in extended_paths:
						paths.append(p)
		except KeyError:
			# this is necessary because of one edge case:
			# an internal exon whose terminal neighbor was removed because it was deemed
			# functionally irrelevant
			# this internal exon would not be in the splicegraph as a result of the removal
			# of its neighboring terminal exon
			pass
		return paths

def parseRef(gtf):
	'''Parse the GTF format reference annotation
	Make sure file includes 'start_codon' lines if 'coding only' option is set True'''

	gtf = open(gtf).readlines()
	genes = dict()
	i = 0
	gene_id = gtf[0].rstrip('\n').split('\t')[8].split('gene_id')[1].split('"')[1]
	while i < len(gtf):
		gene_id = gtf[i].rstrip('\n').split('\t')[8].split('gene_id')[1].split('"')[1]
		fields = gtf[i].rstrip('\n').split('\t')
		chr = fields[0]
		strand = fields[6]
		genes[gene_id] = [chr,strand,[]]
		while ((gtf[i].rstrip('\n').split('\t')[8].split('gene_id')[1].split('"')[1]) == gene_id):
			fields = gtf[i].rstrip('\n').split('\t')
			gene_id = gtf[i].rstrip('\n').split('\t')[8].split('gene_id')[1].split('"')[1]
			transcript_id = gtf[i].rstrip('\n').split('\t')[8].split('transcript_id')[1].split('"')[1]
			transcriptInfo = [set([]),set([]),set([])]
			while ((gtf[i].rstrip('\n').split('\t')[8].split('transcript_id')[1].split('"')[1]) == transcript_id):
				transcript_id = gtf[i].rstrip('\n').split('\t')[8].split('transcript_id')[1].split('"')[1]
				fields = gtf[i].rstrip('\n').split('\t')
				if fields[2] == 'exon':
					exon = (int(fields[3]),int(fields[4]))
					transcriptInfo[0].add(exon)
				elif fields[2] == 'start_codon':
					startCodonStart = int(fields[3])
					transcriptInfo[1].add(startCodonStart)
				elif fields[2] == 'stop_codon':
					stopCodonStart = int(fields[3])
					transcriptInfo[2].add(stopCodonStart)
				i += 1
				if i >= len(gtf):
					break
			genes[gene_id][2].append(transcriptInfo)
			if i >= len(gtf):
				break
	del gtf
	return genes

def getNeoAttributes(refInfo, bam, minGeneCov=40, codingOnly=True):
	chr = refInfo[0]
	strand = refInfo[1]
	exons = set([])
	TSSs = set([])
	TTSs = set([])
	terminalExons = set([])
	stopCodonStarts = set([])
	startCodonStarts = set([])
	startCodonTSSdict = dict()

	for transcript in refInfo[2]:
		try:
			startCodon = list(transcript[1])[0]
		except IndexError:
			pass

		exons = exons.union(transcript[0])
		startCodonStarts = startCodonStarts.union(transcript[1])
		stopCodonStarts = stopCodonStarts.union(transcript[2])
		exonStarts = [i[0] for i in transcript[0]]
		exonEnds = [i[1] for i in transcript[0]]
		exonList = list(exons)
		exonList.sort(key=lambda x: x[0])
		if strand == '+':
			TSS = min(exonStarts)
			TSSs.add(TSS)
			TTSs.add(max(exonEnds))
			terminalExons.add(exonList[-1])
			try:
				startCodonTSSdict[startCodon].add(TSS)
			except (KeyError, NameError):
				try:
					startCodonTSSdict[startCodon] = set([TSS])
				except NameError:
					pass

		elif strand == '-':
			TSS = max(exonEnds)
			TSSs.add(TSS)
			TTSs.add(min(exonStarts))
			terminalExons.add(exonList[0])
			try:
				startCodonTSSdict[startCodon].add(TSS)
			except (KeyError, NameError):
				try:
					startCodonTSSdict[startCodon] = set([TSS])
				except NameError:
					pass

	if strand == '+':
		chrStart = min(list(TSSs))
		chrEnd = max(list(TTSs))
	elif strand == '-':
		chrStart = min(list(TTSs))
		chrEnd = max(list(TSSs))

	if codingOnly:
		if len(startCodonStarts) == 0:
			# no codon starts
			return None
	
	try:
		reads = set([read for read in bam.fetch(chr, int(chrStart), int(chrEnd))])
	except ValueError:
		# chromosome name not found in alignment file
		# print 'chr name not found in bam file'
		return None

	# check sufficient coverage
	totalExonKb = 0
	for exon in exons:
		totalExonKb += (exon[1]-exon[0])/float(1000)

	readsPerKb = len(reads)/totalExonKb
	if readsPerKb < minGeneCov:
		return None

	spliced_reads = set([])
	coverage = dict()

	# remove splice junction reads from the wrong strand
	read_list = list(reads)
	i = 0
	for read in read_list:
		if 'XS' in [i[0] for i in read.tags]:
			if ('XS', strand) in read.tags:
				spliced_reads.add(read)
			else:
				reads.remove(read)

		for position in read.positions:
			try:
				coverage[position + 1] = coverage[position + 1] + 1
			except KeyError:
				coverage[position + 1] = 1

	attr = [chr, chrStart, chrEnd, strand, reads, spliced_reads, coverage, exons,terminalExons,stopCodonStarts,TSSs, startCodonStarts, startCodonTSSdict]
	return attr

def readInput(gtf,bam):
	genes = parseRef(gtf)
	bam = pysam.AlignmentFile(bam, 'rb')
	return genes, bam

def makeGene(gene, geneDict, bam, minGeneCov=40, codingOnly=True):
	args = getNeoAttributes(geneDict[gene],bam, minGeneCov=minGeneCov, codingOnly=codingOnly)
	if args == None:
		return None
	args = [gene] + args
	neo = Gene(*args)
	return neo

def runNeograph(gtf, bam, options=[40,98,10,.15,5,45,True]):
	'''Function for running command-line version of neograph'''
	minGeneCov = options[0]
	minPercentCoverage = options[1]
	coverageMin = options[2]
	minMedianFrac = options[3]
	minJunctionCov = options[4]
	wiggle = options[5]
	codingOnly = options[6]

	genes,bam = readInput(gtf,bam)
	for gene in genes:
		neo = makeGene(gene, genes, bam, minGeneCov=minGeneCov, codingOnly=codingOnly)
		try:
			neo.getPaths(minGeneCov=minGeneCov,p_minCov=coverageMin,minJunctionCov=minJunctionCov,
				exonCovFrac=minMedianFrac,exonPercentCovMin=minPercentCoverage,wiggle=wiggle)
		except AttributeError:
			continue
		neo.writeGTF()

def main():
	'''main function for command-line version of neograph'''

	usage = 'usage: python neograph.py --gtf [path/to/ref.gtf] --bam [path/to/indexed/alignment.bam] [options] > [output]'
	parser = OptionParser(usage=usage)
	parser.add_option('-g','--gtf', dest='gtf', help='reference annotation (GTF format) including start/stop codon info')
	parser.add_option('-b','--bam', dest='bam', help='RNA-seq alignment BAM file (make sure it is indexed)')
	parser.add_option('-r','--minGeneCov', dest='minGeneCoverage', help='minimum number of reads per exonic kb to attempt to build splice graph, default = 40 reads/kb',
		default=40)
	parser.add_option('-p', '--minPercentCoverage', dest='minPercentCoverage', help='minimum percentage of covered bases required for an exon NOTE: "covered" is restricted to bases with >= 10 reads (unless otherwise specified by -c option), default = 97',
		default=97)
	parser.add_option('-c', '--coverageMin', dest='covDictMin', default=10, help='minimum coverage for a base to be considered "covered" (see -p option)')
	parser.add_option('-f', '--minMedianFrac', dest='minMedianFrac',default=.15, help='in conjunction with minimum percent coverage criteria, the median coverage of a candidate exon must be at least this fraction of the max coverage of the nearest non-overlapping canonical exon, default = .15')
	parser.add_option('-j', '--minJunctionCov', dest='minJunctionCov', default=5, help='minimum junction coverage for inclusion in the splice graph, default = 5')
	parser.add_option('-w', '--terminalExonWindow', dest='terminalExonWindow', default=45, 
		help='the window of bp upstream/downstream of a canonical terminal exon against which an island is compared to determine if it is novel. See wiki for more details. default = 45')
	parser.add_option('-t', '--codingOnly', dest='codingOnly', default=True, help='only produce splice graphs for protein coding genes, default = True')

	(options, args) = parser.parse_args()
	if options.gtf == None or options.bam == None:
		print usage
		print ''
		print 'You must provide a BAM file and GTF reference annotation'
		print '"python neograph.py -h" for additional options'
		sys.exit()

	gtf_file = options.gtf
	bam_file = options.bam
	minGeneCov = int(options.minGeneCoverage)
	try:
		minPercentCoverage = int(options.minPercentCoverage)
	except ValueError:
		print 'Please enter a percentage for "--minPercentCoverage", i.e. a number between 1 and 100'
		sys.exit()	
	coverageMin = int(options.covDictMin)
	minMedianFrac = float(options.minMedianFrac)
	if minMedianFrac >= 1:
		print 'Please enter a decimal fraction, i.e. a float from 0 to 1, for "--minMedianFrac"'
		sys.exit()
	junctionCov = int(options.minJunctionCov)
	wiggle = int(options.terminalExonWindow)
	codingOnly = options.codingOnly in [True, 'True', 'TRUE','true', '1', 'yes', 'Yes']
	options = [minGeneCov,minPercentCoverage,coverageMin,minMedianFrac,junctionCov, wiggle, codingOnly]

	runNeograph(gtf_file,bam_file,options=options)

if __name__=='__main__':
	main()