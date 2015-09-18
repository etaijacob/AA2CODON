from collections import defaultdict
import re

__author__ = 'etai'

class Xrefs(dict):
	def __init__(self, *args, **kwargs):
		dict.__init__(self, *args, **kwargs)
		self.__dict__ = self

	def getNucleotideRefs(self):
		d = defaultdict(list)
		for k, v in self.iteritems():
			if k == 'EMBL':
				for i in v:
					e = i.split(';')
					if len(e) >= 3 and e[1] != '-':
						d[k + "-" + e[2]].append(e[1])
			if k == 'RefSeq':
				for i in v:
					e = i.split(';')
					if len(e) >= 2 and e[1] != '-':
						d[k].append(e[1])
			if k.startswith("Ensembl"):
				for i in v:
					e = i.split(';')
					#print e
					if len(e) >= 1 and e[0] != '-':
						d[k].append(e[0])
		return d

	def getNucleotideRefsByPriority(self):
		refs = self.getNucleotideRefs()
		#print refs
		reflist = []
		reftype = []
		refseq = [ref for ref in refs['RefSeq'] if ref.startswith("NM_")]
		refseq.extend([ref for ref in refs['RefSeq'] if ref.startswith("XM_")])
		reflist.extend(refseq)
		reftype.extend([ "refseqn" for i in range(0,len(refseq)) ])

		emblm = [ref for ref in refs['EMBL-mRNA']]
		reflist.extend(emblm)
		reftype.extend([ "emblcds" for i in range(0,len(emblm)) ])

		ensembl = [ref for ref in refs['Ensembl']]
		reflist.extend(ensembl)
		reftype.extend([ "ensembltranscript" for i in range(0,len(ensembl)) ])

		ensemblg = []
		for k in refs.keys():
			if re.match("Ensembl\w+", k):
				for i in refs[k]:
					ensemblg.append(i)
		reflist.extend(ensemblg)
		reftype.extend([ "ensemblgenomestranscript" for i in range(0,len(ensemblg)) ])

		emblg = [ref for ref in refs['EMBL-Genomic_DNA']]
		reflist.extend(emblg)
		reftype.extend([ "emblcds" for i in range(0,len(emblg)) ])
		#print reflist
		#print reftype
		return (reflist, reftype)


