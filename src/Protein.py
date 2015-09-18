from collections import defaultdict
import Xrefs

__author__ = 'etai'


class Protein:
    id = None
    name = None
    dbxrefs = None
    seq = None
    cds = None
    cdsid = None
    cdstype = None
    alignseq = None
    alignseqid = None
    tax = None
    organism = None
    organismHost = None
    proteinExistence = None

    def __init__(self, id=None, name=None, dbxrefs=None, seq=None, cds=None, tax=None, organism=None, organismHost=None,
                 proteinExistence=None, xmlrecord=None):
        if xmlrecord:
            dbxrefs = [i.split(':') for i in xmlrecord.dbxrefs]
            crs = defaultdict(list)
            for x in dbxrefs:
                crs[x[0]].append(x[1])
            self.id = xmlrecord.id
            self.name = xmlrecord.name
            self.dbxrefs = Xrefs.Xrefs(crs)
            self.seq = xmlrecord.seq
            self.tax = xmlrecord.annotations.get("taxonomy")
            self.organism = xmlrecord.annotations.get("organism")
            self.organismHost = xmlrecord.annotations.get("organismHost")
            if xmlrecord.annotations.get("proteinExistence") is not None:
                self.proteinExistence = xmlrecord.annotations.get("proteinExistence")[0]
        else:
            self.id = id
            self.name = name
            self.dbxrefs = xrefs(dbxrefs)
            self.seq = seq
            self.cds = cds
            self.tax = tax
            self.organism = organism
            self.organismHost = organismHost
            self.proteinExistence = proteinExistence

    def __str__(self):
        return "id = " + str(self.id) + "\nname = " + str(self.name) + "\nnucleotide dbxrefs = " + str(
            self.dbxrefs.getNucleotideRefs()) + "\ntaxonomy = " + '<-'.join(self.tax) + "\norganism = " + str(
            self.organism) + "\norganismHost = " + str(self.organismHost) + "\nproteinExistence = " + str(
            self.proteinExistence) + "\nprotein seq = " + str(self.seq) + "\nprotein cds = " + str(
            self.cds) + "\ncds type = " + str(self.cdstype) + "\nalignseq = " + str(self.alignseq)
