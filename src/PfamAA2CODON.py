import os
import sys
from Bio import AlignIO, SwissProt, SeqIO
import re
from collections import defaultdict
from MyUniprotXML import UniprotIterator
import subprocess
import time, json
from ConfigParser import RawConfigParser, ConfigParser
import FetchUtils, Xrefs
from Protein import Protein
from sys import platform as _platform

__author__ = 'etai'

#############################
# External tools
#############################
#
# tranalign:
# ---------
#User should change tranalign location according to needs:
#For windows download and execute ftp://emboss.open-bio.org/pub/EMBOSS/windows/mEMBOSS-6.5.0.0-setup.exe.
#For linux/mac download and extract ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.5.7.tar.gz.
# curl:
# ----
# On windows you need to download and install curl http://curl.haxx.se/download.html.
# The installation configures the curl executable path so no additional configurations are needed.


if _platform == "linux" or _platform == "linux2":
    tranalignCmd = "tranalign"
elif _platform == "darwin":
    tranalignCmd = "tranalign"
elif _platform == "win32":
    tranalignCmd = "C:/mEMBOSS/tranalign.exe"


# Output level
outputLevel = 1
# Debug level
debugLevel = 0
fetchme = True


def main(argv=None):
    global fetchme

    if argv is None:
        argv = sys.argv

    pfam = argv[1]

    if len(argv) > 2:
        if argv[2] == "-noFetch":
            fetchme = False
        else:
            fetchme = True

    if fetchme:
        print "Fetching data.."
    else:
        print "No fetching.."

    FetchUtils.soap_setup()

    proteins = None

    if fetchme:
        if os.path.isfile(pfam + ".uniprot"):
            print "File: " + pfam + ".uniprot data already exists"
            accs = getACCsFromSTHfile(pfam)
        else:
            print "File: " + pfam + ".uniprot data does not exist.\nRetrieving data from EBI.."
            accs = fetchUniprotEntries(pfam)

    proteins = getDataFromUniprot(pfam, getACCsFromSTHfile(pfam))
    # for id, protein in proteins.iteritems():
    # print protein, "\n"

    proteins = poolProteinCDSs(proteins, pfam)
    assignAlignseqsToProteins(proteins, pfam)


def tranalign(pfam):
    alignment = AlignIO.read(open(pfam), "stockholm")
    print dir(alignment)
    for record in alignment:
        print dir(record)
        print record.id, record.seq


# 1. RefSeq: NM_,	(mRNA, first one available on the list)
# 2. RefSeq: XM_, (predicted mRNA, first one available on the list)
# 3. EMBL: ACCESSION_NUMBER, EMBL|ACCESSION_NUMBER|PROTEIN_ID|-|mRNA, (mRNA)
# 4. Ensembl: T suffix, Ensembl|..T|..P|..G, first one available on the list
# 5. EMBL: ACCESSION_NUMBER, EMBL|ACCESSION_NUMBER|PROTEIN_ID|-|Genomic_DNA/RNA <--- check if has CDS in db
# 6. RefSeq: NC_, (Complete genomic molecules) <-- check if has CDS
def poolProteinCDSs(proteins, fname_suffix):
    if not os.path.isfile(fname_suffix + ".config"):
        writeNewConfigFile(fname_suffix + ".config")
    proteins, pool, n, rip, ncds = assignCDSsToProteins(proteins, fname_suffix)
    nprev = n
    prevrip = rip
    print "There are additional " + str(n) + " protein cdss to retrieve."

    if n > 0 and fetchme:
        conf = getConfigData(fname_suffix + ".config")
        prioritypool = int(conf['all']['prioritypool'])
        for reftype in pool.keys():
            updateConfigFile(fname_suffix + ".config", 'all', 'prioritypool', prioritypool + 1)
        for reftype, refs in pool.iteritems():
            print "\n" + str(reftype) + ":"
            print "---------------"
            print reftype + " Pool length = " + str(len(refs))
            if len(refs) > 0:
                fetchEntries(reftype, refs.keys(), fname_suffix + "." + reftype)

        proteins, pool, n, rip, ncds = assignCDSsToProteins(proteins, fname_suffix)
        print "Tried to download a total of " + str(nprev) + " proteins. Succeded to download " + str(
            nprev - n - (rip - prevrip))
        print "Total of " + str(ncds) + " proteins with a cds, out of " + str(len(proteins)) + " (" + str(
            rip) + " have no cdss)"
        return proteins
    else:
        print "No more CDSs to download"
        print "Total of " + str(ncds) + " proteins with a cds, out of " + str(len(proteins)) + " (" + str(
            rip) + " have no cdss)"
        return proteins


# TODOS: change it to protein id and not ref list
# Assign cdss and then check if more cdss need to be retrieved
def getPoolData(proteins, fname_suffix):
    pool = {'refseqn': defaultdict(list),
            'emblcds': defaultdict(list),
            'ensembltranscript': defaultdict(list),
            'ensemblgenomestranscript': defaultdict(list)}

    if not os.path.isfile(fname_suffix + ".config"):
        writeNewConfigFile(fname_suffix + ".config")
        for id, protein in proteins.iteritems():
            reflist, reftype = protein.dbxrefs.getNucleotideRefsByPriority()
            pool[reftype[0]][reflist[0]].append(id)
        return pool, {'refseqn': defaultdict(list),
                      'emblcds': defaultdict(list),
                      'ensembltranscript': defaultdict(list),
                      'ensemblgenomestranscript': defaultdict(list)}
    else:
        conf = getConfigData(fname_suffix + ".config")
        prioritypool = int(conf['all']['prioritypool'])
        # print conf
        cdss = {}
        for reftype in pool.keys():
            cdss[reftype] = read_cdss(fname_suffix + "." + reftype)
        adds = 0
        for id, protein in proteins.iteritems():
            reflist, reftype = protein.dbxrefs.getNucleotideRefsByPriority()
            if len(reflist) > prioritypool:
                # print prioritypool, reftype
                if reflist[prioritypool] not in cdss[reftype[prioritypool]].keys():
                    if len(reflist) > prioritypool + 1:
                        pool[reftype[prioritypool + 1]][reflist[prioritypool + 1]].append(id)
                        adds = adds + 1
        if adds > 0:
            print "Retrieving prioritypool = " + str(prioritypool + 1) + " .Total of " + str(adds) + ".."
            for reftype in pool.keys():
                updateConfigFile(fname_suffix + ".config", 'all', 'prioritypool', prioritypool + 1)
        else:
            print "No additional CDSs to retrieve (prioritypool = " + str(prioritypool) + ")."
        return pool, cdss


def assignAlignseqsToProteins(proteins, pfam):
    alignment = AlignIO.read(open(pfam), "stockholm")
    print dir(alignment)
    alignseqs = defaultdict(dict)
    for record in alignment:
        acc = re.sub("\.\d$", '', record.annotations["accession"])
        alignseqs[acc]["seq"] = record.seq
        alignseqs[acc]["id"] = record.id
    # print acc, record.id
    cnt = 0
    with open(pfam + ".cdss", "w") as fc, open(pfam + ".aligndprots", "w") as fp:
        for id, protein in proteins.iteritems():
            print protein
            if proteins[id].cds != None and id in alignseqs.keys():
                cnt = cnt + 1
                print id, str(proteins[id].name)

                # if re.match("[A|C|G|T]*", str(proteins[id].cds)):
                fp.write(">" + str(alignseqs[id]["id"]) + "\n")
                fp.write(str(alignseqs[id]["seq"] + "\n"))

                fc.write(">" + str(alignseqs[id]["id"]) + "\t" + str(proteins[id].cdsid) + "\n")
                fc.write(str(proteins[id].cds) + "\n")
                # else:
                # print "ERROR!!! " + str(alignseqs[id]["id"])

                proteins[id].alignseq = alignseqs[id]["seq"]
                proteins[id].alignseqid = alignseqs[id]["id"]
            else:
                print "No cds!"

    fn_cdss = pfam + ".cdss"
    fn_alignprots = pfam + ".aligndprots"
    fn_tranalign = pfam + ".tranalignout"
    if os.path.isfile(fn_cdss) and os.path.isfile(fn_alignprots):
        print "cdss and alignprots files exist. Doing alignment projection.."

        #with open(fn_tranalign, "w") as f:
        cmd = tranalignCmd + " -auto -asequence " + fn_cdss + " -bsequence " + fn_alignprots + " -outseq " + fn_tranalign  # + " -table 11"
        p = subprocess.Popen(cmd, shell=True, stdout=False, stderr=subprocess.STDOUT)
        ret_code = p.wait()
    print cnt


def assignCDSsToProteins(proteins, fname_suffix):
    retypes = ['refseqn', 'emblcds', 'ensembltranscript', 'ensemblgenomestranscript']
    cdss = {}
    pool = {'refseqn': defaultdict(list),
            'emblcds': defaultdict(list),
            'ensembltranscript': defaultdict(list),
            'ensemblgenomestranscript': defaultdict(list)}

    conf = getConfigData(fname_suffix + ".config")
    prioritypool = int(conf['all']['prioritypool'])
    print "Pooling priority number: " + str(prioritypool)
    n = 0
    rip = 0
    ncds = 0
    for reftype in retypes:
        cdss[reftype] = read_cdss(fname_suffix + "." + reftype)
    for id, protein in proteins.iteritems():
        has_cds = False
        reflist, reftype = protein.dbxrefs.getNucleotideRefsByPriority()
        for priority in range(0, len(reflist)):
            if reflist[priority] in cdss[reftype[priority]].keys():
                # print "Protein " + id + " has a cds: " + reflist[priority]
                proteins[id].cds = cdss[reftype[priority]][reflist[priority]]
                proteins[id].cdsid = reflist[priority]
                proteins[id].cdstype = reftype[priority]
                has_cds = True
                ncds = ncds + 1
                break
        if not has_cds:
            print "#Protein " + id + " has no cds."
            if len(reflist) > prioritypool + 1:
                pool[reftype[prioritypool + 1]][reflist[prioritypool + 1]].append(id)
                n = n + 1
            else:
                rip = rip + 1
                print "#No more cdss for " + id

    return proteins, pool, n, rip, ncds


def getConfigData(cfgfname):
    Config = RawConfigParser()
    Config.read(cfgfname)

    def ConfigSectionMap():
        dict1 = defaultdict(dict)
        sections = Config.sections()
        for section in sections:
            options = Config.options(section)
            for option in options:
                try:
                    dict1[section][option] = Config.get(section, option)
                    if dict1[section][option] == -1:
                        DebugPrint("skip: %s" % option)
                except:
                    print("exception on %s!" % option)
                    dict1[section][option] = None
        return dict1

    return ConfigSectionMap()


def writeNewConfigFile(cfgfname):
    sections = ['all', 'refseqn', 'emblcds', 'ensembltranscript', 'ensemblgenomestranscript']
    print "No configuration file yet. Producing one.."
    with open(cfgfname, 'w') as f:
        Config = ConfigParser()
        for section in sections:
            Config.add_section(section)
            Config.set(section, "prioritypool", -1)
            Config.set(section, "date", time.strftime("%d/%m/%Y %H:%M:%S"))
        Config.write(f)


def updateConfigFile(cfgfname, section, varname, value):
    print "Updating config. file: " + cfgfname + " " + section, varname, value
    Config = RawConfigParser()
    Config.read(cfgfname)
    Config.set(section, varname, value)
    Config.set(section, "date", time.strftime("%d/%m/%Y %H:%M:%S"))
    with open(cfgfname, 'w') as f:
        Config.write(f)
    return


def fetchEntries(reftype, pool, fname, add=True):
    if os.path.isfile(fname) and add == False:
        print "File: " + fname + " data already exists. Reading sequences from file.."
        return 0  # read_cdss(fname)
    else:
        print "File: " + fname + " data does not exist or partial.\nRetrieving relevant data from EBI.."
        l = range(0, len(pool))
        chunks = FetchUtils.getChunks(l, 100)
        print reftype
        if reftype == 'refseqn':
            for i in range(0, len(chunks)):
                print i
                ids = [pool[j] for j in chunks[i]]
                # fetchBatch <dbName> <id1,id2,...> [formatName [styleName]] [options...]
                print ','.join(ids)
                try:
                    entries = None
                    entries = FetchUtils.returnFetchBatch(reftype, str(','.join(ids)), 'default', 'default')
                    with open(fname, "a") as f:
                        if len(entries) > 0:
                            f.write(''.join(str(v) for v in entries))
                        else:
                            print "Nothing to write"
                except FetchUtils.suds.WebFault as detail:
                    print "MY DETAILS!!!!!!!!!!!!!!!!!!!!:"
                    print detail

                # entries = returnFetchBatch(reftype, ','.join(ids), 'default', 'default')
                # with open(fname, "a") as f:
                # f.write(entries)
            with open(fname, "r") as f, open(fname + "2", "w") as fo:
                lastline = None
                for line in f:
                    match = re.match("(.+)<a href=.+>(\w+)</a>(.+)", line)
                    if match:
                        fo.write(match.group(1) + match.group(2) + match.group(3) + "\n")
                    else:
                        fo.write(line)
                        lastline = line
                if lastline != "//":
                    fo.write("//")
                fo.flush()
                os.fsync(fo)
            if os.path.isfile(fname + "2"):
                os.rename(fname + "2", fname)
        elif reftype == 'emblcds':
            for i in range(0, len(chunks)):
                print i
                ids = [pool[j] for j in chunks[i]]
                # fetchBatch <dbName> <id1,id2,...> [formatName [styleName]] [options...]
                print ','.join(ids)
                try:
                    entries = None
                    entries = FetchUtils.returnFetchBatch(reftype, str(','.join(ids)), 'default', 'default')
                    with open(fname, "a") as f:
                        if len(entries) > 0:
                            f.write(''.join(str(v) for v in entries))
                        else:
                            print "Nothing to write"
                except FetchUtils.suds.WebFault as detail:
                    print "MY DETAILS!!!!!!!!!!!!!!!!!!!!:"
                    print detail

                # entries = returnFetchBatch(reftype, ','.join(ids), 'default', 'default')
                # with open(fname, "a") as f:
                # f.write(entries)

        #the original code used for the research article utilized the beta version
        elif reftype == 'ensembltranscript':
            l = range(0, len(pool))
            chunks = FetchUtils.getChunks(l, 50)
            print reftype
            DEVNULL = open(os.devnull, 'wb')
            for i in range(0, len(chunks)):
                print "retrieveing: " + str(i) + "\n"
                ids = [pool[j] for j in chunks[i]]
                if _platform == "linux" or _platform == "linux2":
                    sep = '"'
                elif _platform == "darwin":
                    sep = '"'
                elif _platform == "win32":
                    sep = '\\\"'
                myids = (', '.join(sep + item + sep for item in ids))
                if _platform == "linux" or _platform == "linux2":
                    curlCmd = "curl -H'Accept: application/json' -H'Content-type: application/json' -XPOST --data '{\"ids\" : [" + myids + "] }' 'http://rest.ensembl.org/sequence/id?object_type=transcript&type=cds'"
                elif _platform == "darwin":
                    curlCmd = "curl -H'Accept: application/json' -H'Content-type: application/json' -XPOST --data '{\"ids\" : [" + myids + "] }' 'http://rest.ensembl.org/sequence/id?object_type=transcript&type=cds'"
                elif _platform == "win32":
                    curlCmd = "curl -H\"Accept: application/json\" -H\"Content-type: application/json\" -XPOST --data \"{\\\"ids\\\" : [" + myids + "] }\" \"http://rest.ensembl.org/sequence/id?object_type=transcript&type=cds\""
                print curlCmd
                with open(fname + "tmp.dnld2", "w") as f:
                    # print "Writing " + gene + " to " + fname
                    p = subprocess.Popen(
                        curlCmd,
                        shell=True, stdout=f, stderr=DEVNULL)

                    ret_code = p.wait()
                    with open(fname + "tmp.dnld2", "r") as fo2:
                        myres = json.loads(fo2.read())
                        print myres
                        if (len(myres) > 0):
                            with open(fname, "a") as f:
                                for i_id in range(0, len(myres)):
                                    f.write(">" + myres[i_id]["id"] + "\n" + myres[i_id]["seq"] + "\n")
            DEVNULL.close()
        elif reftype == 'ensemblgenomestranscript':
            l = range(0, len(pool))
            chunks = FetchUtils.getChunks(l, 50)
            print reftype
            DEVNULL = open(os.devnull, 'wb')
            for i in range(0, len(chunks)):
                print "retrieveing: " + str(i) + "\n"
                ids = [pool[j] for j in chunks[i]]
                if _platform == "linux" or _platform == "linux2":
                    sep = '"'
                elif _platform == "darwin":
                    sep = '"'
                elif _platform == "win32":
                    sep = '\\\"'
                myids = (', '.join(sep + item + sep for item in ids))
                if _platform == "linux" or _platform == "linux2":
                    curlCmd = "curl -H'Accept: application/json' -H'Content-type: application/json' -XPOST --data '{\"ids\" : [" + myids + "] }' 'http://rest.ensemblgenomes.org/sequence/id?object_type=transcript&type=cds'"
                elif _platform == "darwin":
                    curlCmd = "curl -H'Accept: application/json' -H'Content-type: application/json' -XPOST --data '{\"ids\" : [" + myids + "] }' 'http://rest.ensemblgenomes.org/sequence/id?object_type=transcript&type=cds'"
                elif _platform == "win32":
                    curlCmd = "curl -H\"Accept: application/json\" -H\"Content-type: application/json\" -XPOST --data \"{\\\"ids\\\" : [" + myids + "] }\" \"http://rest.ensemblgenomes.org/sequence/id?object_type=transcript&type=cds\""
                print curlCmd
                with open(fname + "tmp.dnld2", "w") as f:
                    # print "Writing " + gene + " to " + fname
                    p = subprocess.Popen(
                        curlCmd,
                        shell=True, stdout=f, stderr=DEVNULL)  # subprocess.STDOUT)
                    ret_code = p.wait()
                    with open(fname + "tmp.dnld2", "r") as fo2:
                        myres = json.loads(fo2.read())
                        print myres
                        if (len(myres) > 0):
                            with open(fname, "a") as f:
                                for i_id in range(0, len(myres)):
                                    f.write(">" + myres[i_id]["id"] + "\n" + myres[i_id]["seq"] + "\n")
            DEVNULL.close()
        else:
            print reftype + " No such reftype."

    return read_cdss(fname)


def read_cdss(fname):
    cds = {}
    handle = None
    if os.path.isfile(fname):
        handle = open(fname, "rU")
    else:
        return cds

    if fname.endswith('emblcds'):
        for record in SeqIO.parse(handle, "embl"):
            cds[record.id] = record.seq
    if fname.endswith('refseqn'):
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    loc = re.match("\[(?P<from>\d+):(?P<to>\d+)\]\((?P<direction>\S)\)", str(feature.location))
                    if (loc):
                        cds[record.id] = record.seq[int(loc.group("from")):int(loc.group("to"))]
                        break
    if fname.endswith('ensembltranscript'):
        for record in SeqIO.parse(handle, "fasta"):
            # print record.id, record.seq
            cds[record.id] = record.seq
    if fname.endswith('ensemblgenomestranscript'):
        for record in SeqIO.parse(handle, "fasta"):
            # print record.id, record.seq
            cds[re.sub("-tr", "tr", record.id)] = record.seq

    handle.close()

    return cds


def getACCsFromSTHfile(pfam):
    pfam_ids = {}
    alignment = AlignIO.read(open(pfam), "stockholm")
    for record in alignment:
        acc = re.sub("\.\d+", '', record.annotations["accession"])
        pfam_ids[acc] = record.id

    print "There are " + str(len(pfam_ids)) + " Aligned entries in total (original sth file)"
    return pfam_ids.keys()


def fetchUniprotEntries(pfam, newaccs=None):
    print "Fetching uniprot entries as follows.."
    uniprot_file = pfam + ".uniprot"
    already_exists = False
    if os.path.isfile(uniprot_file):
        already_exists = True
        lastline = tail(uniprot_file, 1)
        lastline = lastline.rstrip()
        if lastline == "</uniprot>":
            cmd = "head -n -1 " + uniprot_file + " > " + uniprot_file + ".2; " + "mv " + uniprot_file + ".2 " + uniprot_file
            print cmd
            os.system(cmd)

    pfam_ids = newaccs
    if newaccs == None:
        pfam_ids = getACCsFromSTHfile(pfam)
    print pfam_ids
    if pfam_ids:
        l = range(0, len(pfam_ids))
        chunks = FetchUtils.getChunks(l, 50)
        for i in range(0, len(chunks)):
            print "Fetching " + str(i) + " chunk out of " + str(len(chunks))
            sys.stdout.flush()
            accs = [pfam_ids[j] for j in chunks[i]]
            # fetchBatch <dbName> <id1,id2,...> [formatName [styleName]] [options...]
            try:
                uniprots = FetchUtils.returnFetchBatch('uniprot', ','.join(accs), 'xml', 'default')
                with open(pfam + ".uniprot", "a") as f:
                    if not i == 0 or already_exists:
                        uniprots = re.sub(r"\<uniprot.+", '', uniprots)
                        uniprots = re.sub(r"\<\?xml version=.+\>", '', uniprots)
                    if not i == len(chunks) - 1:
                        uniprots = re.sub(r"\</uniprot.*>", '', uniprots)
                    f.write(uniprots)
            except:
                print "Coundnot fetch this chunk"
                continue
            else:
                print "Wow.."
                continue
    return (pfam_ids)


def getDataFromUniprot(pfam, accs):
    global fetchme
    uniprot_file = pfam + ".uniprot"
    print uniprot_file
    lastline = tail(uniprot_file, 1)
    if lastline:
        lastline = lastline.rstrip()
        print "::::::::::: " + lastline
        if lastline != "</uniprot>":
            with open(uniprot_file, "a") as f:
                f.write("</uniprot>\n")

    handle = open(uniprot_file)
    cntr = 0
    proteins = {}
    for record in UniprotIterator(handle):
        proteins[record.id] = Protein(xmlrecord=record)
        cntr = cntr + 1
    # if cntr > 1000:
    # break
    if accs != None:
        accs2get = [acc for acc in accs if acc not in proteins.keys()]
        if len(accs2get) > 0:
            print "fetching: "
            print accs2get
            if fetchme:
                fetchUniprotEntries(pfam, accs2get)
            return getDataFromUniprot(pfam, None)

    print "There are total of " + str(len(proteins)) + " entries mapped in uniprot"
    return proteins


def tail(f, n=1):
    if os.path.isfile(f):
        stdin, stdout = os.popen2("tail -n " + str(n) + " " + f)
        stdin.close()
        lines = stdout.readlines();
        stdout.close()
        if lines:
            return lines[0]
    return None


if __name__ == "__main__":
	sys.exit(main())