__author__ = 'mikeknowles'
"""
This script will download genomes for the orgamism provided
Below is are the required modules to run
"""

from Bio import Entrez
import time
from urllib2 import URLError, HTTPError

start = 0

import urllib, os
from threading import Thread
from Queue import Queue
from argparse import ArgumentParser

dqueue = Queue()  # Queue for multithreading using pthreads


def parser(dqueue):
    """
    :param dqueue: Queue for multiple downloads at once
    :return: fasta WGS or complete chromosome
    """
    while True:  # this is a loop to help multithreading
        gi, count, path = dqueue.get()  # retrieve tuple from the queue
        connected = False  # Start not connected
        while not connected:  # Loop to overcome connection issues with a terrible network
            try:
                summary = Entrez.esummary(db="sra",
                                          id=gi,)
                connected = True
            except IndexError:
                print "[%s] Unable to download" % (time.strftime("%H:%M:%S"))
                connected = True
            except ValueError:
                print "[%s] No records found in handle" % (time.strftime("%H:%M:%S"))
                connected = True
            except HTTPError:
                print "[%s] Unable to connect trying again in 5 seconds" % (time.strftime("%H:%M:%S"))
                time.sleep(5)
        record = Entrez.read(summary)  # reads the fasta data from NCBI
        summary.close()  # Good practise
        acc = record[0]["Runs"][10:19]
        filepath = os.path.join(path, acc)
        print "[%s] Downloading and converting #%i %s to fastq format..." % (time.strftime("%H:%M:%S"), count, acc)
        os.system("fastq-dump -I --split-files -O %s %s" % (filepath, acc))
        dqueue.task_done()

def dlthreads(email, organism, path):
    organism = organism.replace('_', '+')
    count = 0
    if not os.path.isdir(path):
        os.mkdir(path)
    path = os.path.join(path, '')
    Entrez.email = email
    searchterm = "(%s)" \
                 % (organism)
    search = Entrez.esearch(db="sra",
                            term=searchterm,
                            retmax=10000)
    print search.url
    search = Entrez.read(search)
    print "[%s] Found %s genome records" % (time.strftime("%H:%M:%S"), search['Count'])
    # for gi in search["IdList"]:
    #     summary = Entrez.esummary(db="sra",
    #                               id=gi,)
    #     summary = Entrez.read(summary)
    #     print summary
    #     acc = summary[0]["Runs"][10:19]
    #     # /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
    #     print "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra" % (acc[:3], acc[:6], acc, acc)

    for i in range(3):
        threads = Thread(target=parser, args=(dqueue,))
        threads.setDaemon(True)
        threads.start()
    for i in search["IdList"]:
        count += 1
        dqueue.put((i, count, path))
    dqueue.join()
'''
Parser for arguments test
'''
parse = ArgumentParser(description='Download genomes for organism')
parse.add_argument('--version', action='version', version='%(prog)s v0.3')
parse.add_argument('-q', '--query', required=True, help='Query for genome database separated by plus sign(s)')
parse.add_argument('-e', '--email', required=True, help='A valid email address is required')
parse.add_argument('-o', '--output', required=True, help='Specify output directory')

args = vars(parse.parse_args())
dlthreads(args['email'], args['query'], args['output'])
