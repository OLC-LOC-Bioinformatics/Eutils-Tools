__author__ = 'mikeknowles'
"""
This script will download genomes for the orgamism provided
Below is are the required modules to run
"""

from Bio import SeqIO, Entrez
import time
from urllib2 import URLError, HTTPError

start = 0

import urllib, os, sys
from re import sub
from threading import Thread
from Queue import Queue
from argparse import ArgumentParser
from StringIO import StringIO
from gzip import GzipFile

dqueue = Queue() # Queue for multithreading using pthreads

class KeyboardInterruptError(Exception):
    pass


def parser(dqueue):
    """
    :param dqueue: Queue for multiple downloads at once
    :return: fasta WGS or complete chromosome
    """
    while True:  # this is a loop to help multithreading
        try:
            gi, count, path = dqueue.get()  # retrieve tuple from the queue
            connected = False  # Start not connected
            while not connected:  # Loop to overcome connection issues with a terrible network
                try:
                    fetch = Entrez.efetch(db="nuccore",
                                          id=gi,
                                          rettype='fasta',
                                          retmode='text')
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
            record = SeqIO.read(fetch, "fasta")  # reads the fasta data from NCBI
            fetch.close()  # Good practise to close unused handles
            ecoli = record.description.replace(record.id, "")[1:].split(",")[0]  # string split
            print "[%s] Processing %s GI:%s" % (time.strftime("%H:%M:%S"), ecoli, gi)
            filename = sub('[^A-Za-z0-9]+str[ain.]*[^A-Za-z0-9]*|[^A-Za-z0-9]+subsp\.*[^A-Za-z0-9]*'
                           '|[^A-Za-z0-9]+serovar\.*[^A-Za-z0-9]*|[^A-Za-z0-9]+', '_', ecoli)
            name = filename.replace("_complete_genome", "")
            # Create filename for a cataloging and removing unneccesary bits
            filename = path + name + ".fasta"
            fasta = open(filename, 'w')
            if "_" in record.id:
                urlid = record.id.split("_")[-1]
            else:
                urlid = record.id
            # faster method to retrive genomes
            url = "https://www.ncbi.nlm.nih.gov/Traces/wgs/?download=" + urlid[:5] + "1.1.fsa_nt.gz"
            try:
                source = urllib.urlopen(url).read()
                compressedFile = StringIO()
                compressedFile.write(source)
                #
                # Set the file's current position to the beginning
                # of the file so that gzip.GzipFile can read
                # its contents from the top.
                #
                compressedFile.seek(0)

                decompressedFile = GzipFile(fileobj=compressedFile, mode='rb')
                # with open(filename, 'w') as outfile:
                #     outfile.write(decompressedFile.read())

                fasta.write(decompressedFile.read())
                print "[%s] Downloading and unzipping #%i %s..." % (time.strftime("%H:%M:%S"), count, name)
            except IOError:
                fasta = open(filename, "w")
                print "[%s] Downloading #%i %s..." % (time.strftime("%H:%M:%S"), count, name)
                SeqIO.write(record, fasta, "fasta")
            fasta.close()
        except KeyboardInterrupt:
            raise KeyboardInterruptError()
        dqueue.task_done()


def dlthreads(email, organism, path, length, arg=''):
    organism = organism.replace('_', '+')
    count = 0
    lengthrange = length.split("-")
    if not os.path.isdir(path):
        os.mkdir(path)
    path = os.path.join(path, '')
    Entrez.email = email
    searchterm = "({0:s})+AND+\"{1:d}\"[SLEN]:\"{2:d}\"[SLEN]+srcdb+refseq[prop]+" \
        .format(organism, (int(lengthrange[0]) * 10 ** 6), (int(lengthrange[1]) * 10 ** 6)) + arg
    search = Entrez.esearch(db="nuccore",
                            term=searchterm,
                            retmax=10000)
    # print search.url
    search = Entrez.read(search)
    print "[%s] Found %s genome records" % (time.strftime("%H:%M:%S"), search['Count'])
    for i in range(3):
        threads = Thread(target=parser, args=(dqueue,))
        threads.setDaemon(True)
        threads.start()
    try:
        for i in search["IdList"]:
            count += 1
            dqueue.put((i, count, path))
        dqueue.join()
    except KeyboardInterrupt:
        print "[{0:s}] Got ^C while pool mapping, terminating the pool".format(time.strftime("%H:%M:%S"))
        dqueue.empty()
        print '[{0:s}] pool is terminated'
        sys.exit(127)
    except Exception, e:
        print "[{0:s}] Got exception: {1!r:s}, terminating the pool".format(time.strftime("%H:%M:%S"), e)
        dqueue.empty()
        print "[{0:s}] Pool is terminated".format(time.strftime("%H:%M:%S"))
        sys.exit(127)
'''
Parser for arguments test
'''
parse = ArgumentParser(description='Download genomes for organism')
parse.add_argument('--version', action='version', version='%(prog)s v0.3')
parse.add_argument('-q', '--query', required=True, help='Query for genome database separated by plus sign(s)')
parse.add_argument('-e', '--email', required=True, help='A valid email address is required')
parse.add_argument('-o', '--output', required=True, help='Specify output directory')
parse.add_argument('-l', '--length', required=True, help='The range of length for the full genome, the default is 4-7 Mb for E.coli. The default a range in megabases')
parse.add_argument('-c', '--chromosome', action='store_true', help='Specify additional search parameters')

args = parse.parse_args()
if args.chromosome:
    args = vars(args)
    dlthreads(args['email'], args['query'], args['output'], args['length'], 'gene+in+chromosome[prop]')
else:
    args = vars(args)
    dlthreads(args['email'], args['query']+'[Organism]', args['output'], args['length'])
