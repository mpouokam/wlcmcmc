import wlcmc
import os
import numpy as np
import itertools
import sys
import tempfile

diameter = 4.5
ONLYONETRANSLATIONVECTOR = True;

homflydict = {}
homflydict['z^-1.( -a^-1 +a )']='0_1'

homflydict['z^-1.( -a^-3 +a^-1 ) + z.a^-1']='2_1'
homflydict['z^-1.( -a +a^3 ) -z.a']='2_1s'

homflydict['z^-1.( -a^3 +a^5 ) + z.( -a-a^3 )']='4_1'
homflydict['z^-1.( -a^-5 +a^-3 ) + z.( a^-3 +a^-1 )']='4_1s'


def writeknotstosingle3dcfile(arrayofvertices,filename):
    with open(filename,"w") as file: #opens/closes file
        for vertices in arrayofvertices: #each element is a list of vertices
            file.write("%d "%len(vertices)) #start each line with the number of vertices
            np.savetxt(file, vertices.reshape(1,-1)) #reshape as a rowvector to print on one line

def getrandomunknots(numberoftrials= 20,mcsteps = 50000):

    wlc = wlcmc.wormlikechain(numberofkuhnlengths=13,segmentsperkuhnlength=10,diameter=4.5)

    wlc.domontecarlosteps(50000)

    # wlc.testforautocorr(numberofmcsteps=20000,numberofsamples=100)

    wlcs_vertices = []

    for i in range(numberoftrials):
        wlc.domontecarlosteps(mcsteps)
        wlc.centeratorigin()
        wlcs_vertices.append(np.copy(wlc.vertices))


    #TODO: Move this all to a KnotPlot interaction library. from vertices as numpy array get knottypes

    sequencefile = tempfile.mktemp()
    writeknotstosingle3dcfile(wlcs_vertices,sequencefile)

    gaussfile = tempfile.mktemp(suffix='.egc') # because knotplot adds the .egc extension
    relaxamount = 0
    forknotplot = """
    gauss open %s
    seq open %s
    until seqend \"ago %d;gauss;seq +\"
    seq close
    """%(gaussfile,sequencefile,relaxamount)
    writestringtofile("/tmp/forknotplot.kps",forknotplot,flags = 'w')

    os.system("cat /tmp/forknotplot.kps | knotplot -nog")
    knottypes = subprocess.check_output("xinger -znb +r %s | homfly -- | idknot -k -- "%(gaussfile),shell=True).split()

    os.remove(gaussfile)
    os.remove(sequencefile)

    i = 0
    for line in knottypes:
        if line == '0_1':
            i+=1
        else:
            wlcs_vertices.pop(i)
            print "%d is a knot"%i

    print "Found %d unknots"%i    


    return wlcs_vertices

def distancebetweenlinks(m,n):
    import wlcmclib
    biglist = np.concatenate((m,n))
    output = 1000000 #a big number that is larger than everything
    for i in range(len(m)):
        for j in range(len(n)):
            # print i,(i+1)%len(m),j+len(m),(j+1)%len(n)+len(m)
            distance = wlcmclib.distancebetweenedges(biglist,i,(i+1)%len(m),j+len(m),(j+1)%len(n)+len(m))
            if distance<output:
                output = distance
    # print "Distance: %f"%output,
    return output

def writestringtofile(filename,string,flags = 'w'):
    file = open(filename, flags)
    file.write(string)
    file.close()



def writelinktofile(filename,m,n):
    with open(filename, 'w') as file:
        np.savetxt(file, m)
        file.write("\n")
        np.savetxt(file, n)

import random

import subprocess
def gethomflypolynomial(infilename,relaxamount=500):
    with tempfile.NamedTemporaryFile() as gaussfile:
        gausscode = subprocess.check_output("echo 'load %s; ago %d; gauss' | knotplot -nog | grep '|'"%(infilename,relaxamount),shell=True)[:-1]
        if gausscode == '|':
            return 'z^-1.( -a^-1 +a )'
        else:
            return subprocess.check_output("echo '%s' | homfly --"%(gausscode),shell=True)[:-1]


def getlinktype(m,n):
    with tempfile.NamedTemporaryFile() as linkforknotplotfp:
        writelinktofile(linkforknotplotfp.name,m,n)
        homfly = gethomflypolynomial(linkforknotplotfp.name)

    # print homfly
    if homfly in homflydict:
        return homflydict[homfly]
    else:
        return homfly

def getrotationmatricies():
    rotationmatrices = []

    for i in range(2):
        for j in range(2):
            for k in range(2):
                idenwithreflections = np.matrix(((np.identity(3),-np.identity(3),)[i][0],(np.identity(3),-np.identity(3),)[j][1],(np.identity(3),-np.identity(3),)[k][2]))
                # print idenwithreflections
                for p in itertools.permutations(idenwithreflections):
                    rotationmatrices.append(np.concatenate(p))

    return rotationmatrices

def getxlationvectors():
    translationvectors = []
    for row in np.identity(3):
        translationvectors.append(row)
    for row in -np.identity(3):
        translationvectors.append(row)

    if ONLYONETRANSLATIONVECTOR ==True:
        return translationvectors[:1]

    return translationvectors



def getdata(wlcunknots,distances=np.arange(0,400,50)):
    import itertools 

    totals = np.zeros(distances.shape)
    overlap = np.zeros(distances.shape)
    link2_1 = np.zeros(distances.shape)
    link2_1s = np.zeros(distances.shape)
    link4_1 = np.zeros(distances.shape)
    link4_1s = np.zeros(distances.shape)
    link0_1 = np.zeros(distances.shape)

    unitvectors=getxlationvectors() #TODO: add more here
    rotations=getrotationmatricies()

    for d in range(len(distances)):
    # for d in range(1):
        print "Running for distance %f"%distances[d]
        with open("output/dataforcatsim%f.txt"%distances[d],"a") as file:
            for vector in (unitvectors[:1],unitvectors)[d!=0]: #in the case that d==0, only use one xlation vector
                for rotation in rotations:
                    # print vector*distances[d]
                    for (m,n) in itertools.combinations(wlcunknots,2): #iterate over all pairs of elements
                        # print distances[d]
                        nadjusted = n*rotation+vector*distances[d]
                        totals[d] += 1 #increment count of total
                        if distancebetweenlinks(m,nadjusted)<diameter:
                            print "overlap",
                            file.write("overlap\n")
                            overlap[d]+=1
                        else:
                            linktype = getlinktype(m,nadjusted)
                            print linktype,
                            file.write("%s\n"%linktype)
                            if linktype == '2_1':
                                link2_1[d]+=1
                            if linktype == '2_1s':
                                link2_1s[d]+=1
                            elif linktype == '4_1':
                                link4_1[d]+=1
                            elif linktype == '4_1s':
                                link4_1s[d]+=1
                            elif linktype == '0_1':
                                link0_1[d]+=1

    print "Totals: ",totals
    print "Overlap: ",overlap
    print "link0_1: ",link0_1
    print "link2_1: ",link2_1
    print "link2_1s: ",link2_1s
    print "link4_1: ",link4_1
    print "link4_1s: ",link4_1s
    return (totals,overlap,link0_1,link2_1,link2_1s,link4_1,link4_1s)


if __name__ == '__main__':
    wlcunknots = getrandomunknots()
    # print wlcunknots
    getdata(wlcunknots,distances=np.arange(00,400,50))
