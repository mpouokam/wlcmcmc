"""
Created on Oct 14, 2012

@author: brian
"""



from __future__ import print_function

USE_WLCMCLIB = True




import numpy as np
import random
import wlcmclib #shared library of c functions shared via SWIG
from timing import print_timing
import pickle
import scipy.io
import os

# random.seed(37)

import time

class wormlikechain(object):
    """
    A class to run Monte-Carlo simulations for the Worm-Like Chain
    """


    def __init__(self,  numberofsegments=81,
                        segmentsperkuhnlength=10,
                        kuhnlength=100,
                        diameter=5,
                        verbose=True,
                        thetamax=1.3,
                        numberofkuhnlengths=None
                        ):
        """
        Constructor
        """
        if numberofkuhnlengths!=None: #TODO: check if its an integer
            numberofsegments=numberofkuhnlengths*segmentsperkuhnlength

        self.edgelength = kuhnlength / segmentsperkuhnlength 
        
        self.vertices = np.zeros((numberofsegments,3))
        self.vertices[:,0] = np.cos(np.arange(0.0,numberofsegments)*2*np.pi/numberofsegments)
        self.vertices[:,1] = np.sin(np.arange(0.0,numberofsegments)*2*np.pi/numberofsegments)

        multiplyfactor = self.edgelength / np.linalg.norm(self.vertices[1]-self.vertices[0])
        
        self.vertices = self.vertices * multiplyfactor
        
        self.setbendingrigidity(segmentsperkuhnlength)
        self.verbose = verbose
        self.thetamax = thetamax
        self.diameter = diameter
        self.setrandomseedtoclocktime()

    def setbendingrigidity(self,segmentsperkuhnlength=10):
        bendingconstantlookup = [0,0,0.67177455656782,1.2321073947562,1.7603474154772,2.2762843289484,2.7864843908854,3.2935187787792,3.7987547657586,4.3028122687772,4.805740366969,5.3083172041781,5.8105314820542,6.3125693425482,6.81405838005,7.3153574515189,7.81644929059,8.3171957614662,8.8180812466788,9.3191325979456,9.8198472616918,10.320494752223,10.722746712959,11.321787546002,11.822260560226,12.322561101274,12.822967408898,13.323339240011,13.823679599268,14.323991098953,14.823688750868,15.323832765895,15.824597536992,16.324837866329,16.825061784839,17.325269645617,17.825461712154,18.325638213409,18.825799387457,19.325945514385,19.917128129371,20.428022704132,20.805555323082,21.200525295434,21.83364687309,22.333622628985,22.752201204177,23.260810835515,23.768574540077,24.357488916313,24.854837302086,25.352386077399,25.850130007911,26.348061834622,26.84617281159,27.344453160118,27.842892447538,28.341479898832,28.840204649088,29.339055944416,29.838023298442,30.337096610917,30.836266254368,31.335523134081,31.834858726093,32.33426509729,32.833734911161,33.333261422224,33.83283846173,34.332460416797,34.832122204791,35.331819244435,35.831547424863,36.331303073609,36.83108292428,37.33088408454,37.830704004865,38.330540448398,38.830391462174,39.330255349853,39.830130646099,40.330016092625,40.829910615937,41.329813306754,41.829723401059,42.329640262725,42.829563367643,43.329492289287,43.829426685607,44.329366287193,44.8293108866,45.329260328769,45.829214502455,46.329173332568,46.698354370428,47.329104802519,47.995036723518,48.329054621765,46.998827327809,49.329022892324,49.829014007815,50.329009812812,50.829010332326,51.329015587475,51.829025593954,52.3290403608,52.829059889412,53.329084172799,53.829113195031,54.329146930877,54.829185345588,55.329228394835,55.829276024757,56.32932817212,56.829384764568,57.329445720952,57.829510951712,58.329580359335,58.829653838843,59.329731278319,59.829812559461,60.329897558164,60.829986145106,61.330078186345,61.830173543922,62.330272076456,62.830373639751,63.330478087365,63.830585271189,64.330695042009,64.83080725004,65.330921745448,65.831038378857,66.331157001813,66.831277467252,67.331399629921,67.831523346792,68.33164847743,68.831774884361,69.331902433392,69.832030993915,70.332160439192,70.8322906466,71.332421497867,71.832552879269,72.332684681821,72.832816801428,73.332949139022,73.430007648813,73.430007648813,74.827904746871,75.326169871732,75.882429400826,76.382624378312,76.882771764302,77.382873004282,77.882929557404,79.587227064149,78.882914484272,80.160667045551,80.160667045552,80.435354094847,80.935199734245,81.434967600473,81.934658577543,82.434273571277,82.933813507449,83.433279330069,83.932671999784,84.431992492405,84.931241797538,85.430420917351,85.929530865396,86.428572665575,86.927547351167,87.426455963942,87.92529955337,88.424079175879,88.922795894201,89.421450776772,89.920044897195,90.418579333751,90.917055168958,91.415473489191,91.691925388724,91.691925388724,92.870813875949,93.370154941521,93.869490415722,94.368821070265,94.868147647388,95.367470860223,95.86679139326,96.366109902783,96.865427017394,97.364743338508,97.864059440923,98.363375873362,98.862693159065,99.362011796392,99.86133225941,100.36065499855,100.85998044118,101.35930899232,101.8586410352,102.35797693195,102.85731702425,103.35666163393,103.85601106366,104.35536559757,104.85472550189,105.35409102557,105.85346240092,106.35283984425,106.85222355642,107.35161372349,107.85101051733,108.35041409614,108.84982460508,109.34924217682,109.84866693209,110.34809898019,110.84753841957,111.34698533833,111.84643981468,112.34590191751,112.8453717068,113.3448492341,113.84433454303,114.34382766966,114.84332864297,115.34283748523,115.84235421247,116.34187883478,116.84141135674,117.34095177779,117.84050009252,118.34005629109,118.83962035947,119.33919227982,119.83877203076,120.33835958767,120.83795492293,121.33755800624,121.83716880486,122.33678728383,122.83641340622,123.33604713335,123.83568842501,124.33533723967,124.83499353462,125.33465726623,125.83432839007,126.33400686112,126.83369263385,127.33338566247,127.83308590098,128.33279330336,128.83250782367,129.87202208151,130.08960516159,130.34122337028,130.8409843761,131.34074960201,131.84051895298,132.34029233479,132.8400696541,133.33985081856,133.83963573686,134.33942431879,134.83921647533,135.33901211868,135.8388111623,136.33861352102,136.83841911103,137.33822784992,137.83803965672,138.33785445199,138.83767215773,139.33749269754,139.83731599653,140.3371419814,140.83697058045,141.33680172357,141.83663534226,142.33647136965,142.8363097405,143.33615039118,143.83599325971,144.33583828574,144.83568541052,145.33553457695,145.83538572953,146.33523881437,146.83509377917,147.33495057321,147.83480914736,148.33466945404,148.8345314472,149.33439508235,149.83426031647,150.33412710805,150.83399541707,151.33386520494,151.8337364345,152.33360907002,152.83348307716,153.33335842294,153.83323507571,154.33311300517,154.8329921823,155.33287257936,155.83275416986,156.33263692854,156.83252083133,157.33240585534,157.83229197883,158.33217918119,158.83206744292,159.33195674555,159.83184707172,160.33173840504,160.83163073015,161.33152403266,161.83141829911,162.33131351699,162.83120967466,163.33110676137,163.8310047672,164.33090368308,164.83080350073,165.33070421261,165.83060581199,166.33050829281,166.83041164976,167.33031587819,167.83022097411,168.33012693415,168.83003375558,169.32994143623,169.82984997453,170.32975936945,170.82966962046,171.32958072755,171.82949269122,172.32940551242,172.82931919252,173.32923373333,173.82914913709,174.32906540639,174.82898254422,175.32890055389,175.8288194391,176.32873920379,176.82865985225,177.32858138903,177.82850381895,178.32842714711,178.82835137879,179.32827651952,179.82820257504,180.32812955125,180.82805745425,181.32798629029,181.82791606577,182.32784678722,182.82777846131,183.32771109477,183.82764469448,184.32757926737,184.82751482044,185.32745136078,185.82738889552,186.32732743178,186.82726697679,189.96470104051,187.82714912181,188.32709173625,188.82703538827,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,189.96470104051,197.32624297614,193.777390616,191.17246708324,198.82613696379,193.777390616,192.46092268199,200.32604150491,200.82601205416,201.3259837958,201.82595673433,202.32593087408,202.82590621929,203.32588277392,203.82586054181,204.32583952662,204.82581973185,205.32580116079,205.82578381653,206.32576770207,206.82575282008,207.32573917317,207.82572676368,208.32571559383,208.82570566559,209.32569698077,209.82568954097,210.32568334763,210.82567840196,211.325674705,211.8256722576,212.3256710604,212.82567111384,213.3256724182,213.82567497357,214.32567877979,214.82568383655,215.32569014337,215.82569769953,216.32570650415,216.82571655617,217.32572785433,217.82574039714,218.32575418302,218.82576921011,219.32578547639,219.82580297972,220.32582171771,220.82584168782,221.32586288732,221.82588531328,222.32590896265,222.82593383219,223.32595991845,223.82598721786,224.32601572665,224.82604544088,225.32607635648,225.8261084692,226.3261417746,226.82617626813,227.32621194505,227.8262488005,228.32628682942,228.82632602665,229.28955692358,229.28955692358,230.3456732775,230.84801310054,231.35034556804,231.85267042685,232.35498742813,232.85729632717,233.35959688338,233.8618888604,234.36417202598,234.866446152,235.36871101452,235.8709663937,236.37321207381,230.29477865032,237.37767349443,237.87988882397,238.3820936324,238.88428772441,239.38647090867,239.88864299786,240.39080380872,240.89295316186,241.39509088189,241.89721679741,242.3993307409,242.90143254871,243.4035220611,243.90559912218,244.40766357984,244.90971528586,245.41175409573,245.91377986875,246.41579246792,246.91779175993,247.41977761519,247.92174990774,248.42370851527,248.92565331902,249.42758420383,249.92950105816,250.43140377381,250.93329224624,251.43516637422,251.93702606009,252.43887120947,252.94070173136,253.44251753821,253.9443185456,254.44610467254,254.94787584115,256.20393975099,256.20393975099,256.58463951643,257.08729659833,257.58993648072,258.0925590222,258.59516408441,259.09775153193,259.60032123227,260.10287305599,260.60540687645,261.10792256989,261.61042001543,262.11289909502,262.61535969338,263.11780169801,263.6202249992,264.12262948986,264.62501506571,265.12738162503,265.62972906878,266.13205730052,266.63436622642,267.13665575517,267.63892579803,268.14117626868,268.64340708338,269.14561816078,269.64780942195,270.14998079035,270.65213219183,271.15426355456,271.65637480901,272.15846588801,272.66053672654,273.16258726193,273.6646174336,274.16662718322,274.66861645463,275.17058519373,275.6725333486,276.17446086933,276.67636770812,277.17825381912,277.68011915854,278.18196368457,278.68378735727,279.1855901387,279.6873719928,280.18913288538,280.69087278408,281.19259165836,281.69428947954,282.19596622065,282.69762185647,283.19925636354,283.7008697201,284.20246190604,284.70403290295,285.20558269398,285.707111264,286.20861859933,286.71010468797,287.21156951943,287.71301308465,288.21443537621,288.71583638805,289.21721611562,289.71857455579,290.21991170683,290.7212275684,291.2225221415,291.72379542856,292.22504743322,292.72627816054,293.22748761677,293.72867580949,294.22984274744,294.7309884407,295.23211290044,295.73321613915,296.23429817032,296.73535900874,297.23639867022,297.73741717177,298.23841453141,298.73939076828,299.2403459026,299.74127995556,300.24219294943,300.74308490746,301.24395585386,301.74480581389,302.24563481369,302.74644288034,303.24723004189,303.74799632726,304.24874176627,304.74946638963,305.25017022883,305.75085331629,306.25151568526,306.75215736971,307.2527784045,307.75337882524,308.25395866831,308.75451797084,309.25505677071,309.75557510652,310.25607301761,310.75655054397,311.25700772631,311.75744460605,312.25786122518,312.75825762643,313.2586338531,313.75898994921,314.25932595922,314.75964192836,315.25993790243,315.76021392768,316.26047005104,316.76070631995,317.26092278249,317.76111948708,318.26129648283,318.76145381932,319.26159154661,319.76170971524,320.26180837624,320.76188758117,321.26194738194,321.76198783107,322.26200898135,322.76201088613,323.26199359909,323.76195717446,324.2619016667,324.7618271308,325.2617336221,325.76162119636,326.26148990961,326.76133981837,327.26117097941,327.76098344994,328.26077728747,328.76055254984,329.26030929523,329.76004758212,330.25976746935,330.75946901597,331.25915228144,331.75881732553,332.25846420815,332.7580929896,333.25770373043,333.75729649145,334.25687133375,334.75642831865,335.25596750778,335.75548896293,336.25499274614,336.75447891977,337.25394754626,337.75339868843,338.25283240926,338.75224877187,339.25164783965,339.75102967623,340.25039434534,340.749741911,341.2490724374,341.7483859888,342.2476826298,342.74696242508,343.2462254395,343.74547173814,344.24470138617,344.74391444898,345.24311099209,345.74229108117,346.24145478194,346.74060216053,347.23973328289,347.73884821541,348.23794702428,348.73702977615,349.23609653763,349.73514737542,350.23418235645,350.73320154769,351.23220501624,351.73119282937,352.23016505438,352.72912175868,353.22806300989,353.72698887563,354.22589942366]
        self.bendingrigidity = bendingconstantlookup[segmentsperkuhnlength]

    def writetofilein3dcoordsformat(self,filename='myfile.3dcoords',flags='ab'):
        f = open(filename,flags)
        mystr = '%d '%len(self.vertices)
        for i in self.vertices.reshape(-1):
            mystr += '%.9f '%i

        print (mystr, file = f)
        f.close()

        # print (str,file = f)

    def savestateasmatlabfile(self,filename='currentstate.mat'):
        scipy.io.savemat(filename, {    'vertices':self.vertices,
                                        'thetamax':self.thetamax,
                                        'verbose':self.verbose,
                                        'diameter':self.diameter,
                                        'edgelength':self.edgelength,
                                        'bendingrigidity':self.bendingrigidity,
                                        },oned_as='row')

    def loadstatefrommatlabfile(self,filename='currentstate.mat'):
        tmp = scipy.io.loadmat(filename)
        self.vertices = np.array(tmp['vertices'],order='C')
        self.thetamax = tmp['thetamax'][0,0] # loadmat loads 3.4 as array([[3.4]]) for some reason, this gets rid of the brackets
        self.verbose = [False,True][int(tmp['verbose'][0,0])] # reads 0 as False and 1 as True
        self.diameter = tmp['diameter'][0,0]
        self.edgelength = tmp['edgelength'][0,0]
        self.bendingrigidity = tmp['bendingrigidity'][0,0]


    if USE_WLCMCLIB:    #do crankshaft    
        def docrankshaft(self,m,n,theta):
            wlcmclib.docrankshaft(self.vertices,m,n,theta)
    else:
        def docrankshaft(self,m,n,theta):
            n = (m+n)%len(self.vertices)

            if m==n or m-n==1 or n-m==1:
                return

            axisofrotation = self.vertices[m] - self.vertices[n]
            axisofrotation /= np.linalg.norm(axisofrotation)
            
            if m<n:
                self.vertices[m:n]=vrodrot(self.vertices[m:n],self.vertices[m],axisofrotation,theta)
            else:
                self.vertices[m:]=vrodrot(self.vertices[m:],self.vertices[m],axisofrotation,theta)
                self.vertices[:n]=vrodrot(self.vertices[:n],self.vertices[m],axisofrotation,theta)
    
    if USE_WLCMCLIB:        
        def dorandomcrankshaft(self):
            wlcmclib.dorandomcrankshaft(self.vertices,self.thetamax)
    else:
        def dorandomcrankshaft(self):
            m = random.randint(0,len(self.vertices)-1)
            n = random.randint(2,len(self.vertices)-2)
            theta = random.uniform(-self.thetamax,self.thetamax)
            self.docrankshaft(m,n,theta)
    
    def dosinglemontecarlo(self):
        oldenergy = self.getenergy()
        oldvertices = np.copy(self.vertices)
        
        self.dorandomcrankshaft()
        
        newenergy = self.getenergy()
        
        if self.iscollision():
            # reject
            del self.vertices #avoid memory leaks
            self.vertices = oldvertices
            return False
        if (newenergy < oldenergy):
            # accept new conformation. do nothing
            return True
        elif random.random() < np.exp(-(newenergy-oldenergy)):
            # also accept with this probability
            return True
        else:
            # reject
            del self.vertices #avoid memory leaks
            self.vertices = oldvertices
            return False
    
    if USE_WLCMCLIB:
        @print_timing
        def domontecarlosteps(self,numberofsteps):
            passes = wlcmclib.domontecarlosteps(self.vertices,numberofsteps,self.bendingrigidity,self.thetamax,self.diameter, self.edgelength)
            print ("ratio: %d:%d"%(passes,numberofsteps-passes))
    else:
        @print_timing
        def domontecarlosteps(self,numberofsteps):
            acceptances = 0
            rejections = 0
            for i in xrange(numberofsteps):
                if self.dosinglemontecarlo():
                    acceptances += 1
                else:
                    rejections += 1
        
            print ("ratio: %d:%d"%(acceptances,rejections))
    
    
    if USE_WLCMCLIB: 
        def getenergy(self):
            return wlcmclib.getwlcenergy(self.vertices,self.bendingrigidity)
            # NOTE: I compared this c function call and the python version, and they seem to be similar.
            # In [83]: wlc.domontecarlosteps(2000); wlcmc.wlcmclib.getwlcenergy(wlc.vertices,wlc.bendingrigidity)-wlc.getenergy()
            # ratio: 756:1244
            # domontecarlosteps took 2279.741 ms
            # Out[83]: 1.3244516594568267e-11
    else:
        def getenergy(self):
            edges = np.zeros((len(self.vertices),3))
            edges[:-1] = self.vertices[1:]-self.vertices[:-1] 
            edges[-1] = self.vertices[0]-self.vertices[-1]
            
            angles = np.zeros((len(self.vertices)))
            angles[:-1] = vdot(edges[1:],edges[:-1]) 
            angles[-1] = np.dot(edges[0],edges[-1])
            
            angles = angles / self.edgelength / self.edgelength
            
            angles = np.arccos(angles)
            
            return 0.5*self.bendingrigidity * (angles**2).sum()
        
    
    if USE_WLCMCLIB:
        def iscollision(self):
            return wlcmclib.iscollision(self.vertices,self.diameter,self.edgelength)

    def showme(self):
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pylab as plt
        
        fig = plt.figure()
        ax = Axes3D(fig)

        tograph = np.concatenate((self.vertices,[self.vertices[0]]))

        ax.plot(    tograph[:,0],
                    tograph[:,1],
                    tograph[:,2], 'b-o')
        plt.show()

    def showmemayavi(self):
        import mayavi.mlab

        tograph = np.concatenate((self.vertices,[self.vertices[0]]))

        mayavi.mlab.plot3d( tograph[:,0],
                            tograph[:,1],
                            tograph[:,2], tube_radius=self.diameter/2, colormap='Spectral')
        mayavi.mlab.show()

    def centeratorigin(self):
        self.vertices -= self.vertices.sum(0)/len(self.vertices)

    def getradiusofgyration(self):
        #SQRT ( (1/n)*sum (rk - rmean)^2 ) 
        return np.sqrt(((self.vertices - self.vertices.sum(0)/len(self.vertices))**2).sum()/len(self.vertices))
    def testforautocorr(self,numberofmcsteps=12000,numberofsamples=40):
        output = []
        for i in range(numberofsamples):
            self.domontecarlosteps(numberofmcsteps)

            self.centeratorigin() # it's nice to be in the center!
            output += [self.getradiusofgyration(),]

        mystr = ''
        for item in output:
            mystr += '%.10f\n'%item

        f = open('/tmp/deletemenow.txt','w')
        f.write(mystr)
        f.close()

        from subprocess import call
        call(["autocorr","/tmp/deletemenow.txt"])

        return output

    def setrandomseedtoclocktime(self):
        wlcmclib.setrandomseedtoclocktime()


                
def vdot(v1,v2):
    return v1[:,0]*v2[:,0]+v1[:,1]*v2[:,1]+v1[:,2]*v2[:,2]
def vnorm(input):
    return np.sqrt((input**2).sum(1))
def vrodrot(vectors, referencepoint, axis, theta):
    vectors = vectors - referencepoint
    costh = np.cos(theta)
    sinth = np.sin(theta)
    for i in range(len(vectors)):
        vectors[i] = vectors[i]*costh+np.cross(axis,vectors[i])*sinth+axis*(np.dot(axis,vectors[i]))*(1-costh)

    vectors = vectors + referencepoint
    return vectors

#  import wlcmc; self = wlcmc.wormlikechain()

def generateknotprobtestdata(numberofkuhnlengths=16,segmentsperkuhnlength=5,samples=10000,filenameformat="output/knotprob_n%d.3dc"):
    wlc = wormlikechain(numberofkuhnlengths=numberofkuhnlengths,segmentsperkuhnlength=segmentsperkuhnlength,diameter=0)
    filename = filenameformat%numberofkuhnlengths
    wlc.domontecarlosteps(100000) #warmup
    for i in range(samples):
        print ("Writing sample %d/%d to %s"%(i,samples,filename))
        print (len(wlc.vertices),wlc.bendingrigidity, wlc.diameter)
        wlc.domontecarlosteps(1000)
        wlc.centeratorigin()
        wlc.writetofilein3dcoordsformat(filename)
    return wlc

def processknotprobtestdata(sequencefile,relaxamount=10,firstnmany=None):
    import os

    gaussfile = os.path.splitext(sequencefile)[0]+".egc"
    reducedgaussfile = os.path.splitext(sequencefile)[0]+".egcr"
    homflyfile = os.path.splitext(sequencefile)[0]+".homfly"
    knotdistfile = os.path.splitext(sequencefile)[0]+".knotdist"
    tmpfilename = "/tmp/wlcmcforknotplot"

    if firstnmany == None:
        forknotplot = """gauss open %s\;seq open %s\; until seqend \\"ago %d\; gauss\; seq +\\"\;seq close"""%(gaussfile,sequencefile,relaxamount)
    # else:
    #     forknotplot = """gauss open %s\;seq open %s\; alias heya \\"ago %d\; gauss\; seq +\\"\;#%d \\"heya\\"\;\nseq close"""%(gaussfile,sequencefile,relaxamount,firstnmany)
    # # print(forknotplot)
    # os.system("echo %s > %s"%(forknotplot,tmpfilename))
    # os.system("knotplot -nog -stdin < %s"%(tmpfilename))
    os.system("echo %s | knotplot -nog -stdin"%forknotplot)
    os.system("xinger +r -znb %s > %s"%(gaussfile,reducedgaussfile))
    os.system("homfly %s > %s"%(reducedgaussfile,homflyfile))
    os.system("idknot -P %s > %s"%(homflyfile,knotdistfile))


def forvictor(diameter = 0, segmentsperkuhnlength = 10):
    self = wormlikechain(numberofsegments=int(8.1*segmentsperkuhnlength),
                        segmentsperkuhnlength=segmentsperkuhnlength,
                        kuhnlength=100,
                        diameter=diameter,)
    self.testforautocorr(numberofmcsteps=35000,numberofsamples=10) #show that autocorr is .5

    os.system("mkdir output")
    os.chdir("output")
    
    directory = "wlc_%d-segments_%d-segmentsperkuhnlength_%f-diameter"%(int(8.1*segmentsperkuhnlength),segmentsperkuhnlength,diameter)

    os.system("mkdir %s"%(directory))
    os.chdir(directory)
    
    print ("Making directory: %s"%(directory))

    # for i in range(70,80):
    for i in range(1000000):
        print ('Writing %d'%i)
        self.domontecarlosteps(15000)
        self.centeratorigin()
        # self.writetofilein3dcoordsformat(filename='%08d.txt'%j)
        np.savetxt('%08d.txt'%i,self.vertices)

if __name__ == '__main__':
    forvictor()

    # self = wormlikechain(numberofsegments=81,
    #                     segmentsperkuhnlength=10,
    #                     kuhnlength=100,
    #                     diameter=5,)
    # # self.bendingrigidity = 4
    # # print (self.bendingrigidity)
    # # processknotprobtestdata("output/knotprob_n60_10000.3dc",firstnmany=2)
    # self.testforautocorr(numberofmcsteps=35000,numberofsamples=100) #show that autocorr is .5


    # for i in range(200):
    #     print ('Writing %d of 100'%i)
    #     self.domontecarlosteps(15000)
    #     self.centeratorigin()
    #     self.writetofilein3dcoordsformat(filename='output/3dcoords.txt')
    
    # self.showmemayavi()
    pass
    #yo


#1.89891188677
#ratio: 357:643
#domontecarlosteps took 1721.683 ms
#43.3809542302