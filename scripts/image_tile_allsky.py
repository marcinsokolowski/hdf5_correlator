from __future__ import print_function
from time import gmtime, strftime
import optparse

# create command line arguments
parser=optparse.OptionParser()
parser.set_usage("Usage: casapy [--nologger] -c imaging_script.py [options] <filename>.ms")

casa_index=sys.argv.index('-c')+2  # remove CASA options
(options,args)=parser.parse_args(sys.argv[casa_index:])
vis=args[-1]  # path to calibrated CASA measurement set

#uvdist in m
uvranges=[]
uvdist_low=[1,85,35,1]
uvdist_high=[1e4,240,1e4,500]
for i in range(0,len(uvdist_low)):
    uvranges.append('%s~%s'%(uvdist_low[i],uvdist_high[i]))
uvranges=[''] #ovveride previous by including all

# beamsize ~200 arcmin
imsize=168#1440 for 5 arcmin image
cellside='50arcmin'
#imsize=1400*5#1440 for 5 arcmin image
#cellside='1arcmin' #'0.75arcmin'
#weights=['briggs','natural','uniform']
weights=['briggs']
robusts=[0.5]
stokeses=['XX', 'YY']
#stokeses=['YY']
#antennas=['Tile079', 'Tile076','']
antennas=['']
cyclefactor=30#10 #Was 1.5
threshold='600mJy' #Was 30mJy
psfmodes=['clark'] #Was also hogbom 
niter=0 # was 40000

weight='briggs'
for weight in weights:
    print('Imaging %s'%(vis))
    name_base=vis 
    for antenna in antennas:
        for stokes in stokeses:
            for robust in robusts:
                for psfmode in psfmodes:
                    for uvrange in uvranges:
                        if antenna!='':
                            cyclefactor=cyclefactor
                        else:
                            cyclefactor=10#1.5 #For full MWA have lower cyclefactor
                            threshold='30mJy'
                        print('Cleaning %s stokes %s with weight %s (robust=%s), threshold %s,  cyclefactor %s, psfmode %s' % (antenna, stokes, weight, robust, threshold, cyclefactor, psfmode))
                        print('Time now is %s' % strftime("%Y-%m-%d %H:%M:%S"))

                        imagename=name_base+'_'+weight+str(robust)+'_TH'+threshold+'_CF'+str(cyclefactor)\
                                +'_'+psfmode+'_'+antenna+'_'+stokes+'_'+cellside+'_'+str(imsize)+'px'+'_UV'+uvrange
                        print('Imagename: %s' % imagename)
                        clean (vis=vis, imagename=imagename, gridmode ='widefield', psfmode=psfmode,
                               robust=robust, weighting=weight, imagermode ='csclean', wprojplanes =1, 
                               facets =1, niter =niter, imsize =[imsize,imsize], cell =[cellside,cellside], 
                               threshold=threshold, stokes=stokes, mode ='mfs', selectdata =True, 
                               uvrange=uvrange, antenna=antenna, cyclefactor=cyclefactor,usescratch=False) # 2019-10-31 

                        exportfits(imagename+'.image', imagename+'.fits')

print('Imaging completed at %s' % strftime("%Y-%m-%d %H:%M:%S"))
