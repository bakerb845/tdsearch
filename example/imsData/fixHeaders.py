#!/usr/bin/python
import glob
import string
from obspy.io.sac import SACTrace
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.utcdatetime import UTCDateTime

oriTime = '2016-09-09:00:30:01.440'
oTime = UTCDateTime(oriTime)
khole = '00'
sacfls = glob.glob('*sac')
for sacfl in sacfls:
   #print sacfl
   fnameSplit = string.split(sacfl,'.') 
   netw = fnameSplit[0].upper()
   stat = fnameSplit[1].upper()
   chan = fnameSplit[2].upper()
   chan = string.split(chan,'_')[0]
   #print sacfl, netw, stat, chan, khole
   sac = SACTrace.read(sacfl)
   #DPRK
   sac.evla = 41.287
   sac.evlo = 129.078
   sac.evdp = 1.0 
   sac.cmpinc = 0.0
   sac.cmpaz = 0.0
   sac.a =-12345.0
   sac.ka = '-12345\0'
   sac.t0 =-12345.0
   sac.kt0 = '-12345\0\0'
   sac.ko = 'O\0' 
   sac.knetwk = netw + '\0'
   sac.kstnm = stat + '\0' 
   sac.kcmpnm = chan + '\0'
   sac.khole = khole + '\0'
   sac.kevnm = '\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0'
   #sac.gcarc =-12345.
   sac.dist =-12345.
   sac.az =-12345.
   sac.baz =-12345.
   #print dir(sac)
   sTime = UTCDateTime(year=sac.nzyear, julday=sac.nzjday,
                       hour=sac.nzhour, minute=sac.nzmin,
                       second=sac.nzsec) + sac.nzmsec/100.0 + sac.b
   o = oTime - sTime # origin time relative to start time
   sac.o = o
   #gps2dist_azimuth(sac.evla, sac.evlo, sac.stla, sac.stlo) 
   oflnm = netw + '.' + stat + '.' + chan + '.' + khole + '.SAC'
   #print oflnm
   sac.write(oflnm)
