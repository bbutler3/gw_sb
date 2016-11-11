import numpy as np
from astropy.io import fits
import os.path
import timeit
import sys
from glob import glob
start_time = timeit.default_timer()
elapsed0 = 0
elapsed00 = 0 

#INPUT
expnums = [475906, 475907]  #or iterate through range(expnum1,expnum2)
#nites = [20150302, 20150917, 20150920, 20150926, 20151007, 20151227]

#final = '/data/des41.b/data/rbutler/final_tables/475909_z_final.fits'
#if os.path.isfile(final):
#    sys.exit("ERROR: There's an existing file with the name designated at the " 
#             "end of the code. Either delete that or rename the output file.")

for exp in expnums:
    last = 8 #CCD
    first = 6 #CCD
    
    f = 0
    g = 0
    
    band = []
    truemag = []
    trueflux = []
    obsmag = []
    obsmagerr = []
    nite = []
    chips = []
    expnum = []
    CCD = []
    mjd = []
    season = []
    
    _id = []
    ra = []
    dec = []
    mag = []
    magerr = []
    ccdnum = []
    fakeid = []
    
    fakelist = []
    
    sb = []
    sb_err = []
    std = []
    rad = []
    err = []                                                                              
    for index in range(last-first+1):
        f = f+1
        #chip = 'z_'
        number = str("%02d" % (first + f - 1,))
        chip = number
        part1 = '/pnfs/des/persistent/gw/exp/*'
        part2 = '/dp107/*'
        part3 = '/*+fakeSN_filterObj.out'
        filename = part1 + '/' + str(exp) + part2 + chip + part3
        #print "23:  ", timeit.default_timer() - start_time
	#print glob(filename)
	globfil = glob(filename)
	#print globfil
        if len(globfil)>0 and os.path.isfile(globfil[0]):
            #print "25:  ", timeit.default_timer() - start_time
            g = g+1
            lineglobfil = open(globfil[0])
            line1 = lineglobfil.readline()
            line2 = lineglobfil.readline()
            
            ID,MJD2,CCDNUM2,reject,RA,DEC,MAG,MAGERR,x,y,SN_FAKEID = np.genfromtxt(\
                globfil[0], delimiter=' ',skip_header=18, skip_footer=5,\
                usecols=(1,3,5,6,8,9,10,11,12,13,45),unpack=True)
            
            x = np.around(x)
            y = np.around(y)
            
            #chip2 = 'z_'
            number2 = str("%02d" % (first + f - 1,))
            chip2 = number2
            
            part1a = '/pnfs/des/persistent/gw/exp/*'
            part2a = '/dp107/*'
            part3a = '/*template.fits'
            image = part1a + '/' + str(exp) + part2a + chip + part3a
            imagefil = glob(image)[0]
                
            fits1 = fits.open(imagefil)
            
            data1 = fits1[0].data
            
            scale = 0.2623 # ds9 ; header says .27 # arcsec/pixel in image
            
            template_zp = 31.1982 
            filterObj_zp = 31.4
            
            avg_bg = 0.
            
            a = 0
            b = 0
            c = 0 
            d = 0
            #print "79:  ", timeit.default_timer() - start_time
            for p in range(len(ID)):
                if reject[p]==0.:
                    b=b+1
            
            print ' '
            print 'total candidates =',len(ID), "for chip", number2
            print 'non-rejects =',b
            
            for q in range(len(ID)):
                if SN_FAKEID[q]!=0.:
                    c=c+1
            
            print 'fakes =',c
    #        print "88:  ", timeit.default_timer() - start_time
            for i in range(len(ID)):
		    #print "h1"
                if reject[i]==0:
                    #print "h2"
		    if SN_FAKEID[i]==0:
   	                #print "h3"
    #                    print "96:  ", timeit.default_timer() - start_time
                        # image coordinates
                        x1=0
                        y1=0
                        x_low, x_up, y_low, y_up = 0,0,0,0
                        
                        x_low, y_up = x[i]-20, y[i]+20 
                        x_up, y_low = x[i]+20, y[i]-20 
                        
                        ###
                        
                        #coordinate adjustment for python starting at 0
                        x_low, x_up, y_low, y_up = x_low-1, x_up-1,\
                            y_low-1, y_up-1
                        
                        x1 = x[i]-1
                        y1 = y[i]-1
                        
                        ###
                        
                        ### be aware: in these cuts, the x-coordinate comes SECOND ###
                        
                        box = data1[y_low:y_up,x_low:x_up]
                        box = box.ravel()
                        
                        # subtract out average background
                        box = box - avg_bg
                        
                        #for i in range(len(box_r1)):
                        #    if box_r1[i] <= 0:
                        #        box_r1[i] = 0.00001
                        
                        row, col= np.indices((4096,2048)) 
                        
                        new_row = row - y[i]
                        new_col = col - x[i]
                        
                        row_x = new_row[y_low:y_up,x_low:x_up]
                        col_y = new_col[y_low:y_up,x_low:x_up]
                        
                        row_x1 = row_x.ravel()
                        col_y1 = col_y.ravel()
                        
                        #image_x = row_x1 + 5374
                        #image_y = col_y1 + 2057
                        
                        radius = scale*np.sqrt(col_y1**2+row_x1**2)
                        
                        ### loop that specifies area used in calculations ###
                        
                        z_radcut = 2.    #radius of circle around center
                        #z_isocut = 1400.   #pixel intensity lower cutoff
                        
                        z_circle = []
                        error = []
    #                    print "151: ", timeit.default_timer() - start_time
                        for i2 in range(len(box)):
                            if ((radius[i2] <= z_radcut)): # and 
                            #((box[i]+avg_bg) >=z_isocut)):
                                z_circle.append(box[i2])
                        
                        #since the gain in this particular image's header is ~0.98, I just used
                        # Poisson square root error in each pixel, and then averaged those 
                        # errors later.
                        for i3 in range(len(box)):
                            if ((radius[i3] <= z_radcut)):
                                poisson = np.sqrt(abs(box[i3]))
                                error.append(poisson)
                        #print "avg", np.mean(error), "new", np.sqrt(np.sum(np.square(error)))/\
                            len(z_circle)
                        ######################################################
                
                        #rmag = template_zp -2.5*numpy.log10(box_r1/(scale**2))
                        #gmag = template_zp -2.5*numpy.log10(box_g1/(scale**2))
                        #g_r = gmag-rmag
                        
                        ### surface brightness and other calculations ###
                        
                        #length = len(z_circle)
                        #print 'circle radius ("):  ', z_radcut
                        #print "# of pixels:        ", length
                        
                        sum_circle = np.sum(z_circle)
                        med = np.median(z_circle)
                        errmean = np.mean(error)
                        err.append(errmean)
                        stdev = np.std(z_circle)
                        std.append(stdev)
                        #med_and_stdev = ufloat(med, stdev)
                        
                        #mean_circle = np.mean(z_circle)
                        #print "total flux:         ", sum_circle
                        #print "median intensity:   ", med
                        #print "mean intensity:     ", mean_circle
                        
                        #set flat background to 30 mag/arcsec^2
                        if med<=.207432394746636:
                            med = .207432394746636
                            
                #     print "med", med, "errmean", errmean, "stdev", stdev
                            
                        
                        #sb_sum_circle = template_zp -2.5*numpy.log10(sum_circle/(scale**2))
                        #print "SB sum", sb_sum_circle
                        
                        sb_med_circle = template_zp -2.5*np.log10(med/(scale**2))
                        #sb_med_circle = template_zp -2.5*unumpy.log10(med_and_stdev/(scale**2))
                        #dsb_plus = template_zp - 2.5*np.log10(1+stdev/med)
                        #dsb_minus = template_zp - 2.5*np.log10(1-stdev/med)
                        
                        sn_g1 = 1.086*(errmean/med)
                        #sn_g1 = 1.086*((np.sqrt(np.sum(np.square(error)))/len(z_circle))/med)
                        
                        #print "SB w/ median:       ", round(sb_med_circle,3)
                    
                        p = med/(scale**2)
                        sp = stdev/(scale**2)
                        errcalc = template_zp -2.5*.434*(sp/p)
                        
                        #sb_mean_circle = template_zp -2.5*numpy.log10(mean_circle/(scale**2))
                        #print "SB w/ mean:         ", round(sb_mean_circle,3)
                        #print sb_med_circle
                        sb.append(sb_med_circle)
                        sb_err.append(sn_g1)
                #     print p, sp, sb_med_circle, sn_g1 #dsb_plus, dsb_minus
                        #if we use mean instead of median, the errors go down to ~2 mag/arcsec**2
                        #using poisson errors, since gain is 0.98, we can reduce it much more
                        #(took mean of poisson errors for each pixel
                    
                        band.append(line1.split()[1])
                        nite.append(int(line2.split()[1]))
                        chips.append(number2)
                        expnum.append(exp)
                        season.append(107)
                        rad.append(2.)
                        truemag.append(0.)
                        trueflux.append(0.)
            #            meds.append(med)
    #                    print "233: ", timeit.default_timer() - start_time
            ############################################################################
            ###############################    FAKES    ################################
            ############################################################################
                    
                    else: 
    #                    print "fake?"
                        d=d+1 
                        part1b = '/data/des41.b/data/rbutler/10_19_2016/'
                        part2b = 'fake_list_dp107_sorted_uniq.tab'      
                        table1 = part1b + part2b
            
                        FAKEID,EXPNUM,CCDNUM1,TRUEMAG,TRUEFLUX,FLUXCOUNT,NITE,MJD1 = \
                            np.genfromtxt(table1, delimiter=' ',\
                            skip_header=1,usecols=(0,1,2,3,4,5,7,8),unpack=True)
            
                        BAND = np.genfromtxt(table1,dtype = 'string',delimiter=' ',
                            skip_header=1,usecols=(6,),unpack=True)
                            
                        MJD1 = np.around(MJD1,decimals=3)
                                                
                        for j in range(len(FAKEID)):
    #                        print j
    #                        if SN_FAKEID[i]==FAKEID[j]:
    #                            print "ID"
    #                            print MJD2[i], MJD1[j]
    #                            print CCDNUM2[i], CCDNUM1[j]
    #                            print BAND[j]
    #                        if CCDNUM2[i]==CCDNUM1[j]:    
    #                            print "CCD"
    #                        if (MJD2[i] <= (MJD1[j]+.002)) and \
    #                            (MJD2[i] >= (MJD1[j]-.002)):
    #                            print "MJD"
    #                        if BAND[j]=='i':
    #                            print "BAND" 
    #                        print " "
                            if (SN_FAKEID[i] == FAKEID[j]) and (CCDNUM2[i] == CCDNUM1[j]) \
                                and (MJD2[i] <= (MJD1[j]+.002)) and (MJD2[i] >= (MJD1[j]-.002)) and (BAND[j] == line1.split()[1]):
    #                                print "next loop"
                                    a = a+1
                                    if SN_FAKEID[i] not in fakelist:
                                        print int(SN_FAKEID[i])
                                        #ID.append(SN_FAKEID[i])
                                        band.append(BAND[j])
                                        truemag.append(TRUEMAG[j])
                                        trueflux.append(TRUEFLUX[j])
                                        #obsmag.append(MAG[i])
                                        #obsmagerr.append(MAGERR[i])
                                        nite.append(NITE[j])
                                        chips.append(number2)
                                        expnum.append(EXPNUM[j])
                                        #CCD.append(CCDNUM1[j])
                                        #mjd.append(MJD1[j])
                                    fakelist.append(SN_FAKEID[i])
                                    #print fakelist
                        # image coordinates
                        x1=0
                        y1=0
                        x_low, x_up, y_low, y_up = 0,0,0,0
                        
                        x_low, y_up = x[i]-20, y[i]+20 
                        x_up, y_low = x[i]+20, y[i]-20 
                        
                        #print x_low
                        # center of object/area 
                        
                        ###
                        
                        #coordinate adjustment for python starting at 0
                        x_low, x_up, y_low, y_up = x_low-1, x_up-1,\
                            y_low-1, y_up-1
                        
                        x1 = x[i]-1
                        y1 = y[i]-1
                        
                        ###
                        
                        ### be aware: in these cuts, the x-coordinate comes SECOND ###
                        
                        box = data1[y_low:y_up,x_low:x_up]
                        box = box.ravel()
                        
                        # subtract out average background
                        box = box - avg_bg
                        
                        #for i in range(len(box_r1)):
                        #    if box_r1[i] <= 0:
                        #        box_r1[i] = 0.00001
                        
                        row, col= np.indices((4096,2048)) 
                        
                        new_row = row - y[i]
                        new_col = col - x[i]
                        
                        row_x = new_row[y_low:y_up,x_low:x_up]
                        col_y = new_col[y_low:y_up,x_low:x_up]
                        
                        row_x1 = row_x.ravel()
                        col_y1 = col_y.ravel()
                        
                        #image_x = row_x1 + 5374
                        #image_y = col_y1 + 2057
                        
                        radius = scale*np.sqrt(col_y1**2+row_x1**2)
                        
                        ### loop that specifies area used in calculations ###
                        
                        radcut = 2.    #radius of circle around center
                        #isocut = 1400.   #pixel intensity lower cutoff
                        
                        circle = []
                        error = []
                        
                        for i2 in range(len(box)):
                            if ((radius[i2] <= radcut)): # and 
                            #((box[i]+avg_bg) >=isocut)):
                                circle.append(box[i2])
                        
                        #since the gain in this particular image's header is ~0.98, I just used
                        # Poisson square root error in each pixel, and then averaged those 
                        # errors later.
                        for i3 in range(len(box)):
                            if ((radius[i3] <= radcut)):
                                poisson = np.sqrt(abs(box[i3]))
                                error.append(poisson)
                        #print "avg", np.mean(error), "new", np.sqrt(np.sum(np.square(error)))/\
                            len(circle)
                        ######################################################
                
                        #rmag = template_zp -2.5*numpy.log10(box_r1/(scale**2))
                        #gmag = template_zp -2.5*numpy.log10(box_g1/(scale**2))
                        #g_r = gmag-rmag
                        
                        ### surface brightness and other calculations ###
                        
                        #length = len(circle)
                        #print 'circle radius ("):  ', radcut
                        #print "# of pixels:        ", length
                        
                        sum_circle = np.sum(circle)
                        med = np.median(circle)
                        errmean = np.mean(error)
                        err.append(errmean)
                        stdev = np.std(circle)
                        std.append(stdev)
                        season.append(107)
            
                        #med_and_stdev = ufloat(med, stdev)
                        
                        #mean_circle = np.mean(circle)
                        #print "total flux:         ", sum_circle
                        #print "median intensity:   ", med
                        #print "mean intensity:     ", mean_circle
                        
                        #set flat background to 30 mag/arcsec^2
                        if med<=.207432394746636:
                            med = .207432394746636
                            
                #     print "med", med, "errmean", errmean, "stdev", stdev
                            
                        
                        #sb_sum_circle = template_zp -2.5*numpy.log10(sum_circle/(scale**2))
                        #print "SB sum", sb_sum_circle
                        
                        sb_med_circle = template_zp -2.5*np.log10(med/(scale**2))
                        #sb_med_circle = template_zp -2.5*unumpy.log10(med_and_stdev/(scale**2))
                        #dsb_plus = template_zp - 2.5*np.log10(1+stdev/med)
                        #dsb_minus = template_zp - 2.5*np.log10(1-stdev/med)
                        
                        sn_g1 = 1.086*(errmean/med)
                        #sn_g1 = 1.086*((np.sqrt(np.sum(np.square(error)))/len(circle))/med)
                        
                        #print "SB w/ median:       ", round(sb_med_circle,3)
                    
                        p = med/(scale**2)
                        sp = stdev/(scale**2)
                        errcalc = template_zp -2.5*.434*(sp/p)
                        
                        #sb_mean_circle = template_zp -2.5*numpy.log10(mean_circle/(scale**2))
                        #print "SB w/ mean:         ", round(sb_mean_circle,3)
                        #print sb_med_circle
                        sb.append(sb_med_circle)
                        sb_err.append(sn_g1)
                #     print p, sp, sb_med_circle, sn_g1 #dsb_plus, dsb_minus
                        #if we use mean instead of median, the errors go down to ~2 mag/arcsec**2
                        #using poisson errors, since gain is 0.98, we can reduce it much more
                        #(took mean of poisson errors for each pixel
                    
                #        band.append('z')
                #        nite.append(20150917)
                #        expnum.append(475915)
                #        season.append(44)
                        rad.append(2.)
                #        meds.append(med)    
                else:
                    sb.append(0.)
                    sb_err.append(0.)
                    std.append(0.)
                    err.append(0.)
                    band.append(line1.split()[1])
                    nite.append(int(line2.split()[1]))
                    chips.append(number2)
                    expnum.append(exp)
                    season.append(107)
                    rad.append(0.)
                    truemag.append(0.)
                    trueflux.append(0.)
            #        meds.append(0.)
            _id.extend(ID)
            ra.extend(RA)
            dec.extend(DEC)
            mag.extend(MAG)
            magerr.extend(MAGERR)
            ccdnum.extend(CCDNUM2)
            fakeid.extend(SN_FAKEID)           
            fits1.close()
            print 'a =',a
            #print sb
            #print errmean
            #print sb_err
            #print sb_err[181]
            #print meds[230]
            print 'd =',d
            #print MJD2
            #print MJD1
            
            
            
            
            elapsed1 = timeit.default_timer() - start_time - elapsed0
            elapsed0 = timeit.default_timer() - start_time
            
            m, s = divmod(elapsed0, 60)
            h, m = divmod(m, 60)
            
            print 'time elapsed:',round(elapsed1,2),'seconds for chip',number2
            print '             ',"%d:%02d:%02d" % (h, m, s), 'overall'
#        print glob(filename)
#        if len(glob(filename))>0 and os.path.isfile(glob(filename)[0]): 
    mag = [alpha + 31.4 for alpha in mag]
    
    MAG_ERR_div_MAG = [100*(beta/gamma) for beta,gamma in zip(magerr,mag)]
    
    tbhdu1 = fits.BinTableHDU.from_columns(
        [fits.Column(name='ID', format='K', array=_id),
        fits.Column(name='RA', format='E', array=ra),
        fits.Column(name='DEC', format='E', array=dec),
        fits.Column(name='radius', format='E', array=rad),
        fits.Column(name='band', format='1A', array=band),
        fits.Column(name='mag', format='E', array=mag),
        fits.Column(name='mag_err', format='E', array=magerr),
        fits.Column(name='percent mag_err', format='E', array=MAG_ERR_div_MAG),
        fits.Column(name='true_mag', format='E', array=truemag),
        fits.Column(name='true_flux', format='E', array=trueflux),
        fits.Column(name='NITE', format='K', array=nite),
        fits.Column(name='EXPNUM', format='K', array=expnum),
        fits.Column(name='CCDNUM', format='K', array=ccdnum),
    #    fits.Column(name='CHIP', format='K', array=chips),
        fits.Column(name='SEASON', format='K', array=season),
        fits.Column(name='SB', format='E', array=sb),
        fits.Column(name='SB_err', format='E', array=sb_err),
        fits.Column(name='fakeID', format='K', array=fakeid)])
    
    #chip3 = 'z_'
    #number3 = str("%02d" % (first + f - 1,))
    #chip3 = chip3 + number3
    
    #newfile = '475908_710--660_'+chip3+'out.fits'
    final1 = '/data/des41.b/data/rbutler/final_tables/'
    final2 = '_final.fits'
    final = final1 + str(exp) + final2
    #tbhdu1.writeto(final)
    #print int(exp),'written'
    
    elapsed01 = timeit.default_timer() - start_time - elapsed00
    elapsed00 = timeit.default_timer() - start_time
    
    m0, s0 = divmod(elapsed01, 60)
    h0, m0 = divmod(m0, 60)
    
    print ' '
    print "%d:%02d:%02d" % (h0, m0, s0),'elapsed for exposure',int(exp)
         
elapsed = timeit.default_timer() - start_time
mt, st = divmod(elapsed, 60)
ht, mt = divmod(mt, 60)
print ' ' 
print 'time elapsed:',"%d:%02d:%02d" % (ht, mt, st), "total"
