# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:50:04 2018

@author: emilio
"""
def rgc(gal):

    from astropy.coordinates import ICRS, Distance, Angle, SkyCoord
    from astropy import units as u
    import numpy as np
    from bin_functions import read_gc, read_pne
    
    def correct_rgc(coord, glx_ctr=ICRS('00h42m44.33s +41d16m07.5s'),
            glx_PA=Angle('37d42m54s'),
            glx_incl=Angle('77.5d'),
            glx_dist=Distance(783, unit=u.kpc)):
                        
        """Computes deprojected galactocentric distance.
        Inspired by: http://idl-moustakas.googlecode.com/svn-history/
            r560/trunk/impro/hiiregions/im_hiiregion_deproject.pro
        Parameters
        ----------
        coord : :class:`astropy.coordinates.ICRS`
            Coordinate of points to compute galactocentric distance for.
            Can be either a single coordinate, or array of coordinates.
        glx_ctr : :class:`astropy.coordinates.ICRS`
            Galaxy center.
        glx_PA : :class:`astropy.coordinates.Angle`
            Position angle of galaxy disk.
        glx_incl : :class:`astropy.coordinates.Angle`
            Inclination angle of the galaxy disk.
        glx_dist : :class:`astropy.coordinates.Distance`
            Distance to galaxy.
        Returns
        -------
        obj_dist : class:`astropy.coordinates.Distance`
            Galactocentric distance(s) for coordinate point(s).
        """
        # distance from coord to glx centre
        sky_radius = glx_ctr.separation(coord)
        avg_dec = 0.5 * (glx_ctr.dec + coord.dec).radian
        x = (glx_ctr.ra - coord.ra) * np.cos(avg_dec)
        y = glx_ctr.dec - coord.dec
        if x==Angle('0d00m00s'):
            x=Angle('0d00m0.0001s')
        # azimuthal angle from coord to glx  -- not completely happy with this
        phi = glx_PA - Angle('90d') \
                + Angle(np.arctan(y.arcsec / x.arcsec), unit=u.rad)
       
        # convert to coordinates in rotated frame, where y-axis is galaxy major
        # ax; have to convert to arcmin b/c can't do sqrt(x^2+y^2) when x and y
        # are angles
        xp = (sky_radius * np.cos(phi.radian)).arcmin
        yp = (sky_radius * np.sin(phi.radian)).arcmin
    
        # de-project
        ypp = yp / np.cos(glx_incl.radian)
        obj_radius = np.sqrt(xp ** 2 + ypp ** 2)  # in arcmin
        obj_dist = Distance(Angle(obj_radius, unit=u.arcmin).radian * glx_dist, unit=glx_dist.unit)
    
        # Computing PA in disk (unused)
    #    obj_phi = Angle(np.arctan(ypp / xp), unit=u.rad)
        # TODO Zero out very small angles, i.e.
        # if np.abs(Angle(xp, unit=u.arcmin)) < Angle(1e-5, unit=u.rad):
        #     obj_phi = Angle(0.0)
    
        return xp, obj_dist
        
    def get_r(RA, DEC):
        R = np.zeros(len(RA))
        ys = np.zeros(len(RA))
        for i in range(0, len(RA)):
            aux = SkyCoord(RA[i], DEC[i])
            aux = aux.icrs
            ra = str(int(aux.ra.hms.h))+'h'+str(int(aux.ra.hms.m))+'m'+str(int(aux.ra.hms.s))+'s'
            dec = str(int(aux.dec.dms.d))+'d'+str(abs(int(aux.dec.dms.m)))+'m'+str(abs(int(aux.dec.dms.s)))+'s'
            coordinates = ra+' '+dec
            coord = ICRS(coordinates)    
            ygc, r = correct_rgc(coord, galcenter, Angle(str(inps[3])), Angle(str(inps[2])), Distance(inps[4], unit=u.kpc))  
            R[i]=r.value
            ys[i]=ygc
        return R, ys    
        
    RAgc, DECgc, inps = read_gc(gal)
    RApne, DECpne = read_pne(gal)
    
    pa = inps[3]
    inc = inps[2]
    d = inps[4].to(u.kpc)
    
    pa = str(pa.value)+'d'
    inc = str(inc.value)+'d'
    
    pa = Angle(pa)
    inc = Angle(inc)
    d = Distance(d.value, unit=u.kpc)
    
    RAg, DECg = inps[0], inps[1]
    galcenter = RAg+' '+DECg
    galcenter = ICRS(galcenter)

    Rgc, ysgc = get_r(RAgc, DECgc)
    Rpne, ysp = get_r(RApne, DECpne)
    
    return Rgc, Rpne, ysgc, ysp