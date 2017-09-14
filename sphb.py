import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.mlab import griddata
import bisect
import sys

if sys.version[0] == '2':
    sys.path.append('./tracer/')
    import tracer
else:
    from tracer import tracer

class SphB:
      """General class describing a magnetic field in spherical shell.
      The field is on a uniform spherical grid. The idea is that these 
      routines can be used for any magnetic field on a spherical grid,
      e.g. analytical, potential field, or magneto-frictional simulation.
      ---
      modified 28/3/14 --ary
      added plot_PlaneOfSky 02/09/14 SJE
      """

      def __init__(self, r, th, ph):
            self.nr = 0
            self.nth = 0
            self.nph = 0
            self.r0 = 0
            self.r1 = 0
            self.th0 = 0
            self.th1 = 0
            self.ph0 = 0
            self.ph1 = 0
            self.r = []
            self.th = []
            self.ph = []
            self.brg = []
            self.bthg = []
            self.bphg = []
            self.bbg = []
            self.jrg = []
            self.jthg = []
            self.jphg = []
            self.jjg = []
            self.arg = []
            self.athg = []
            self.aphg = []
            self.ssReverse = []
            self.flipIndex = []
            self.flipIndexTop=[]
            # Note: need to use a particular subclass to fill relevant values
            # Note: need to run ssInfo to generate ssReverse and flipIndex

      def ssInfo(self):
            rr=self.r

            t1=self.th

            p1=self.ph
            pp,tt=np.meshgrid(p1,t1)
            print("shape=",tt.shape, pp.shape)
            ss1=self.brg[self.nr-1,:,:]

            ssReverse=np.absolute(ss1)
            self.ssReverse=ssReverse
            self.flipIndex=np.where(ss1 < 0)
            return ss1

      def br(self, r, th, ph):
            """Return Br by interpolation of uniform grid."""
            sh = r.shape
            r = r.flatten()
            br = r*0
            th = th.flatten()
            ph = ph.flatten()
            for i in range(0,len(r)):
                  try:
                        i_r, dr = self.index_frac(self.r, r[i])
                  except IndexError:
                        continue
                  i_th, dth = self.index_frac(self.th, th[i])
                  i_ph, dph = self.index_frac(self.ph, ph[i])
                  br[i] = self.interpgrid(self.brg,i_r,dr,i_th,dth,i_ph,dph)
            br = np.reshape(br, sh)
            return br

      def bth(self, r, th, ph):
            """Return Bth by interpolation of uniform grid."""
            sh = r.shape
            r = r.flatten()
            bth = r*0
            th = th.flatten()
            ph = ph.flatten()
            for i in range(0,len(r)):
                  try:
                        i_r, dr = self.index_frac(self.r, r[i])
                  except IndexError:
                        continue
                  i_th, dth = self.index_frac(self.th, th[i])
                  i_ph, dph = self.index_frac(self.ph, ph[i])
                  bth[i] = self.interpgrid(self.bthg,i_r,dr,i_th,dth,i_ph,dph)
            bth = np.reshape(bth, sh)
            return bth

      def bph(self, r, th, ph):
            """Return Bph by interpolation of uniform grid."""
            sh = r.shape
            r = r.flatten()
            bph = r*0
            th = th.flatten()
            ph = ph.flatten()
            for i in range(0,len(r)):
                  try:
                        i_r, dr = self.index_frac(self.r, r[i])
                  except IndexError:
                        continue
                  i_th, dth = self.index_frac(self.th, th[i])
                  i_ph, dph = self.index_frac(self.ph, ph[i])
                  bph[i] = self.interpgrid(self.bphg,i_r,dr,i_th,dth,i_ph,dph)
            bph = np.reshape(bph, sh)
            return bph


      def plot2d(self, r0, cmpt):
            """Plot on spherical surface at fixed r.
            ---
            modified 28/3/14 --ary
            """

            # Initialise grid in phi and theta:
            th1 = self.th[1:-1]
            ph2, th2 = np.meshgrid(self.ph[1:-1], th1[::-1])

            # Ensure that cmpt is valid and get array of required component:
            try:
                  func = getattr(self, cmpt)
                  print(func)
            except:
                  print('Error: syntax is b.plot2d(r0, cmpt) where cmpt is ')
                  print('either br, bth or bph.')
                  sys.exit(1)
            bc = func(th2*0 + r0, th2, ph2)
            
            # Plot 2d map:
            plt.figure()
            plt.imshow(bc, interpolation='nearest', \
                       cmap=cm.jet, aspect='auto', origin='lower', \
                       extent=[self.ph0*180/np.pi, self.ph1*180/np.pi,
                               self.th1*180/np.pi, self.th0*180/np.pi])
            plt.colorbar()
            plt.title(cmpt)
            plt.xlabel('Longitude')
            plt.ylabel('Colatitude')
            plt.show()

            return ph2,th2,bc


            #################Not Working yet####################################
      def plot_PlaneOfSky(self,phi0,cmpt,min=-1,max=1):
            
            r1=self.r[1:-1]
            r2,th2=np.meshgrid(r1,self.th[1:-1])
            x2=np.multiply(r2,np.sin(th2))
            y2=np.multiply(r2,np.cos(th2))
           
            

            numcols,numrows=300,300
            
             # Ensure that cmpt is valid and get array of required component:
            try:
                  func = getattr(self, cmpt)
                  print(func)
            except:
                  print('Error: syntax is b.plot2d(r0, cmpt) where cmpt is ')
                  print('either br, bth or bph.')
                  sys.exit(1)
                  
            print(np.amin(r2), np.amax(r2), np.amin(th2),np.amax(th2), np.amin(r2*0+phi0))
            bc = func(r2, th2, r2*0+phi0)
            bc2=func(r2,th2,r2*0+((phi0+np.pi) % (np.pi*2)))
            if cmpt=='bph':
                  bc2=-bc2
            
            xxi=np.linspace(0.1,self.r1,numcols)
            yyi=np.linspace(-self.r1,self.r1,numrows)
            
            x3=np.ravel(x2)
            y3=np.ravel(y2)
            bc3=np.ravel(bc)
            bc23=np.ravel(bc2)
            xxi,yyi=np.meshgrid(xxi,yyi)
            x,y,z=x3,y3,bc3
            bci=griddata(x,y,z,xxi,yyi)
            bci2=griddata(-x,y,bc23,-xxi,yyi)

           
            im=plt.contourf(xxi,yyi,bci,np.linspace(min,max,200))
            im2=plt.contourf(-xxi,yyi,bci2,np.linspace(min,max,200))
            limbx,limby=np.sin(np.linspace(0,np.pi,100)),np.cos(np.linspace(0,np.pi,100))
            #ax.scatter(x2,y2,s=10)
            #ax.scatter(limbx,limby,c='white',s=10)
            #plt.scatter(limbx,limby,c='white',s=10)
            #fig.colorbar(im)
            #plt.axis('equal')
           # plt.show()

            return r2,th2,bc
            ###############################################
            

      def fieldlines(self, r0, th0, ph0, cubic):

            fl = []
            nmax = 1000
            maxError = 0.003**2
            minB = 1e-4
            if (cubic):
                  r1, th1, ph1 = tracer.fieldline_cubic(self.brg, self.bthg, self.bphg, \
                                                  self.r, self.th, self.ph, r0, th0, ph0, \
                                                  nmax, maxError, minB)
            else:
                  r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                                                  self.r, self.th, self.ph, r0, th0, ph0, \
                                                  nmax, maxError, minB)
            # Pack fieldlines into same structure as before (and strip
            # blank entries):
            for i in range(0, len(r0)):
                  th2 = th1[i,:]
                  ph2 = ph1[i,:]
                  r2 = r1[i,:]
                  th2 = th2[r2 > 0.0]
                  ph2 = ph2[r2 > 0.0]
                  r2 = r2[r2 > 0.0]
                  x2 = r2*np.sin(th2)*np.cos(ph2)
                  y2 = r2*np.sin(th2)*np.sin(ph2)
                  z2 = r2*np.cos(th2)
                  fl.append(np.array([x2,y2,z2]))
            return fl


      def fieldlines_sph(self, r0, th0, ph0):

            fl = []
            nmax = 1000
            maxError = 0.003**2
            minB = 1e-4
            r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                        self.r, self.th, self.ph, r0, th0, ph0, \
                        nmax, maxError, minB)
            # Pack fieldlines into same structure as before (and strip
            # blank entries):
            for i in range(0, len(r0)):
                  th2 = th1[i,:]
                  ph2 = ph1[i,:]
                  r2 = r1[i,:]
                  th2 = th2[r2 > 0.0]
                  ph2 = ph2[r2 > 0.0]
                  r2 = r2[r2 > 0.0]
                  fl.append(np.array([r2,th2,ph2]))
            return fl

      def fieldlines_jsq(self, r0, th0, ph0):
        # 
        fl = []
        nmax = 1000
        maxError = 0.003**2
        minB = 1e-4
        r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                        self.r, self.th, self.ph, r0, th0, ph0, \
                        nmax, maxError, minB)
        jsq = tracer.jsq(self.jjg**2, r1, th1, ph1, self.r, self.th, self.ph)
        # Pack fieldlines into same structure as before (and strip
        # blank entries):
        for i in range(0, len(r0)):
              th2 = th1[i,:]
              ph2 = ph1[i,:]
              r2 = r1[i,:]
              th2 = th2[r2 > 0.0]
              ph2 = ph2[r2 > 0.0]
              r2 = r2[r2 > 0.0]
              #x2 = r2*np.sin(th2)*np.cos(ph2)
              #y2 = r2*np.sin(th2)*np.sin(ph2)
              #z2 = r2*np.cos(th2)
              fl.append(np.array([r2,th2,ph2]))
        return jsq, fl

      def fieldlines_jpar(self, r0, th0, ph0):
        # 
        fl = []
        nmax = 1000
        maxError = 0.003**2
        minB = 1e-4
        r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                        self.r, self.th, self.ph, r0, th0, ph0, \
                        nmax, maxError, minB)
	## Here, we can use the existing routine to average |J|^2 along
	##   the field line to instead average the parallel current.
	##   Just pass along a computed array of the parallel current.
        jparg = ((self.brg * self.jrg) + (self.bthg * self.jthg) + (self.bphg * self.jphg)) / (self.bbg)**2
        jpar = tracer.njsq(jparg, r1, th1, ph1, self.r, self.th, self.ph)
        # Pack fieldlines into same structure as before (and strip
        # blank entries):
        for i in range(0, len(r0)):
              th2 = th1[i,:]
              ph2 = ph1[i,:]
              r2 = r1[i,:]
              th2 = th2[r2 > 0.0]
              ph2 = ph2[r2 > 0.0]
              r2 = r2[r2 > 0.0]
              fl.append(np.array([r2,th2,ph2]))
        return jpar, fl

#      def fieldlinespy(self, r0, th0, ph0):
#
#            fl = []
#            for i in range(0, len(r0)):
#                  print i, 'of', len(r0)
#                  r1, th1, ph1 = tracerpy.fieldline(self.brg, self.bthg, \
#                              self.bphg, self.r, self.th, self.ph, \
#                              r0[i], th0[i], ph0[i])
#                  x1 = r1*np.sin(th1)*np.cos(ph1)
#                  y1 = r1*np.sin(th1)*np.sin(ph1)
#                  z1 = r1*np.cos(th1)
#                  fl.append(np.array([x1,y1,z1]))
#
#            return fl

      def fieldlines_hlcy(self, r0, th0, ph0):
        # 
        fl = []
        nmax = 1000
        maxError = 0.003**2
        minB = 1e-4
        r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                        self.r, self.th, self.ph, r0, th0, ph0, \
                        nmax, maxError, minB)
	## Here, we can compute the field-line helicity by integrating
	##   A dot B / |B| along each fieldline
        hlarg = ((self.arg * self.brg) + (self.athg * self.bthg) + (self.aphg * self.bphg)) / (self.bbg)
        hlcy = tracer.njsq(hlarg, r1, th1, ph1, self.r, self.th, self.ph)
        # Pack fieldlines into same structure as before (and strip
        # blank entries):
        for i in range(0, len(r0)):
              th2 = th1[i,:]
              ph2 = ph1[i,:]
              r2 = r1[i,:]
              th2 = th2[r2 > 0.0]
              ph2 = ph2[r2 > 0.0]
              r2 = r2[r2 > 0.0]
              fl.append(np.array([r2,th2,ph2]))
        return hlcy, fl

      def fieldlines_hlcy_jpar(self, r0, th0, ph0):
        # 
        fl = []
        nmax = 1000
        maxError = 0.003**2
        minB = 1e-4
        r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                        self.r, self.th, self.ph, r0, th0, ph0, \
                        nmax, maxError, minB)
	## Here, we can compute the field-line helicity by integrating
	##   A dot B / |B| along each fieldline
        hlarg = ((self.arg * self.brg) + (self.athg * self.bthg) + (self.aphg * self.bphg)) / (self.bbg)
        hlcy = tracer.njsq(hlarg, r1, th1, ph1, self.r, self.th, self.ph)
        jparg = ((self.brg * self.jrg) + (self.bthg * self.jthg) + (self.bphg * self.jphg)) / (self.bbg)**2
        jpar = tracer.njsq(jparg, r1, th1, ph1, self.r, self.th, self.ph)
        # Pack fieldlines into same structure as before (and strip
        # blank entries):
        for i in range(0, len(r0)):
              th2 = th1[i,:]
              ph2 = ph1[i,:]
              r2 = r1[i,:]
              th2 = th2[r2 > 0.0]
              ph2 = ph2[r2 > 0.0]
              r2 = r2[r2 > 0.0]
              fl.append(np.array([r2,th2,ph2]))
        return hlcy, jpar, fl

      def fieldlines_hlcy_jpar_len(self, r0, th0, ph0):
        # 
        fl = []
        nmax = 1000
        maxError = 0.003**2
        minB = 1e-4
        r1, th1, ph1 = tracer.fieldline(self.brg, self.bthg, self.bphg, \
                        self.r, self.th, self.ph, r0, th0, ph0, \
                        nmax, maxError, minB)
	## Here, we can compute the field-line helicity by integrating
	##   A dot B / |B| along each fieldline
        hlarg = ((self.arg * self.brg) + (self.athg * self.bthg) + (self.aphg * self.bphg)) / (self.bbg)
        hlcy = tracer.jsq(hlarg, r1, th1, ph1, self.r, self.th, self.ph)
        jparg = ((self.brg * self.jrg) + (self.bthg * self.jthg) + (self.bphg * self.jphg)) / (self.bbg)**2
        jpar = tracer.jsq(jparg, r1, th1, ph1, self.r, self.th, self.ph)
        # Pack fieldlines into same structure as before (and strip
        # blank entries):
        for i in range(0, len(r0)):
              th2 = th1[i,:]
              ph2 = ph1[i,:]
              r2 = r1[i,:]
              th2 = th2[r2 > 0.0]
              ph2 = ph2[r2 > 0.0]
              r2 = r2[r2 > 0.0]
              fl.append(np.array([r2,th2,ph2]))
        return hlcy, jpar, fl

      def index_frac(self, x, x0):
            """Get cell index and fraction for interpolation of x0 in
            regular array x.
            """
            index = bisect.bisect(x, x0) - 1
            if index < 0: raise IndexError
            if index > len(x)-2: raise IndexError
            frac = (x0 - x[index]) / (x[index+1] - x[index])
            return index, frac

      def interpgrid(self, a, i_r, dr, i_th, dth, i_ph, dph):
            """Fast 3D, linear interpolation on an integer grid.
            
            i_r, i_th, i_ph are integer indices of array a corresponding 
            to the cell of the array in which the point lies.
            
            dr,dth,dph are the fractional indices corresponding to the location
            within the cell. These are typically between 0 and 1 (not
            required.)
            """      
            a000 = a[i_r, i_th, i_ph]
            a001 = a[i_r, i_th, i_ph+1]
            a010 = a[i_r, i_th+1, i_ph]
            a100 = a[i_r+1, i_th, i_ph]
            a011 = a[i_r, i_th+1, i_ph+1]
            a101 = a[i_r+1, i_th, i_ph+1]
            a110 = a[i_r+1, i_th+1, i_ph]
            a111 = a[i_r+1, i_th+1, i_ph+1]
            ai = (1-dth)*(1-dph)*((1-dr)*a000 + dr*a100) \
                 + (1-dth)*dph*((1-dr)*a001 + dr*a101) \
                 + dth*(1-dph)*((1-dr)*a010 + dr*a110) \
                 + dth*dph*((1-dr)*a011 + dr*a111)
            return ai

