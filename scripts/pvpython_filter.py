import numpy

alpha_0 = inputs[0].PointData['alpha_0'].flatten()
alpha_1 = inputs[0].PointData['alpha_0'].flatten()

bg_field = 1
pdi = self.GetInputDataObject(0,0)
pdo = self.GetOutputDataObject(0)
numPoints = pdi.GetNumberOfPoints()
import time

ind=self.GetInput()
iin=ind.GetInformation()
timestep = time.clock()

dx = 100.0
dy = 100.0

for (j,alpha) in enumerate(alpha_0):
    x1,y1,z1 = pdi.GetPoint(j)[:3]
    x2,y2,z2 = pdi.GetPoint(j+1)[:3]
    if abs(x2-x1) <= dx and abs(x2-x1)>0:
        dx = x2 - x1
    if abs(y2-y1) <= dy and abs(y2-y1)>1e-6:
        dy = y2 - y1
    #pdb.set_trace()

    #print "y1-y2 = ", dy
    #print "x1-x2 = ", dx
counter = 0
last_step = -1.0
prev_x = 0.0
efold = numpy.zeros(numPoints)
y_first = -1.0

for (j,alpha) in enumerate(alpha_0):
    x,y,z = pdi.GetPoint(j)
    #print "y = ", y
    #if x==prex_x+dx or x==0.0:
    if x<100.0 and x>0.0 and y==0.0:
        #pdb.set_trace()
        y_first = []
        for (i,alpha1) in enumerate(alpha_1):
            x1,y1,z1 = pdi.GetPoint(i)
            if x1 == x:
                if y1 not in y_first:
                    efold[j] += exp(alpha1)*dy
                    #print "alpha1=", dy
                    y_first.append(y1)
                    #print x1,y1,x,y,efold[j],j
        #break
        if x!=last_step:
            counter = counter + 1
            print "y_integrated_density(",floor(timestep),",",counter,") = ", efold[j]/62.83, ";"
            print "x(",counter,") = ", x, ";"
            last_step = x
